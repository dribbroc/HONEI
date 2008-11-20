/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * This file is part of the Utility C++ library. LibUtil is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibUtil is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/backends/cell/interface.hh>
#include <honei/backends/cell/ppe/spe_error.hh>
#include <honei/backends/cell/ppe/spe_kernel.hh>
#include <honei/backends/cell/ppe/spe_manager.hh>
#include <honei/util/lock.hh>
#include <honei/util/log.hh>
#include <honei/util/mutex.hh>
#include <honei/util/thread.hh>

#include <string>
#include <tr1/functional>
#include <fstream>

#include <libspe2.h>
#include <pthread.h>

using namespace honei;


struct SPE::Implementation
{
    /// Our libspe context.
    spe_context_ptr_t context;

    /// Our device id.
    const DeviceId device;

    /// Our mutex.
    Mutex * const mutex;

    /// Our SPE (blocking) thread.
    pthread_t * thread;

    /// Our thread's attributes.
    pthread_attr_t * attr;

    /// Our current SPE kernel.
    SPEKernel * kernel;

    /// Return the next device id.
    inline DeviceId next_device_id()
    {
        static DeviceId result(0);

        return result++;
    }

    static inline std::string stop_reason(const spe_stop_info_t & stop_info)
    {
        std::string result("SPE: Reason: ");

        switch (stop_info.stop_reason)
        {
            case SPE_EXIT:
                result += "SPE exited with exit code '" + stringify(stop_info.result.spe_exit_code) + "'";
                break;

            case SPE_STOP_AND_SIGNAL:
                result += "SPE stopped and signaled code '" + stringify(stop_info.result.spe_signal_code) + "'";
                break;

            case SPE_RUNTIME_ERROR:
                result += "SPE runtime error";
                break;

            case SPE_RUNTIME_EXCEPTION:
                result += "SPE runtime exception due to ";
                switch (stop_info.result.spe_runtime_exception)
                {
                    case SPE_DMA_ALIGNMENT:
                        result += "DMA alignment error";
                        break;

                    case SPE_DMA_SEGMENTATION:
                        result += "DMA segmentation error";
                        break;

                    case SPE_DMA_STORAGE:
                        result += "DMA storage error";
                }
                break;

            case SPE_RUNTIME_FATAL:
                result += "SPE runtime fatal error";
                break;

            case SPE_CALLBACK_ERROR:
                result += "SPE callback error";
                break;

            default:
                result += "unknown error (" + stringify(stop_info.stop_reason) + ")";
        }

        return result;
    }

    static void * thread_function(void * argument)
    {
        SPE * spe(static_cast<SPE *>(argument));
        Implementation * imp(spe->_imp.get());

        CONTEXT("When running execution thread for SPE #"
                + stringify(imp->device) + " with kernel '"  + imp->kernel->kernel_info().name + "':");

        {
            Lock l(*imp->mutex);
            imp->kernel->load(*spe);
        }

        try
        {
            cell::LocalStoreAddress entry_point = { SPE_DEFAULT_ENTRY };
            signed int retval(0);
            spe_stop_info_t stop_info;

            do
            {
                retval = spe_context_run(imp->context, &entry_point.value, 0, imp->kernel->argument(),
                        imp->kernel->environment(), &stop_info);
            }
            while (retval > 0);

            LOGMESSAGE(ll_minimal, "SPE #" + stringify(imp->device) + " stopped at " + stringify(entry_point));
            LOGMESSAGE(ll_minimal, stop_reason(stop_info));

            if (retval < 0)
            {
                spe->dump("spu-" + stringify(spe->id()) + "-dump-final");
                throw SPEError("spe_context_run", errno);
            }
        }
        catch (Exception & e)
        {
            LOGMESSAGE(ll_minimal, e.message());
            throw;
        }

        LOGMESSAGE(ll_minimal, "SPE: Thread exiting");

        pthread_exit(0);
    }

    /// Constructor.
    Implementation() :
        mutex(new Mutex),
        context(spe_context_create(SPE_EVENTS_ENABLE, 0)),
        device(next_device_id()),
        thread(0),
        attr(0),
        kernel(0)
    {
        if (! context)
        {
            throw SPEError("spe_context_create", errno);
        }

        LOGMESSAGE(ll_minimal, "SPE Implementation (" + stringify(this) + ") created.");
    }

    /// Destructor.
    ~Implementation()
    {
        CONTEXT("When destroying SPE:");
        LOGMESSAGE(ll_minimal, "SPE: Destroying SPE::Implementation!");

        pthread_join(*thread, 0);

        // Kill the thread.
        delete thread;

        pthread_attr_destroy(attr);

        delete attr;

        // Delete mutex.
        delete mutex;

        // Release the context.
        int retval(spe_context_destroy(context)), counter(5);

        while ((retval == -1) && (errno == EAGAIN) && (counter > 0))
        {
            LOGMESSAGE(ll_minimal, "SPE '" + stringify(device)
                    + "' has running threads on destruction.");
            usleep(1000);
            --counter;
        }

        if (retval == -1)
        {
            LOGMESSAGE(ll_minimal, "spe_context_destroy failed, " + stringify(strerror(errno)));
        }
    }
};

SPE::SPE() :
    _imp(new Implementation)
{
    LOGMESSAGE(ll_minimal, "SPE(" + stringify(this) + ") created with imp pointing to " + stringify(_imp.get()));
}

SPE::SPE(const SPE & other) :
    _imp(other._imp)
{
    LOGMESSAGE(ll_minimal, "SPE(" + stringify(this) + ") created from SPE(" + stringify(&other) + "), imp = " + stringify(_imp.get()));
}

SPE::~SPE()
{
    CONTEXT("When destroying SPE:");
    LOGMESSAGE(ll_minimal, "SPE (" + stringify(this) + ") destroyed.");
}

spe_context_ptr_t
SPE::context() const
{
    return _imp->context;
}

void
SPE::dump(const std::string & filename) const
{
    CONTEXT("When dumping local store memory of SPE #" + stringify(_imp->device) + ":");
    Lock l(*_imp->mutex);

    std::fstream file(filename.c_str(), std::ios::out);
    char * ls_area(static_cast<char *>(spe_ls_area_get(_imp->context)));

    for (char * c(ls_area), * c_end(ls_area + 256 * 1024) ; c != c_end ; ++c)
    {
        file << *c;
    }

    LOGMESSAGE(ll_minimal, "SPE: Dumped LS content to file '" + filename + "'");
}

void
SPE::run(const SPEKernel & kernel)
{
    CONTEXT("When loading and running kernel in SPE(" + stringify(this) + ")");
    Lock l(*_imp->mutex);

    LOGMESSAGE(ll_minimal, "Loading Kernel '" + stringify(kernel.kernel_info().name) + "' into SPE #" + stringify(_imp->device));

    if (_imp->kernel)
        delete _imp->kernel;

    if (_imp->thread)
    {
        pthread_join(*_imp->thread, 0);
        delete _imp->thread;
        _imp->thread = 0;
    }

    if (_imp->attr)
    {
        pthread_attr_destroy(_imp->attr);
        delete _imp->attr;
        _imp->attr = 0;
    }

    _imp->kernel = new SPEKernel(kernel);

    // Release the context.
    int retval(spe_context_destroy(_imp->context)), counter(5);

    while ((retval == -1) && (errno == EAGAIN) && (counter > 0))
    {
        LOGMESSAGE(ll_minimal, "SPE '" + stringify(_imp->device)
                + "' has still running threads on destruction.");
        usleep(1000);
        --counter;
    }

    if (retval == -1)
    {
        LOGMESSAGE(ll_minimal, "spe_context_destroy failed, " + stringify(strerror(errno)));
    }

    // Create new context.
    _imp->context = spe_context_create(SPE_EVENTS_ENABLE, 0);;

    _imp->thread = new pthread_t;

    _imp->attr = new pthread_attr_t;

    if (0 != (retval = pthread_attr_init(_imp->attr)))
        throw PThreadError("pthread_attr_init", retval);

    if (0 != (retval = pthread_attr_setdetachstate(_imp->attr, PTHREAD_CREATE_JOINABLE)))
        throw PThreadError("pthread_attr_setdetachstate", retval);

    if (0 != (retval = pthread_create(_imp->thread, _imp->attr, &Implementation::thread_function, new SPE(*this))))
        throw ExternalError("libpthread", "pthread_create failed, " + stringify(strerror(retval)));

}

DeviceId
SPE::id() const
{
    Lock l(*_imp->mutex);
    return _imp->device;
}

bool
SPE::idle() const
{
    return false; /// \todo Implement bool Thread::finished()?
}

SPEKernel *
SPE::kernel() const
{
    Lock l(*_imp->mutex);
    return _imp->kernel;
}

void
SPE::send_mail(unsigned int mail) const
{
    Lock l(*_imp->mutex);

    spe_in_mbox_write(_imp->context, &mail, 1, SPE_MBOX_ANY_BLOCKING);
}

void
SPE::signal() const
{
    Lock l(*_imp->mutex);

    static unsigned signal_value(0x1234);
    spe_signal_write(_imp->context, SPE_SIG_NOTIFY_REG_1, signal_value); /// \todo remove hardcoded numbers
    ++signal_value;
}

