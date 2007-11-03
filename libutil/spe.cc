/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#include <libutil/lock.hh>
#include <libutil/log.hh>
#include <libutil/mutex.hh>
#include <libutil/spe_kernel.hh>
#include <libutil/spe_manager.hh>
#include <libutil/thread.hh>

#include <string>
#include <tr1/functional>

#include <libspe2.h>
#include <libwrapiter/libwrapiter_forward_iterator.hh>

using namespace honei;


struct SPE::Implementation
{
    /// Our libspe context.
    spe_context_ptr_t context;

    /// Our device id.
    const DeviceId device;

    /// Our mutex.
    Mutex * const mutex;

    /// Our Thread::Function object that encapsulates spe_thread_function.
    Thread::Function thread_function;

    /// Our SPE (blocking) thread.
    Thread * thread;

    /// Our current SPE kernel.
    SPEKernel * kernel;

    /// Return the next device id.
    inline DeviceId next_device_id()
    {
        static DeviceId result(0);

        return result++;
    }

    static void * spe_thread_function(void * argument)
    {
        Implementation * imp(static_cast<Implementation *>(argument));
        CONTEXT("When running execution thread for SPE #" + stringify(imp->device) + ":");

        unsigned int entry_point(SPE_DEFAULT_ENTRY);
        spe_stop_info_t stop_info;
        signed int retval(0);

        imp->kernel->load(SPE(imp));

        try
        {
            do
            {
                retval = spe_context_run(imp->context, &entry_point, 0, imp->kernel->argument(),
                        imp->kernel->environment(), &stop_info);
            }
            while (retval > 0);

            if (retval < 0)
            {
                Log::instance()->message(ll_minimal, "SPE #" + stringify(imp->device) + " stopped "
                        "at " + stringify(reinterpret_cast<void *>(entry_point)));

                std::string msg("Reason: ");
                switch (stop_info.stop_reason)
                {
                    case SPE_RUNTIME_ERROR:
                        msg += "SPE runtime error";
                        break;

                    case SPE_RUNTIME_EXCEPTION:
                        msg += "SPE runtime exception due to ";
                        switch (stop_info.result.spe_runtime_exception)
                        {
                            case SPE_DMA_ALIGNMENT:
                                msg += "DMA alignment error";
                                break;

                            case SPE_DMA_SEGMENTATION:
                                msg += "DMA segmentation error";
                                break;

                            case SPE_DMA_STORAGE:
                                msg += "DMA storage error";

                            case SPE_INVALID_DMA:
                                msg += "; Invalid DMA error reported!";
                        }
                        break;

                    case SPE_RUNTIME_FATAL:
                        msg += "SPE runtime fatal error";
                        break;

                    case SPE_CALLBACK_ERROR:
                        msg += "SPE callback error";
                        break;

                    default:
                        msg += "unknown error (" + stringify(stop_info.stop_reason) + ")";
                }
                Log::instance()->message(ll_minimal, msg);
                throw SPEError("spe_context_run", errno);
            }
        }
        catch (Exception & e)
        {
            Log::instance()->message(ll_minimal, e.message());
            throw;
        }
    }

    /// Constructor.
    Implementation() :
        context(spe_context_create(SPE_EVENTS_ENABLE, 0)),
        device(next_device_id()),
        mutex(new Mutex),
        thread_function(std::tr1::bind(spe_thread_function, this)),
        thread(0),
        kernel(0)
    {
        if (! context)
        {
            throw SPEError("spe_context_create", errno);
        }
    }

    /// Destructor.
    ~Implementation()
    {
        // Kill the thread.
        delete thread;

        // Delete mutex.
        delete mutex;

        // Release the context.
        int retval(spe_context_destroy(context)), counter(5);

        while ((retval == -1) && (errno == EAGAIN) && (counter > 0))
        {
            Log::instance()->message(ll_minimal, "SPE '" + stringify(device)
                    + "' has still running threads on destruction.");
            usleep(1000);
            --counter;
        }

        if (retval == -1)
        {
            Log::instance()->message(ll_minimal, "spe_context_destroy failed, " + stringify(strerror(errno)));
        }
    }
};

SPE::SPE() :
    _imp(new Implementation)
{
}

SPE::SPE(Implementation * imp) :
    _imp(imp)
{
}

SPE::~SPE()
{
}

spe_context_ptr_t
SPE::context() const
{
    return _imp->context;
}

void
SPE::run(const SPEKernel & kernel)
{
    CONTEXT("When loading and running kernel:");
    Lock l(*_imp->mutex);

    if (_imp->kernel)
        delete _imp->kernel;

    if (_imp->thread)
        delete _imp->thread;

    _imp->kernel = new SPEKernel(kernel);
    _imp->thread = new Thread(_imp->thread_function);
}

DeviceId
SPE::id() const
{
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

    spe_signal_write(_imp->context, SPE_SIG_NOTIFY_REG_1, 0x1234); /// \todo remove hardcoded numbers
}

