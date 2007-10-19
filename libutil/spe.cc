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

#include <cerrno>
#include <string>
#include <tr1/functional>

#include <libspe2.h>
#include <libwrapiter/libwrapiter_forward_iterator.hh>

using namespace honei;

SPEError::SPEError(const std::string & msg, const std::string & reason) :
    ExternalError("libspe2", msg + ", reason is " + reason + ".")
{
}

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
        signed int retval(0);

        imp->kernel->load(SPE(imp));

        do
        {
            retval = spe_context_run(imp->context, &entry_point, 0, imp->kernel->argument(),
                    imp->kernel->environment(), 0);
        }
        while (retval > 0);

        if (retval < 0)
        {
            throw ExternalError("libspe2", "spe_context_run failed, " + stringify(strerror(errno)));
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
            std::string reason;
            switch (errno)
            {
                case ENOMEM:
                    reason = "lack of system resources (ENOMEM)";
                    break;

                case EPERM:
                    reason = "insufficient permissions to create context (EPERM)";
                    break;

                case EFAULT:
                    reason = "internal error (EFAULT)";
                    break;

                default:
                    reason = "unknown error (errno = " + stringify(errno) + ")";
            }

            throw SPEError("Could not create SPE context", reason);
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
            std::string reason;
            switch (errno)
            {
                case ESRCH:
                    reason = "invalid context (ESRCH)";
                    break;

                case EFAULT:
                    reason = "internal error (EFAULT)";
                    break;

                default:
                    reason = "unknown error (errno = " + stringify(errno) + ")";
            }

            Log::instance()->message(ll_minimal, "Could not destroy SPE context " + reason + ".");
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

