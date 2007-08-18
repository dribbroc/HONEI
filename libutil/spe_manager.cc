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
#include <libutil/spe_manager.hh>

#include <list>
#include <string>

#include <libspe2.h>

#include <libwrapiter/libwrapiter_forward_iterator.hh>

using namespace pg512;

SPEError::SPEError(const std::string & msg, const std::string & reason) :
    ExternalError("libspe2", msg + ", reason is " + reason + ".")
{
}

struct SPE::Implementation
{
    typedef std::list<SPETask *> TaskList;

    /// Our mutex.
    Mutex * const mutex;

    std::tr1::function<void () throw ()> function;

    /// Our thread.
    Thread * thread;

    /// Our queued tasks.
    TaskList task_list;

    /// Our libspe context.
    spe_context_ptr_t context;

    /// Our device id.
    const DeviceId device;

    /// Enqueue an SPETask.
    inline void enqueue(SPETask * task)
    {
        Lock l(*mutex);

        task_list.push_back(task);
    }

    /// Dequeue an SPETask and return it or return 0.
    inline SPETask * dequeue()
    {
        Lock l(*mutex);
        SPETask * result(0);

        if (! task_list.empty())
        {
            result = task_list.front();
            task_list.pop_front();
        }

        return result;
    }

    /// Return the next device id.
    inline DeviceId next_device_id()
    {
        static DeviceId result(0);

        return result++;
    }

    /// Constructor.
    Implementation() :
        mutex(new Mutex),
        context(spe_context_create(0, 0)),
        device(next_device_id())
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

        thread = new Thread(std::tr1::bind(&Implementation::spe_thread, this));
    }

    /// Destructor.
    ~Implementation()
    {
        // Kill the thread.
        delete thread;

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

    static void * spe_thread_function(Implementation * imp)
    {
        do
        {
            SPETask * task(0);

            // Dequeue task
            {
                Lock dequeue_lock(*imp->mutex);
                task = imp->dequeue();
            }

            if (task)
            {
                // task->run();
            }
            else
            {
                // Wait or sleep?
            }

            // Update exit condition
            {
                Lock query_lock(*imp->mutex);
                /// \todo
            }
        }
        while (true);
    }
};

SPE::SPE() :
    _imp(new Implementation)
{
}

SPE::SPE(const SPE & other) :
    _imp(other._imp)
{
}

SPE::~SPE()
{
}

DeviceId
SPE::id()
{
    return _imp->device;
}

struct SPEManager::Implementation
{
    typedef std::list<SPE> SPEList;

    /// Our number of available SPEs.
    unsigned spe_count;

    /// Out list of SPEs.
    SPEList spe_list;

    /// Constructor.
    Implementation() :
        spe_count(spe_cpu_info_get(SPE_COUNT_USABLE_SPES, -1)) /// \todo USABLE or PHYSICAL?
    {
    }
};

SPEManager::SPEManager() :
    _imp(new Implementation)
{
    unsigned count(_imp->spe_count);
    while(count-- > 0)
    {
        _imp->spe_list.push_back(SPE());
    }
}

SPEManager::~SPEManager()
{
}

SPEManager *
SPEManager::instance()
{
    static SPEManager result;

    return &result;
}

SPEManager::Iterator
SPEManager::begin() const
{
    return Iterator(_imp->spe_list.begin());
}

SPEManager::Iterator
SPEManager::end() const
{
    return Iterator(_imp->spe_list.end());
}
