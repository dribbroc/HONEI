/* vim: set sw=4 sts=4 et foldmethod=syntax : */
/*
 * Copyright (c) 2009, 2010, 2011 Sven Mallach <mallach@honei.org>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/backends/multicore/atomic_slist-impl.hh>
#include <honei/backends/multicore/thread_function.hh>
#include <honei/util/attributes.hh>
#include <honei/util/exception.hh>
#include <honei/util/lock.hh>
#include <honei/util/log.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <honei/util/stringify.hh>

#include <sys/syscall.h>

namespace honei
{
    /* AffinityThreadFunction is the type of thread function
     * assigned to a pool thread if affinity is enabled. */

    template <> struct Implementation<mc::AffinityThreadFunction>
    {
        /// Our internal ID.
        const unsigned pool_id;

        /// Our Thread ID (given by the operating system).
        unsigned thread_id;

        /* The logical processor this thread is bound to, in fact a
         * scheduler id which is only used with affinity enabled. */
        const unsigned sched_lpu;

        /// The thread pool's mutex (for grabbing work).
        Mutex * const pool_mutex;

        /// The task list administrated by the thread pool.
        std::list<mc::ThreadTask *> * const tasklist;

        /// ConditionVariable for work-status.
        ConditionVariable * const global_barrier;

        /// Flag if the Thread shall stop.
        volatile bool terminate;

        Implementation(mc::PoolSyncData * const psync, std::list<mc::ThreadTask *> * const list,
                unsigned pid, unsigned sched) :
            pool_id(pid),
            thread_id(0),
            sched_lpu(sched),
            pool_mutex(psync->mutex),
            tasklist(list),
            global_barrier(psync->barrier),
            terminate(false)
            {
            }

            ~Implementation()
            {
            }
    };

    /* SimpleThreadFunction is the type of thread function
     * assigned to a pool thread if affinity is disabled.
     * This allows for a simplified task list implementation
     * which can save substantial locking delays */

    // Tell the compiler that we want to instantiate AtomicSList
    // for a ThreadTask pointer
    template class AtomicSList<mc::ThreadTask *>;

    template <> struct Implementation<mc::SimpleThreadFunction>
    {
        /// Our internal ID.
        const unsigned pool_id;

        /// Our Thread ID (given by the operating system).
        unsigned thread_id;

        /// The thread pool's mutex (for grabbing work).
        Mutex * const pool_mutex;

        /// ConditionVariable for work-status.
        ConditionVariable * const global_barrier;

        /// The task list administrated by the thread pool.
        AtomicSList<mc::ThreadTask *> * const tasklist;

        /// Flag if the Thread shall stop.
        volatile bool terminate;

        Implementation(mc::PoolSyncData * const psync, AtomicSList<mc::ThreadTask *> * const list,
                unsigned pid) :
            pool_id(pid),
            thread_id(0),
            pool_mutex(psync->mutex),
            global_barrier(psync->barrier),
            tasklist(list),
            terminate(false)
            {
            }

            ~Implementation()
            {
            }
    };
}

using namespace honei::mc;

ThreadFunctionBase::~ThreadFunctionBase()
{
}

AffinityThreadFunction::AffinityThreadFunction(PoolSyncData * const psync,
        std::list<ThreadTask *> * const list, unsigned pool_id, unsigned sched_id) :
    PrivateImplementationPattern<AffinityThreadFunction, Shared>(new
            Implementation<AffinityThreadFunction>(psync, list, pool_id, sched_id))
{
}

AffinityThreadFunction::~AffinityThreadFunction()
{
    bool loop(true);
    // Loop until the thread really stopped executing this
    // AffinityThreadFunction object.
    do
    {
        loop = _imp->terminate;
    }
    while (loop);
}

void AffinityThreadFunction::stop()
{
    Lock l(*_imp->pool_mutex);
    _imp->terminate = true;
    _imp->global_barrier->broadcast();
}

void AffinityThreadFunction::operator() ()
{
    /// The ThreadTask to be currently executed by the thread.
    mc::ThreadTask * task(0);

    /// A comparison object for mc::ThreadTask objects.
    mc::TaskComp * const comp(new mc::TaskComp(_imp->sched_lpu));

    std::list<mc::ThreadTask *>::iterator i, i_end;

    /* Set thread_id from operating system and use this on the pool
     * side as sign for the thread to be setup. Then let the thread
     * go to sleep until it will be assigned its first task. */

    {
        Lock l(*_imp->pool_mutex);
        _imp->thread_id = syscall(__NR_gettid);
        _imp->global_barrier->wait(*_imp->pool_mutex);
    }

    do
    {
        task = 0;

        // Check work-status and condition again under same mutex protection
        // to avoid race conditions between pick / if and wait
        {
            Lock l(*_imp->pool_mutex); // Try to avoid this lock...

            for (i = _imp->tasklist->begin(), i_end = _imp->tasklist->end() ; i != i_end ; ++i)
            {
                if ((*comp)(*i))
                {
                    task = *i;
                    unsigned & sched_id = task->ticket->sid();
                    sched_id = _imp->sched_lpu;
                    _imp->tasklist->remove(*i);
#ifdef DEBUG
                    std::string msg = "Thread " + stringify(_imp->pool_id) + " on LPU " +
                        stringify(_imp->sched_lpu) + " will execute ticket " +
                        stringify(task->ticket->uid()) + "\n";
                    LOGMESSAGE(lc_backend, msg);
#endif
                    break;
                }
            }

            if (task == 0)
            {
                if (! _imp->terminate)
                    _imp->global_barrier->wait(*_imp->pool_mutex);
                else
                    break;
            }
        }

        if (task != 0)
        {
            (*task->functor)();
            task->ticket->mark();
            delete task;
        }
    }
    while (true);

    delete comp;

    _imp->terminate = false; // Signal AffinityThreadFunction DTOR that we arrived here
}

unsigned AffinityThreadFunction::pool_id() const
{
    return _imp->pool_id;
}

unsigned AffinityThreadFunction::tid() const
{
    return _imp->thread_id;
}

SimpleThreadFunction::SimpleThreadFunction(PoolSyncData * const psync, AtomicSList<ThreadTask *> * const list,
        unsigned pool_id) :
    PrivateImplementationPattern<SimpleThreadFunction, Shared>(new
            Implementation<SimpleThreadFunction>(psync, list, pool_id))
{
}

SimpleThreadFunction::~SimpleThreadFunction()
{
    bool loop(true);
    // Loop until the thread really stopped executing this
    // SimpleThreadFunction object.
    do
    {
        loop = _imp->terminate;
    }
    while (loop);
}

void SimpleThreadFunction::stop()
{
    Lock l(*_imp->pool_mutex);
    _imp->terminate = true;
    _imp->global_barrier->broadcast();
}

void SimpleThreadFunction::operator() ()
{
    /// The ThreadTask to be currently executed by the thread.
    mc::ThreadTask * task(0);

    /* Set thread_id from operating system and use this on the pool
     * side as sign for the thread to be setup. Then let the thread
     * go to sleep until it will be assigned its first task. */

    {
        Lock l(*_imp->pool_mutex);
        _imp->thread_id = syscall(__NR_gettid);
        _imp->global_barrier->wait(*_imp->pool_mutex);
    }

    do
    {
        task = _imp->tasklist->pop_front();

        if (task == 0)
        {
            Lock l(*_imp->pool_mutex);

            if (! _imp->terminate)
                _imp->global_barrier->wait(*_imp->pool_mutex);
            else
                break;
        }

        if (task != 0)
        {
#ifdef DEBUG
            std::string msg = "Thread " + stringify(_imp->pool_id) + " will execute ticket " +
                stringify(task->ticket->uid()) + "\n";
            LOGMESSAGE(lc_backend, msg);
#endif
            (*task->functor)();
            task->ticket->mark();
            delete task;
        }
    }
    while (true);

    _imp->terminate = false; // Signal SimpleThreadFunction DTOR that we arrived here
}

unsigned SimpleThreadFunction::tid() const
{
    return _imp->thread_id;
}
