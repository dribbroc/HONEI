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
#include <honei/backends/multicore/thread_pool.hh>
#include <honei/backends/multicore/ticket.hh>
#include <honei/util/attributes.hh>
#include <honei/util/exception.hh>
#include <honei/util/lock.hh>
#include <honei/util/log.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <honei/util/stringify.hh>
#include <honei/util/thread.hh>

#include <sys/syscall.h>
#include <math.h>

namespace honei
{
    /* TFImplementationBase is a base-class for all concrete implementations
     * of a ThreadFunction and provides the members that all derivatives
     * (have to) share. */
    struct TFImplementationBase
    {
        /// Our internal ID.
        const unsigned pool_id;

        /// Our Thread ID (given by the operating system).
        unsigned thread_id;

        /// The thread pool's mutex (for grabbing work).
        Mutex * const pool_mutex;

        /// ConditionVariable for work-status.
        ConditionVariable * const global_barrier;

        /// Flag if the Thread shall stop.
        volatile bool terminate;

        TFImplementationBase(mc::PoolSyncData * const psync, unsigned pid) :
            pool_id(pid),
            thread_id(0),
            pool_mutex(psync->mutex),
            global_barrier(psync->barrier),
            terminate(false)
        {
        }

        virtual ~TFImplementationBase()
        {
            bool loop(true);
            // Loop until the thread really stopped executing the
            // ThreadFunction object.
            do
            {
                loop = terminate;
            }
            while (loop);
        }

        void stop()
        {
            Lock l(*pool_mutex);
            terminate = true;
            global_barrier->broadcast();
        }
    };

    /* AffinityThreadFunction is the type of thread function
     * assigned to a pool thread if affinity is enabled. */
    template <> struct Implementation<mc::AffinityThreadFunction> :
        public TFImplementationBase
    {
        /* The logical processor this thread is bound to, in fact a
         * scheduler id which is only used with affinity enabled. */
        const unsigned sched_lpu;

        /// The task list administrated by the thread pool.
        std::list<mc::ThreadTask *> * const tasklist;

        Implementation(mc::PoolSyncData * const psync, std::list<mc::ThreadTask *> * const list,
                unsigned pid, unsigned sched) :
            TFImplementationBase(psync, pid),
            sched_lpu(sched),
            tasklist(list)
        {
        }

        ~Implementation()
        {
        }
    };

    template <> struct Implementation<mc::WorkStealingThreadFunction> :
        public TFImplementationBase
    {
        /* The logical processor this thread is bound to, in fact a
         * scheduler id which is only used with affinity enabled. */
        const unsigned sched_lpu;

        /// The task list (local to this thread!)
        std::list<mc::ThreadTask *> tasklist;

        /// Reference to the thread-vector of the thread pool (for stealing)
        const std::vector<std::pair<Thread *, mc::WorkStealingThreadFunction *> > & threads;

        /// Local mutual exclusion
        Mutex * const local_mutex;

        /// The overall number of pooled threads
        const unsigned num_threads;

        volatile bool & global_terminate;

        Implementation(mc::PoolSyncData * const psync, unsigned pid, unsigned sched,
                const std::vector<std::pair<Thread *, mc::WorkStealingThreadFunction *> > & thr, unsigned num_thr, volatile bool & term) :
            TFImplementationBase(psync, pid),
            sched_lpu(sched),
            threads(thr),
            local_mutex(new Mutex),
            num_threads(num_thr),
            global_terminate(term)
        {
        }

        ~Implementation()
        {
            delete local_mutex;
        }

        void enqueue(mc::ThreadTask * task)
        {
            Lock l(*local_mutex);
            tasklist.push_back(task);
        }

        bool steal(std::list<mc::ThreadTask *> & thief_list)
        {
            Lock l(*local_mutex);

            if (tasklist.empty())
                return false;
            else
            {
                int size = tasklist.size();
                std::list<mc::ThreadTask *>::iterator it(tasklist.begin());

                for (int i(0) ; i < size / 2 ; ++i)
                    ++it;

                thief_list.splice(thief_list.end(), tasklist, it, tasklist.end());
                return true;
            }
        }
    };

    /* SimpleThreadFunction is the type of thread function
     * assigned to a pool thread if affinity is disabled.
     * This allows for a simplified task list implementation
     * which can save substantial locking delays */

    // Tell the compiler that we want to instantiate AtomicSList
    // for a ThreadTask pointer
    template class AtomicSList<mc::ThreadTask *>;

    template <> struct Implementation<mc::SimpleThreadFunction> :
        public TFImplementationBase
    {
        /// The task list administrated by the thread pool.
        AtomicSList<mc::ThreadTask *> * const tasklist;

        Implementation(mc::PoolSyncData * const psync, AtomicSList<mc::ThreadTask *> * const list,
                unsigned pid) :
            TFImplementationBase(psync, pid),
            tasklist(list)
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
}

void AffinityThreadFunction::stop()
{
    _imp->stop();
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

    _imp->terminate = false; // Signal Implementation DTOR that we arrived here
}

unsigned AffinityThreadFunction::pool_id() const
{
    return _imp->pool_id;
}

unsigned AffinityThreadFunction::tid() const
{
    return _imp->thread_id;
}

WorkStealingThreadFunction::WorkStealingThreadFunction(PoolSyncData * const psync,
        unsigned pool_id, unsigned sched_id, const std::vector<std::pair<Thread *, mc::WorkStealingThreadFunction *> > & threads, unsigned num_thr, volatile bool & terminate) :
    PrivateImplementationPattern<WorkStealingThreadFunction, Shared>(new
            Implementation<WorkStealingThreadFunction>(psync, pool_id, sched_id, threads, num_thr, terminate))
{
}

WorkStealingThreadFunction::~WorkStealingThreadFunction()
{
}

void WorkStealingThreadFunction::stop()
{
    _imp->stop();
}

void WorkStealingThreadFunction::enqueue(mc::ThreadTask * task)
{
    _imp->enqueue(task);
}

void WorkStealingThreadFunction::operator() ()
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

        {
            Lock l(*_imp->local_mutex);

            if (! _imp->tasklist.empty())
            {
                task = _imp->tasklist.front();
                _imp->tasklist.pop_front();
            }
        }

        if (task == 0)
        {
            const int iter(rand() % _imp->num_threads);
            WorkStealingThreadFunction * const tfunc = _imp->threads[iter].second;

            pthread_mutex_lock(_imp->pool_mutex->mutex());
            pthread_mutex_lock(_imp->local_mutex->mutex());
            bool ok = (_imp->global_terminate ? false : tfunc->_imp->steal(_imp->tasklist));

            if (! ok && _imp->tasklist.empty())
            {
                // Have to make sure that no task has been added after unlocking the local
                // mutex. this is a possible race-condition to circumvent here!
                if (! _imp->terminate)
                {
                    pthread_mutex_unlock(_imp->local_mutex->mutex());
                    _imp->global_barrier->wait(*_imp->pool_mutex);
                    pthread_mutex_unlock(_imp->pool_mutex->mutex());
                }
                else
                {
                    pthread_mutex_unlock(_imp->pool_mutex->mutex());
                    pthread_mutex_unlock(_imp->local_mutex->mutex());
                    break;
                }
            }
            else
            {
                pthread_mutex_unlock(_imp->pool_mutex->mutex());
                task = _imp->tasklist.front();
                _imp->tasklist.pop_front();
                pthread_mutex_unlock(_imp->local_mutex->mutex());
            }
        }

        if (task != 0)
        {
            unsigned & sched_id = task->ticket->sid();
            sched_id = _imp->sched_lpu;
#ifdef DEBUG
            std::string msg = "Thread " + stringify(_imp->pool_id) + " on LPU " +
                stringify(_imp->sched_lpu) + " will execute ticket " +
                stringify(task->ticket->uid()) + "\n";
            LOGMESSAGE(lc_backend, msg);
#endif

            (*task->functor)();
            task->ticket->mark();
            delete task;
        }
    }
    while (true);

    delete comp;

    _imp->terminate = false; // Signal Implementation DTOR that we arrived here
}

unsigned WorkStealingThreadFunction::pool_id() const
{
    return _imp->pool_id;
}

unsigned WorkStealingThreadFunction::tid() const
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
}

void SimpleThreadFunction::stop()
{
    _imp->stop();
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

        {
            Lock l(*_imp->pool_mutex);

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

    _imp->terminate = false; // Signal Implementation DTOR that we arrived here
}

unsigned SimpleThreadFunction::tid() const
{
    return _imp->thread_id;
}
