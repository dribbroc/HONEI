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

#include <honei/backends/multicore/concurrent_list-impl.hh>
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
#include <deque>

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

        virtual void operator() () = 0;

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
        std::deque<mc::ThreadTask *> * const tasklist;

        Implementation(mc::PoolSyncData * const psync, std::deque<mc::ThreadTask *> * const list,
                unsigned pid, unsigned sched) :
            TFImplementationBase(psync, pid),
            sched_lpu(sched),
            tasklist(list)
        {
        }

        ~Implementation()
        {
        }

        void operator() ()
        {
            /// The ThreadTask to be currently executed by the thread.
            mc::ThreadTask * task(0);

            /// A comparison object for mc::ThreadTask objects.
            mc::TaskComp * const comp(new mc::TaskComp(sched_lpu));

            std::deque<mc::ThreadTask *>::iterator i, i_end;

        /* Set thread_id from operating system and use this on the pool
         * side as sign for the thread to be setup. Then let the thread
         * go to sleep until it will be assigned its first task. */

            {
                Lock l(*pool_mutex);
                thread_id = syscall(__NR_gettid);
                global_barrier->wait(*pool_mutex);
            }

            do
            {
                task = 0;

                // Check work-status and condition again under same mutex protection
                // to avoid race conditions between pick / if and wait
                {
                    Lock l(*pool_mutex); // Try to avoid this lock...

                    for (i = tasklist->begin(), i_end = tasklist->end() ; i != i_end ; ++i)
                    {
                        if ((*comp)(*i))
                        {
                            task = *i;
                            unsigned & sched_id = task->ticket->sid();
                            sched_id = sched_lpu;
                            tasklist->erase(i);
#ifdef DEBUG
                            std::string msg = "Thread " + stringify(pool_id) + " on LPU " +
                            stringify(sched_lpu) + " will execute ticket " +
                            stringify(task->ticket->uid()) + "\n";
                            LOGMESSAGE(lc_backend, msg);
#endif
                            break;
                        }
                    }

                    if (task == 0)
                    {
                        if (! terminate)
                            global_barrier->wait(*pool_mutex);
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

            terminate = false; // Signal Implementation DTOR that we arrived here
        }
    };

    template <> struct Implementation<mc::WorkStealingThreadFunction<mc::ConcurrentList<mc::ThreadTask *> > > :
        public TFImplementationBase
    {
        /* The logical processor this thread is bound to, in fact a
         * scheduler id which is only used with affinity enabled. */
        const unsigned sched_lpu;

        /// The task list (local to this thread!)
        ConcurrentList<mc::ThreadTask *> tasklist;

        /// Reference to the thread-vector of the thread pool (for stealing)
        const std::vector<std::pair<Thread *, mc::WorkStealingThreadFunction<mc::ConcurrentList<mc::ThreadTask *> > * > > & threads;

        /// Local mutual exclusion
        Mutex * const local_mutex;

        /// The overall number of pooled threads
        const unsigned num_threads;

        volatile bool & global_terminate;

        Mutex * const steal_mutex;

        Implementation(mc::PoolSyncData * const psync, unsigned pid, unsigned sched,
                const std::vector<std::pair<Thread *, mc::WorkStealingThreadFunction<mc::ConcurrentList<mc::ThreadTask *> > *> > & thr,
                unsigned num_thr, volatile bool & term, Mutex * const stm) :
            TFImplementationBase(psync, pid),
            sched_lpu(sched),
            threads(thr),
            local_mutex(new Mutex),
            num_threads(num_thr),
            global_terminate(term),
            steal_mutex(stm)
        {
        }

        ~Implementation()
        {
            delete local_mutex;
        }

        void operator() ()
        {
            /// The ThreadTask to be currently executed by the thread.
            mc::ThreadTask * task(0);

            /* Set thread_id from operating system and use this on the pool
             * side as sign for the thread to be setup. Then let the thread
             * go to sleep until it will be assigned its first task. */

            {
                Lock l(*pool_mutex);
                thread_id = syscall(__NR_gettid);
                global_barrier->wait(*pool_mutex);
            }

            do
            {
                task = tasklist.pop_front();

                if (task == 0)
                {
                    pthread_mutex_lock(steal_mutex->mutex());
                    pthread_mutex_lock(local_mutex->mutex());

                    if (! global_terminate)
                    {
                        const int iter(rand() % num_threads);
                        mc::WorkStealingThreadFunction<mc::ConcurrentList<mc::ThreadTask *> > * const tfunc = threads[iter].second;
                        tfunc->steal(tasklist);
                        pthread_mutex_unlock(steal_mutex->mutex());

                        task = tasklist.pop_front();

                        if (task == 0)
                        {
                            pthread_mutex_lock(pool_mutex->mutex());
                            pthread_mutex_unlock(local_mutex->mutex());
                            global_barrier->wait(*pool_mutex);
                            pthread_mutex_unlock(pool_mutex->mutex());
                        }
                        else
                            pthread_mutex_unlock(local_mutex->mutex());
                    }
                    else
                    {
                        pthread_mutex_unlock(steal_mutex->mutex());

                        task = tasklist.pop_front();
                        if (task == 0 && terminate)
                        {
                            pthread_mutex_unlock(local_mutex->mutex());
                            break;
                        }

                        pthread_mutex_unlock(local_mutex->mutex());
                    }
                }

                if (task != 0)
                {
                    unsigned & sched_id = task->ticket->sid();
                    sched_id = sched_lpu;
#ifdef DEBUG
                    std::string msg = "Thread " + stringify(pool_id) + " on LPU " +
                    stringify(sched_lpu) + " will execute ticket " +
                    stringify(task->ticket->uid()) + "\n";
                    LOGMESSAGE(lc_backend, msg);
#endif

                    (*task->functor)();
                    task->ticket->mark();
                    delete task;
                }

            }
            while (true);

            terminate = false; // Signal Implementation DTOR that we arrived here
        }

        void enqueue(mc::ThreadTask * task)
        {
            Lock l(*local_mutex);
            tasklist.push_back(task);
        }

        bool steal(ConcurrentList<mc::ThreadTask *> & thief_list)
        {
            if (tasklist.empty())
                return false;
            else
            {
                ThreadTask * t = tasklist.pop_front();
                thief_list.push_back(t);
                return true;
            }
        }
    };

    template <> struct Implementation<mc::WorkStealingThreadFunction<std::deque<mc::ThreadTask *> > > :
        public TFImplementationBase
    {
        /* The logical processor this thread is bound to, in fact a
         * scheduler id which is only used with affinity enabled. */
        const unsigned sched_lpu;

        /// The task list (local to this thread!)
        std::deque<mc::ThreadTask *> tasklist;

        /// Reference to the thread-vector of the thread pool (for stealing)
        const std::vector<std::pair<Thread *, mc::WorkStealingThreadFunction<std::deque<mc::ThreadTask *> > *> > & threads;

        /// Local mutual exclusion
        Mutex * const local_mutex;

        /// The overall number of pooled threads
        const unsigned num_threads;

        volatile bool & global_terminate;

        Mutex * const steal_mutex;

        Implementation(mc::PoolSyncData * const psync, unsigned pid, unsigned sched,
                const std::vector<std::pair<Thread *, mc::WorkStealingThreadFunction<std::deque<mc::ThreadTask *> > *> > & thr,
                unsigned num_thr, volatile bool & term, Mutex * const stm) :
            TFImplementationBase(psync, pid),
            sched_lpu(sched),
            threads(thr),
            local_mutex(new Mutex),
            num_threads(num_thr),
            global_terminate(term),
            steal_mutex(stm)
        {
        }

        ~Implementation()
        {
            delete local_mutex;
        }

        void operator() ()
        {
            /// The ThreadTask to be currently executed by the thread.
            mc::ThreadTask * task(0);

            /* Set thread_id from operating system and use this on the pool
             * side as sign for the thread to be setup. Then let the thread
             * go to sleep until it will be assigned its first task. */

            {
                Lock l(*pool_mutex);
                thread_id = syscall(__NR_gettid);
                global_barrier->wait(*pool_mutex);
            }

            do
            {
                task = 0;

                {
                    Lock l(*local_mutex);

                    if (! tasklist.empty())
                    {
                        task = tasklist.front();
                        tasklist.pop_front();
                    }
                }

                if (task == 0)
                {
                    pthread_mutex_lock(steal_mutex->mutex());
                    pthread_mutex_lock(local_mutex->mutex());

                    if (! tasklist.empty())
                    {
                        task = tasklist.front();
                        tasklist.pop_front();
                    }
                    else if (! global_terminate)
                    {
                        const int iter(rand() % num_threads);
                        mc::WorkStealingThreadFunction<std::deque<mc::ThreadTask *> > * const tfunc = threads[iter].second;
                        if (tfunc->steal(tasklist))
                        {
                            task = tasklist.front();
                            tasklist.pop_front();
                        }
                    }
                    pthread_mutex_unlock(steal_mutex->mutex());

                    if (task == 0 && ! (global_terminate || terminate))
                    {
                        pthread_mutex_lock(pool_mutex->mutex());
                        pthread_mutex_unlock(local_mutex->mutex());
                        global_barrier->wait(*pool_mutex);
                        pthread_mutex_unlock(pool_mutex->mutex());
                    }
                    else if (terminate)
                    {
                        pthread_mutex_unlock(local_mutex->mutex());
                        break;
                    }
                    else
                    {
                        pthread_mutex_unlock(local_mutex->mutex());
                    }
                }

                if (task != 0)
                {
                    unsigned & sched_id = task->ticket->sid();
                    sched_id = sched_lpu;
#ifdef DEBUG
                    std::string msg = "Thread " + stringify(pool_id) + " on LPU " +
                    stringify(sched_lpu) + " will execute ticket " +
                    stringify(task->ticket->uid()) + "\n";
                    LOGMESSAGE(lc_backend, msg);
#endif

                    (*task->functor)();
                    task->ticket->mark();
                    delete task;
                }
            }
            while (true);

            terminate = false; // Signal Implementation DTOR that we arrived here
        }

        void enqueue(mc::ThreadTask * task)
        {
            Lock l(*local_mutex);
            tasklist.push_back(task);
        }

        bool steal(std::deque<mc::ThreadTask *> & thief_list)
        {
            Lock l(*local_mutex);

            if (tasklist.empty())
                return false;
            else
            {
                int size = tasklist.size();

                switch (size)
                {
                    case 1:
                    {
                        mc::ThreadTask * t = tasklist.back();
                        tasklist.pop_back();
                        thief_list.push_back(t);
                        break;
                    }

                    default:
                    {
                        for (int i(0) ; i < (size >> 1) ; ++i)
                        {
                            mc::ThreadTask * t = tasklist.back();
                            tasklist.pop_back();
                            thief_list.push_back(t);
                        }
                    }

                }

                return true;
            }
        }
    };

    /* SimpleThreadFunction is the type of thread function
     * assigned to a pool thread if affinity is disabled.
     * This allows for a simplified task list implementation
     * which can save substantial locking delays */

    // Tell the compiler that we want to instantiate ConcurrentList
    // for a ThreadTask pointer
    template class ConcurrentList<mc::ThreadTask *>;

    template <> struct Implementation<mc::SimpleThreadFunction<ConcurrentList<mc::ThreadTask *> > > :
        public TFImplementationBase
    {
        /// The task list administrated by the thread pool.
        ConcurrentList<mc::ThreadTask *> * const tasklist;

        Implementation(mc::PoolSyncData * const psync, ConcurrentList<mc::ThreadTask *> * const list,
                unsigned pid) :
            TFImplementationBase(psync, pid),
            tasklist(list)
        {
        }

        ~Implementation()
        {
        }

        void operator () ()
        {
            /// The ThreadTask to be currently executed by the thread.
            mc::ThreadTask * task(0);

            /* Set thread_id from operating system and use this on the pool
             * side as sign for the thread to be setup. Then let the thread
             * go to sleep until it will be assigned its first task. */

            {
                Lock l(*pool_mutex);
                thread_id = syscall(__NR_gettid);
                global_barrier->wait(*pool_mutex);
            }

            do
            {
                task = tasklist->pop_front();
                {
                    Lock l(*pool_mutex);

                    if (task == 0)
                    {
                        if (! terminate)
                            global_barrier->wait(*pool_mutex);
                        else
                            break;
                    }
                }

                if (task != 0)
                {
#ifdef DEBUG
                    std::string msg = "Thread " + stringify(pool_id) + " will execute ticket " +
                    stringify(task->ticket->uid()) + "\n";
                    LOGMESSAGE(lc_backend, msg);
#endif
                    (*task->functor)();
                    task->ticket->mark();
                    delete task;
                }
            }
            while (true);

            terminate = false; // Signal Implementation DTOR that we arrived here
        }
    };

    template <> struct Implementation<mc::SimpleThreadFunction<std::deque<mc::ThreadTask *> > > :
        public TFImplementationBase
    {
        /// The task list administrated by the thread pool.
        std::deque<mc::ThreadTask *> * const tasklist;

        Implementation(mc::PoolSyncData * const psync, std::deque<mc::ThreadTask *> * const list,
                unsigned pid) :
            TFImplementationBase(psync, pid),
            tasklist(list)
        {
        }

        ~Implementation()
        {
        }

        void operator () ()
        {
            /// The ThreadTask to be currently executed by the thread.
            mc::ThreadTask * task(0);

            /* Set thread_id from operating system and use this on the pool
             * side as sign for the thread to be setup. Then let the thread
             * go to sleep until it will be assigned its first task. */

            {
                Lock l(*pool_mutex);
                thread_id = syscall(__NR_gettid);
                global_barrier->wait(*pool_mutex);
            }

            do
            {
                {
                    Lock l(*pool_mutex);

                    task = (tasklist->empty() ? 0 : tasklist->front());

                    if (task == 0)
                    {
                        if (! terminate)
                            global_barrier->wait(*pool_mutex);
                        else
                            break;
                    }
                    else
                        tasklist->pop_front();
                }

                if (task != 0)
                {
#ifdef DEBUG
                    std::string msg = "Thread " + stringify(pool_id) + " will execute ticket " +
                    stringify(task->ticket->uid()) + "\n";
                    LOGMESSAGE(lc_backend, msg);
#endif
                    (*task->functor)();
                    task->ticket->mark();
                    delete task;
                }
            }
            while (true);

            terminate = false; // Signal Implementation DTOR that we arrived here
        }
    };
}

using namespace honei::mc;

ThreadFunctionBase::~ThreadFunctionBase()
{
}

AffinityThreadFunction::AffinityThreadFunction(PoolSyncData * const psync,
        std::deque<ThreadTask *> * const list, unsigned pool_id, unsigned sched_id) :
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
    (*_imp)();
}

unsigned AffinityThreadFunction::pool_id() const
{
    return _imp->pool_id;
}

unsigned AffinityThreadFunction::tid() const
{
    return _imp->thread_id;
}

WorkStealingThreadFunction<std::deque<mc::ThreadTask *> >::WorkStealingThreadFunction(PoolSyncData * const psync, unsigned pool_id, unsigned sched_id,
        const std::vector<std::pair<Thread *, mc::WorkStealingThreadFunction<std::deque<mc::ThreadTask *> > *> > & threads,
        unsigned num_thr, volatile bool & terminate, Mutex * const steal_mutex) :
    PrivateImplementationPattern<WorkStealingThreadFunction<std::deque<mc::ThreadTask *> >, Shared>(new
            Implementation<WorkStealingThreadFunction<std::deque<mc::ThreadTask *> > >(psync, pool_id, sched_id, threads, num_thr, terminate, steal_mutex))
{
}

WorkStealingThreadFunction<std::deque<mc::ThreadTask *> >::~WorkStealingThreadFunction()
{
}

void WorkStealingThreadFunction<std::deque<mc::ThreadTask *> >::stop()
{
    _imp->stop();
}

void WorkStealingThreadFunction<std::deque<mc::ThreadTask *> >::enqueue(mc::ThreadTask * task)
{
    _imp->enqueue(task);
}

bool WorkStealingThreadFunction<std::deque<mc::ThreadTask *> >::steal(std::deque<mc::ThreadTask *> & thief_list)
{
    return _imp->steal(thief_list);
}

void WorkStealingThreadFunction<std::deque<mc::ThreadTask *> >::operator() ()
{
    (*_imp)();
}

unsigned WorkStealingThreadFunction<std::deque<mc::ThreadTask *> >::pool_id() const
{
    return _imp->pool_id;
}

unsigned WorkStealingThreadFunction<std::deque<mc::ThreadTask *> >::tid() const
{
    return _imp->thread_id;
}

WorkStealingThreadFunction<mc::ConcurrentList<mc::ThreadTask *> >::WorkStealingThreadFunction(PoolSyncData * const psync, unsigned pool_id, unsigned sched_id,
        const std::vector<std::pair<Thread *, mc::WorkStealingThreadFunction<mc::ConcurrentList<mc::ThreadTask *> > *> > & threads,
        unsigned num_thr, volatile bool & terminate, Mutex * const steal_mutex) :
    PrivateImplementationPattern<WorkStealingThreadFunction<mc::ConcurrentList<mc::ThreadTask *> >, Shared>(new
            Implementation<WorkStealingThreadFunction<mc::ConcurrentList<mc::ThreadTask *> > >(psync, pool_id, sched_id, threads, num_thr, terminate, steal_mutex))
{
}

WorkStealingThreadFunction<mc::ConcurrentList<mc::ThreadTask *> >::~WorkStealingThreadFunction()
{
}

void WorkStealingThreadFunction<mc::ConcurrentList<mc::ThreadTask *> >::stop()
{
    _imp->stop();
}

void WorkStealingThreadFunction<mc::ConcurrentList<mc::ThreadTask *> >::enqueue(mc::ThreadTask * task)
{
    _imp->enqueue(task);
}

bool WorkStealingThreadFunction<mc::ConcurrentList<mc::ThreadTask *> >::steal(mc::ConcurrentList<mc::ThreadTask *> & thief_list)
{
    return _imp->steal(thief_list);
}

void WorkStealingThreadFunction<mc::ConcurrentList<mc::ThreadTask *> >::operator() ()
{
    (*_imp)();
}

unsigned WorkStealingThreadFunction<mc::ConcurrentList<mc::ThreadTask *> >::pool_id() const
{
    return _imp->pool_id;
}

unsigned WorkStealingThreadFunction<mc::ConcurrentList<mc::ThreadTask *> >::tid() const
{
    return _imp->thread_id;
}


SimpleThreadFunction<ConcurrentList<ThreadTask *> >::SimpleThreadFunction(PoolSyncData * const psync,
        ConcurrentList<ThreadTask *> * const list, unsigned pool_id) :
    PrivateImplementationPattern<SimpleThreadFunction<ConcurrentList<ThreadTask *> >, Shared>(new
            Implementation<SimpleThreadFunction<ConcurrentList<ThreadTask *> > >(psync, list, pool_id))
{
}

SimpleThreadFunction<ConcurrentList<ThreadTask *> >::~SimpleThreadFunction()
{
}

void SimpleThreadFunction<ConcurrentList<ThreadTask *> >::stop()
{
    _imp->stop();
}

void SimpleThreadFunction<ConcurrentList<ThreadTask *> >::operator() ()
{
    (*_imp)();
}

unsigned SimpleThreadFunction<ConcurrentList<ThreadTask *> >::tid() const
{
    return _imp->thread_id;
}


SimpleThreadFunction<std::deque<ThreadTask *> >::SimpleThreadFunction(PoolSyncData * const psync,
        std::deque<ThreadTask *> * const list, unsigned pool_id) :
    PrivateImplementationPattern<SimpleThreadFunction<std::deque<ThreadTask *> >, Shared>(new
            Implementation<SimpleThreadFunction<std::deque<ThreadTask *> > >(psync, list, pool_id))
{
}

SimpleThreadFunction<std::deque<ThreadTask *> >::~SimpleThreadFunction()
{
}

void SimpleThreadFunction<std::deque<ThreadTask *> >::stop()
{
    _imp->stop();
}

void SimpleThreadFunction<std::deque<ThreadTask *> >::operator() ()
{
    (*_imp)();
}

unsigned SimpleThreadFunction<std::deque<ThreadTask *> >::tid() const
{
    return _imp->thread_id;
}
