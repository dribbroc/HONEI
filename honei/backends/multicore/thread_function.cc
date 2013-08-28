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

#include <honei/backends/multicore/cas_deque-impl.hh>
#include <honei/backends/multicore/concurrent_deque-impl.hh>
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

#include <deque>
#include <math.h>
#include <sys/syscall.h>
#include <unistd.h>

namespace honei
{
    namespace mc
    {
    // Tell the compiler that we want to instantiate CASDeque
    // and ConcurrentDeque for a ThreadTask pointer
    template class ConcurrentDeque<mc::ThreadTask *>;
    template class CASDeque<mc::ThreadTask *>;
    }

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

        /// Thread-specific data
        ThreadData * const thread_data;

        TFImplementationBase(mc::PoolSyncData * const psync, ThreadData * const tdata, unsigned pid) :
            pool_id(pid),
            thread_id(0),
            pool_mutex(psync->mutex),
            global_barrier(psync->barrier),
            thread_data(tdata)
        {
        }

        virtual ~TFImplementationBase()
        {
        }

        virtual void operator() () = 0;
    };

    /* AffinityThreadFunction is the type of thread function
     * assigned to a pool thread if affinity is enabled. */
    template <> struct Implementation<mc::AffinityThreadFunction> :
        public TFImplementationBase
    {
        /// The logical processor this thread is bound to
        LPU * const lpu;

        /// The task list administrated by the thread pool.
        std::deque<mc::ThreadTask *> * const tasklist;

        Implementation(mc::PoolSyncData * const psync, ThreadData * const tdata, std::deque<mc::ThreadTask *> * const list,
                unsigned pid, LPU * const pu) :
            TFImplementationBase(psync, tdata, pid),
            lpu(pu),
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
            mc::TaskComp * const comp(new mc::TaskComp(lpu));

            std::deque<mc::ThreadTask *>::iterator i, i_end;

            /* Set thread_id from operating system and signal the pool
             * that this thread is online. Then let the thread go to
             * sleep until it will be assigned its first task. */

            {
                Lock l(*pool_mutex);
                thread_id = syscall(__NR_gettid);
                global_barrier->broadcast();
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
                            unsigned & sched_id = task->ticket.sid();
                            sched_id = lpu->sched_id;
                            tasklist->erase(i);

#ifdef DEBUG
                            std::string msg = "Thread " + stringify(pool_id) + " on LPU " +
                            stringify(lpu->sched_id) + " will execute ticket " +
                            stringify(task->ticket.uid()) + "\n";
                            LOGMESSAGE(lc_backend, msg);
#endif
                            break;
                        }
                    }

                    if (task == 0)
                    {
                        if (! thread_data->terminate)
                            global_barrier->wait(*pool_mutex);
                        else
                            break;
                    }
                }

                if (task != 0)
                {
                    (*task->functor)();
                    task->ticket.mark();
                    delete task;
                }
            }
            while (true);

            delete comp;

            // The thread function's DTOR will not be called before this
            // thread has stopped execution here.
        }
    };

    template <> struct Implementation<mc::WorkStealingThreadFunction<mc::CASDeque<mc::ThreadTask *> > > :
        public TFImplementationBase
    {
        /// The logical processor this thread is bound to (if any)
        const unsigned sched_id;

        /// The task list (local to this thread!)
        CASDeque<mc::ThreadTask *> tasklist;

        /// Reference to the thread-vector of the thread pool (for stealing)
        const std::vector<mc::WorkStealingThreadFunction<mc::CASDeque<mc::ThreadTask *> > * > & threads;

        /// The overall number of pooled threads
        const unsigned num_threads;

        volatile bool & global_terminate;

        Mutex * const steal_mutex;

        Implementation(mc::PoolSyncData * const psync, ThreadData * const tdata, unsigned pid, unsigned sid,
                const std::vector<mc::WorkStealingThreadFunction<mc::CASDeque<mc::ThreadTask *> > *> & thr,
                unsigned num_thr, volatile bool & term) :
            TFImplementationBase(psync, tdata, pid),
            sched_id(sid),
            threads(thr),
            num_threads(num_thr),
            global_terminate(term),
            steal_mutex(psync->steal_mutex)
        {
        }

        ~Implementation()
        {
        }

        void operator() ()
        {
            /// The ThreadTask to be currently executed by the thread.
            mc::ThreadTask * task(0);

            /* Set thread_id from operating system and signal the pool
             * that this thread is online. Then let the thread go to
             * sleep until it will be assigned its first task. */
            {
                Lock l(*pool_mutex);
                thread_id = syscall(__NR_gettid);
                global_barrier->broadcast();
            }
            {
                // Now wait until all thread are online - otherwise
                // we might to try to steal from trashed memory...
                Lock ll(*steal_mutex);
            }

            do
            {
                task = tasklist.pop_front();

                if (task == 0)
                {
                    if (! global_terminate)
                    {
                        const int iter(rand() % num_threads);
                        mc::WorkStealingThreadFunction<mc::CASDeque<mc::ThreadTask *> > * const tfunc = threads[iter];
                        tfunc->steal(tasklist);

                        pthread_mutex_lock(pool_mutex->mutex());
                        task = tasklist.pop_front();

                        if (task == 0 && ! global_terminate)
                        {
                            global_barrier->wait(*pool_mutex);
                            pthread_mutex_unlock(pool_mutex->mutex());
                        }
                        else
                            pthread_mutex_unlock(pool_mutex->mutex());
                    }
                    else
                    {
                        task = tasklist.pop_front();
                        if (task == 0 && thread_data->terminate)
                        {
                            break;
                        }
                    }
                }

                if (task != 0)
                {
                    unsigned & tsched_id = task->ticket.sid();
                    tsched_id = sched_id;
#ifdef DEBUG
                    std::string msg = "Thread " + stringify(pool_id) + " on LPU " +
                    stringify(sched_id) + " will execute ticket " +
                    stringify(task->ticket.uid()) + "\n";
                    LOGMESSAGE(lc_backend, msg);
#endif

                    (*task->functor)();
                    task->ticket.mark();
                    delete task;
                }

            }
            while (true);

            // The thread function's DTOR will not be called before this
            // thread has stopped execution here.
        }

        void enqueue(mc::ThreadTask * task)
        {
            tasklist.push_back(task);
        }

        bool steal(CASDeque<mc::ThreadTask *> & thief_list)
        {
            ThreadTask * t = tasklist.pop_back();

            if (t == NULL)
                return false;
            else
            {
                thief_list.push_back(t);
                return true;
            }
        }
    };

    template <> struct Implementation<mc::WorkStealingThreadFunction<mc::ConcurrentDeque<mc::ThreadTask *> > > :
        public TFImplementationBase
    {
        /// The logical processor this thread is bound to (if any)
        const unsigned sched_id;

        /// The task list (local to this thread!)
        ConcurrentDeque<mc::ThreadTask *> tasklist;

        /// Reference to the thread-vector of the thread pool (for stealing)
        const std::vector<mc::WorkStealingThreadFunction<mc::ConcurrentDeque<mc::ThreadTask *> > * > & threads;

        /// The overall number of pooled threads
        const unsigned num_threads;

        volatile bool & global_terminate;

        Mutex * const steal_mutex;

        Implementation(mc::PoolSyncData * const psync, ThreadData * const tdata, unsigned pid, unsigned sid,
                const std::vector<mc::WorkStealingThreadFunction<mc::ConcurrentDeque<mc::ThreadTask *> > *> & thr,
                unsigned num_thr, volatile bool & term) :
            TFImplementationBase(psync, tdata, pid),
            sched_id(sid),
            threads(thr),
            num_threads(num_thr),
            global_terminate(term),
            steal_mutex(psync->steal_mutex)
        {
        }

        ~Implementation()
        {
        }

        void operator() ()
        {
            /// The ThreadTask to be currently executed by the thread.
            mc::ThreadTask * task(0);

            /* Set thread_id from operating system and signal the pool
             * that this thread is online. Then let the thread go to
             * sleep until it will be assigned its first task. */
            {
                Lock l(*pool_mutex);
                thread_id = syscall(__NR_gettid);
                global_barrier->broadcast();
            }
            {
                // Now wait until all thread are online - otherwise
                // we might to try to steal from trashed memory...
                Lock ll(*steal_mutex);
            }

            do
            {
                task = tasklist.pop_front();

                if (task == 0)
                {
                    if (! global_terminate)
                    {
                        const int iter(rand() % num_threads);
                        mc::WorkStealingThreadFunction<mc::ConcurrentDeque<mc::ThreadTask *> > * const tfunc = threads[iter];
                        tfunc->steal(tasklist);

                        pthread_mutex_lock(pool_mutex->mutex());
                        task = tasklist.pop_front();

                        if (task == 0 && ! global_terminate)
                        {
                            global_barrier->wait(*pool_mutex);
                            pthread_mutex_unlock(pool_mutex->mutex());
                        }
                        else
                            pthread_mutex_unlock(pool_mutex->mutex());
                    }
                    else
                    {
                        task = tasklist.pop_front();
                        if (task == 0 && thread_data->terminate)
                        {
                            break;
                        }
                    }
                }

                if (task != 0)
                {
                    unsigned & tsched_id = task->ticket.sid();
                    tsched_id = sched_id;
#ifdef DEBUG
                    std::string msg = "Thread " + stringify(pool_id) + " on LPU " +
                    stringify(sched_id) + " will execute ticket " +
                    stringify(task->ticket.uid()) + "\n";
                    LOGMESSAGE(lc_backend, msg);
#endif

                    (*task->functor)();
                    task->ticket.mark();
                    delete task;
                }

            }
            while (true);

            // The thread function's DTOR will not be called before this
            // thread has stopped execution here.
        }

        void enqueue(mc::ThreadTask * task)
        {
            tasklist.push_back(task);
        }

        bool steal(ConcurrentDeque<mc::ThreadTask *> & thief_list)
        {
            ThreadTask * t = tasklist.pop_back();

            if (t == NULL)
                return false;
            else
            {
                thief_list.push_back(t);
                return true;
            }
        }
    };

    template <> struct Implementation<mc::WorkStealingThreadFunction<std::deque<mc::ThreadTask *> > > :
        public TFImplementationBase
    {
        /// The logical processor this thread is bound to (if any)
        const unsigned sched_id;

        /// The task list (local to this thread!)
        std::deque<mc::ThreadTask *> tasklist;

        /// Reference to the thread-vector of the thread pool (for stealing)
        const std::vector<mc::WorkStealingThreadFunction<std::deque<mc::ThreadTask *> > *> & threads;

        /// The overall number of pooled threads
        const unsigned num_threads;

        volatile bool & global_terminate;

        Mutex * const steal_mutex;

        Implementation(mc::PoolSyncData * const psync, ThreadData * const tdata, unsigned pid, unsigned sid,
                const std::vector<mc::WorkStealingThreadFunction<std::deque<mc::ThreadTask *> > *> & thr,
                unsigned num_thr, volatile bool & term) :
            TFImplementationBase(psync, tdata, pid),
            sched_id(sid),
            threads(thr),
            num_threads(num_thr),
            global_terminate(term),
            steal_mutex(psync->steal_mutex)
        {
        }

        ~Implementation()
        {
        }

        void operator() ()
        {
            /// The ThreadTask to be currently executed by the thread.
            mc::ThreadTask * task(0);

            /* Set thread_id from operating system and signal the pool
             * that this thread is online. Then let the thread go to
             * sleep until it will be assigned its first task. */
            {
                Lock l(*pool_mutex);
                thread_id = syscall(__NR_gettid);
                global_barrier->broadcast();
            }
            {
                // Now wait until all thread are online - otherwise
                // we might to try to steal from trashed memory...
                Lock ll(*steal_mutex);
            }

            do
            {
                task = 0;

                {
                    Lock l(*thread_data->local_mutex);

                    if (! tasklist.empty())
                    {
                        task = tasklist.front();
                        tasklist.pop_front();
                    }
                }

                if (task == 0)
                {
                    pthread_mutex_lock(steal_mutex->mutex());
                    pthread_mutex_lock(thread_data->local_mutex->mutex());

                    if (! tasklist.empty())
                    {
                        task = tasklist.front();
                        tasklist.pop_front();
                    }
                    else if (! global_terminate)
                    {
                        const int iter(rand() % num_threads);
                        mc::WorkStealingThreadFunction<std::deque<mc::ThreadTask *> > * const tfunc = threads[iter];
                        if (tfunc->steal(tasklist))
                        {
                            task = tasklist.front();
                            tasklist.pop_front();
                        }
                    }
                    pthread_mutex_unlock(steal_mutex->mutex());

                    if (task == 0 && ! (global_terminate || thread_data->terminate))
                    {
                        pthread_mutex_lock(pool_mutex->mutex());
                        pthread_mutex_unlock(thread_data->local_mutex->mutex());
                        global_barrier->wait(*pool_mutex);
                        pthread_mutex_unlock(pool_mutex->mutex());
                    }
                    else if (task == 0 && thread_data->terminate)
                    {
                        pthread_mutex_unlock(thread_data->local_mutex->mutex());
                        break;
                    }
                    else
                    {
                        pthread_mutex_unlock(thread_data->local_mutex->mutex());
                    }
                }

                if (task != 0)
                {
                    unsigned & tsched_id = task->ticket.sid();
                    tsched_id = sched_id;
#ifdef DEBUG
                    std::string msg = "Thread " + stringify(pool_id) + " on LPU " +
                    stringify(sched_id) + " will execute ticket " +
                    stringify(task->ticket.uid()) + "\n";
                    LOGMESSAGE(lc_backend, msg);
#endif

                    (*task->functor)();
                    task->ticket.mark();
                    delete task;
                }
            }
            while (true);

            // The thread function's DTOR will not be called before this
            // thread has stopped execution here.
        }

        void enqueue(mc::ThreadTask * task)
        {
            Lock l(*thread_data->local_mutex);
            tasklist.push_back(task);
        }

        bool steal(std::deque<mc::ThreadTask *> & thief_list)
        {
            Lock l(*thread_data->local_mutex);

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

    template <> struct Implementation<mc::SimpleThreadFunction<CASDeque<mc::ThreadTask *> > > :
        public TFImplementationBase
    {
        /// The task list administrated by the thread pool.
        CASDeque<mc::ThreadTask *> * const tasklist;

        Implementation(mc::PoolSyncData * const psync, mc::ThreadData * const tdata, CASDeque<mc::ThreadTask *> * const list,
                unsigned pid) :
            TFImplementationBase(psync, tdata, pid),
            tasklist(list)
        {
        }

        virtual ~Implementation()
        {
        }

        void operator () ()
        {
            /// The ThreadTask to be currently executed by the thread.
            mc::ThreadTask * task(0);

            /* Set thread_id from operating system and signal the pool
             * that this thread is online. Then let the thread go to
             * sleep until it will be assigned its first task. */

            {
                Lock l(*pool_mutex);
                thread_id = syscall(__NR_gettid);
                global_barrier->broadcast();
            }

            do
            {
                task = tasklist->pop_front();

                if (task == 0)
                {
                    Lock l(*pool_mutex);
                    task = tasklist->pop_front();

                    if (task == 0)
                    {
                        if (! thread_data->terminate)
                            global_barrier->wait(*pool_mutex);
                        else
                            break;
                    }
                }

                if (task != 0)
                {
                    unsigned & sched_id = task->ticket.sid();
                    sched_id = 0xFFFF;
#ifdef DEBUG
                    std::string msg = "Thread " + stringify(pool_id) + " will execute ticket " +
                    stringify(task->ticket.uid()) + "\n";
                    LOGMESSAGE(lc_backend, msg);
#endif
                    (*task->functor)();
                    task->ticket.mark();
                    delete task;
                }
            }
            while (true);

            // The thread function's DTOR will not be called before this
            // thread has stopped execution here.
        }
    };

    template <> struct Implementation<mc::SimpleThreadFunction<ConcurrentDeque<mc::ThreadTask *> > > :
        public TFImplementationBase
    {
        /// The task list administrated by the thread pool.
        ConcurrentDeque<mc::ThreadTask *> * const tasklist;

        Implementation(mc::PoolSyncData * const psync, mc::ThreadData * const tdata, ConcurrentDeque<mc::ThreadTask *> * const list,
                unsigned pid) :
            TFImplementationBase(psync, tdata, pid),
            tasklist(list)
        {
        }

        virtual ~Implementation()
        {
        }

        void operator () ()
        {
            /// The ThreadTask to be currently executed by the thread.
            mc::ThreadTask * task(0);

            /* Set thread_id from operating system and signal the pool
             * that this thread is online. Then let the thread go to
             * sleep until it will be assigned its first task. */

            {
                Lock l(*pool_mutex);
                thread_id = syscall(__NR_gettid);
                global_barrier->broadcast();
            }

            do
            {
                task = tasklist->pop_front();

                if (task == 0)
                {
                    Lock l(*pool_mutex);
                    task = tasklist->pop_front();

                    if (task == 0)
                    {
                        if (! thread_data->terminate)
                            global_barrier->wait(*pool_mutex);
                        else
                            break;
                    }
                }

                if (task != 0)
                {
                    unsigned & sched_id = task->ticket.sid();
                    sched_id = 0xFFFF;

#ifdef DEBUG
                    std::string msg = "Thread " + stringify(pool_id) + " will execute ticket " +
                    stringify(task->ticket.uid()) + "\n";
                    LOGMESSAGE(lc_backend, msg);
#endif
                    (*task->functor)();
                    task->ticket.mark();
                    delete task;
                }
            }
            while (true);

            // The thread function's DTOR will not be called before this
            // thread has stopped execution here.
        }
    };

    template <> struct Implementation<mc::SimpleThreadFunction<std::deque<mc::ThreadTask *> > > :
        public TFImplementationBase
    {
        /// The task list administrated by the thread pool.
        std::deque<mc::ThreadTask *> * const tasklist;

        Implementation(mc::PoolSyncData * const psync, mc::ThreadData * const tdata, std::deque<mc::ThreadTask *> * const list,
                unsigned pid) :
            TFImplementationBase(psync, tdata, pid),
            tasklist(list)
        {
        }

        virtual ~Implementation()
        {
        }

        void operator () ()
        {
            /// The ThreadTask to be currently executed by the thread.
            mc::ThreadTask * task(0);

            /* Set thread_id from operating system and signal the pool
             * that this thread is online. Then let the thread go to
             * sleep until it will be assigned its first task. */

            {
                Lock l(*pool_mutex);
                thread_id = syscall(__NR_gettid);
                global_barrier->broadcast();
            }

            do
            {
                {
                    Lock l(*pool_mutex);

                    task = (tasklist->empty() ? 0 : tasklist->front());

                    if (task == 0)
                    {
                        if (! thread_data->terminate)
                            global_barrier->wait(*pool_mutex);
                        else
                            break;
                    }
                    else
                        tasklist->pop_front();
                }

                if (task != 0)
                {
                    unsigned & sched_id = task->ticket.sid();
                    sched_id = 0xFFFF;

#ifdef DEBUG
                    std::string msg = "Thread " + stringify(pool_id) + " will execute ticket " +
                    stringify(task->ticket.uid()) + "\n";
                    LOGMESSAGE(lc_backend, msg);
#endif
                    (*task->functor)();
                    task->ticket.mark();
                    delete task;
                }
            }
            while (true);

            // The thread function's DTOR will not be called before this
            // thread has stopped execution here.
        }
    };
}

using namespace honei::mc;

ThreadFunctionBase::~ThreadFunctionBase()
{
}

AffinityThreadFunction::AffinityThreadFunction(PoolSyncData * const psync, ThreadData * const tdata,
        std::deque<ThreadTask *> * const list, unsigned pool_id, LPU * const lpu) :
    PrivateImplementationPattern<AffinityThreadFunction, Shared>(new
            Implementation<AffinityThreadFunction>(psync, tdata, list, pool_id, lpu))
{
}

AffinityThreadFunction::~AffinityThreadFunction()
{
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

WorkStealingThreadFunction<std::deque<mc::ThreadTask *> >::WorkStealingThreadFunction(PoolSyncData * const psync,
        ThreadData * const tdata, unsigned pool_id, unsigned sched_id,
        const std::vector<mc::WorkStealingThreadFunction<std::deque<mc::ThreadTask *> > *> & threads,
        unsigned num_thr, volatile bool & terminate) :
    PrivateImplementationPattern<WorkStealingThreadFunction<std::deque<mc::ThreadTask *> >, Shared>(new
            Implementation<WorkStealingThreadFunction<std::deque<mc::ThreadTask *> > >(psync, tdata, pool_id, sched_id, threads, num_thr, terminate))
{
}

WorkStealingThreadFunction<std::deque<mc::ThreadTask *> >::~WorkStealingThreadFunction()
{
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

WorkStealingThreadFunction<mc::CASDeque<mc::ThreadTask *> >::WorkStealingThreadFunction(PoolSyncData * const psync,
        ThreadData * const tdata, unsigned pool_id, unsigned sched_id,
        const std::vector<mc::WorkStealingThreadFunction<mc::CASDeque<mc::ThreadTask *> > *> & threads,
        unsigned num_thr, volatile bool & terminate) :
    PrivateImplementationPattern<WorkStealingThreadFunction<mc::CASDeque<mc::ThreadTask *> >, Shared>(new
            Implementation<WorkStealingThreadFunction<mc::CASDeque<mc::ThreadTask *> > >(psync, tdata, pool_id, sched_id, threads, num_thr, terminate))
{
}

WorkStealingThreadFunction<mc::CASDeque<mc::ThreadTask *> >::~WorkStealingThreadFunction()
{
}

void WorkStealingThreadFunction<mc::CASDeque<mc::ThreadTask *> >::enqueue(mc::ThreadTask * task)
{
    _imp->enqueue(task);
}

bool WorkStealingThreadFunction<mc::CASDeque<mc::ThreadTask *> >::steal(mc::CASDeque<mc::ThreadTask *> & thief_list)
{
    return _imp->steal(thief_list);
}

void WorkStealingThreadFunction<mc::CASDeque<mc::ThreadTask *> >::operator() ()
{
    (*_imp)();
}

unsigned WorkStealingThreadFunction<mc::CASDeque<mc::ThreadTask *> >::pool_id() const
{
    return _imp->pool_id;
}

unsigned WorkStealingThreadFunction<mc::CASDeque<mc::ThreadTask *> >::tid() const
{
    return _imp->thread_id;
}

WorkStealingThreadFunction<mc::ConcurrentDeque<mc::ThreadTask *> >::WorkStealingThreadFunction(PoolSyncData * const psync,
        ThreadData * const tdata, unsigned pool_id, unsigned sched_id,
        const std::vector<mc::WorkStealingThreadFunction<mc::ConcurrentDeque<mc::ThreadTask *> > *> & threads,
        unsigned num_thr, volatile bool & terminate) :
    PrivateImplementationPattern<WorkStealingThreadFunction<mc::ConcurrentDeque<mc::ThreadTask *> >, Shared>(new
            Implementation<WorkStealingThreadFunction<mc::ConcurrentDeque<mc::ThreadTask *> > >(psync, tdata, pool_id, sched_id, threads, num_thr, terminate))
{
}

WorkStealingThreadFunction<mc::ConcurrentDeque<mc::ThreadTask *> >::~WorkStealingThreadFunction()
{
}

void WorkStealingThreadFunction<mc::ConcurrentDeque<mc::ThreadTask *> >::enqueue(mc::ThreadTask * task)
{
    _imp->enqueue(task);
}

bool WorkStealingThreadFunction<mc::ConcurrentDeque<mc::ThreadTask *> >::steal(mc::ConcurrentDeque<mc::ThreadTask *> & thief_list)
{
    return _imp->steal(thief_list);
}

void WorkStealingThreadFunction<mc::ConcurrentDeque<mc::ThreadTask *> >::operator() ()
{
    (*_imp)();
}

unsigned WorkStealingThreadFunction<mc::ConcurrentDeque<mc::ThreadTask *> >::pool_id() const
{
    return _imp->pool_id;
}

unsigned WorkStealingThreadFunction<mc::ConcurrentDeque<mc::ThreadTask *> >::tid() const
{
    return _imp->thread_id;
}


SimpleThreadFunction<ConcurrentDeque<ThreadTask *> >::SimpleThreadFunction(PoolSyncData * const psync, ThreadData * const tdata,
        ConcurrentDeque<ThreadTask *> * const list, unsigned pool_id) :
    PrivateImplementationPattern<SimpleThreadFunction<ConcurrentDeque<ThreadTask *> >, Shared>(new
            Implementation<SimpleThreadFunction<ConcurrentDeque<ThreadTask *> > >(psync, tdata, list, pool_id))
{
}

SimpleThreadFunction<ConcurrentDeque<ThreadTask *> >::~SimpleThreadFunction()
{
}

void SimpleThreadFunction<ConcurrentDeque<ThreadTask *> >::operator() ()
{
    (*_imp)();
}

unsigned SimpleThreadFunction<ConcurrentDeque<ThreadTask *> >::tid() const
{
    return _imp->thread_id;
}

SimpleThreadFunction<std::deque<ThreadTask *> >::SimpleThreadFunction(PoolSyncData * const psync, ThreadData * const tdata,
        std::deque<ThreadTask *> * const list, unsigned pool_id) :
    PrivateImplementationPattern<SimpleThreadFunction<std::deque<ThreadTask *> >, Shared>(new
            Implementation<SimpleThreadFunction<std::deque<ThreadTask *> > >(psync, tdata, list, pool_id))
{
}

SimpleThreadFunction<std::deque<ThreadTask *> >::~SimpleThreadFunction()
{
}

void SimpleThreadFunction<std::deque<ThreadTask *> >::operator() ()
{
    (*_imp)();
}

unsigned SimpleThreadFunction<std::deque<ThreadTask *> >::tid() const
{
    return _imp->thread_id;
}

SimpleThreadFunction<CASDeque<ThreadTask *> >::SimpleThreadFunction(PoolSyncData * const psync, ThreadData * const tdata,
        CASDeque<ThreadTask *> * const list, unsigned pool_id) :
    PrivateImplementationPattern<SimpleThreadFunction<CASDeque<ThreadTask *> >, Shared>(new
            Implementation<SimpleThreadFunction<CASDeque<ThreadTask *> > >(psync, tdata, list, pool_id))
{
}

SimpleThreadFunction<CASDeque<ThreadTask *> >::~SimpleThreadFunction()
{
}

void SimpleThreadFunction<CASDeque<ThreadTask *> >::operator() ()
{
    (*_imp)();
}

unsigned SimpleThreadFunction<CASDeque<ThreadTask *> >::tid() const
{
    return _imp->thread_id;
}
