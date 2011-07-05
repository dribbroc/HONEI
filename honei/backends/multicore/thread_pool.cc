/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008, 2009, 2010, 2011 Sven Mallach <mallach@honei.org>
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

#include <honei/backends/multicore/thread_pool.hh>
#include <honei/backends/multicore/thread_function.hh>
#include <honei/backends/multicore/topology.hh>
#include <honei/util/configuration.hh>
#include <honei/util/exception.hh>
#include <honei/util/instantiation_policy-impl.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <honei/util/lock.hh>
#include <honei/util/log.hh>
#include <honei/util/stringify.hh>
#include <honei/util/thread.hh>

#include <errno.h>
#include <sched.h>
#include <stdlib.h>
#include <sys/syscall.h>

namespace honei
{
    template <> struct Implementation<mc::ThreadPool>
    {
        /// Information about processor topology (such as number of processing units)
        mc::Topology * const _topology;

        /// Number of demanded threads (via honeirc)
        const unsigned _demanded_threads;

        /// Number of pooled threads
        const unsigned _num_threads;

        /// List of user POSIX threads
        std::list<Thread *> _threads;

        /// List of thread specific data
        std::list<mc::ThreadData *> _thread_data;

        /// Exchange data between the pool and its threads
        mc::PoolSyncData * const _pool_sync;

        /// Function pointer to any of the default dispatch strategies
        mc::DispatchPolicy (* policy) ();

        Implementation() :
            _topology(mc::Topology::instance()),
            _demanded_threads(Configuration::instance()->get_value("mc::num_threads", _topology->num_lpus())),
            _num_threads(_demanded_threads > _topology->num_lpus() ? _demanded_threads : _topology->num_lpus()),
            _pool_sync(new mc::PoolSyncData)
        {
            CONTEXT("When initializing the basic implementation of thread pool:");

            bool affinity = Configuration::instance()->get_value("mc::affinity", true);

#ifdef DEBUG
            std::string msg = "Affinity is ";
            msg += (affinity ? "enabled.\n" : "disabled.\n");
            msg += "The default DispatchPolicy is configured to be: \n";
#endif
            std::string dis = Configuration::instance()->get_value("mc::dispatch", "anycore");

            if (! affinity || dis == "anycore")
            {
                policy = &mc::DispatchPolicy::any_core;
#ifdef DEBUG
                msg += "arbitrary - the next available thread will execute a task\n";
#endif
            }
            else if (dis == "alternating")
            {
                policy = &mc::DispatchPolicy::alternating_node;
#ifdef DEBUG
                msg += "alternating - tasks shall be assigned alternatingly on available NUMA nodes\n";
#endif
            }
            else if (dis == "linear")
            {
                policy = &mc::DispatchPolicy::linear_node;
#ifdef DEBUG
                msg += "linear - tasks shall be assigned to fill up available NUMA nodes one-by-one\n";
#endif
            }

#ifdef DEBUG
            msg += "Will create " + stringify(_num_threads) + " POSIX worker threads\n";
            LOGMESSAGE(lc_backend, msg);
#endif
        }

        virtual ~Implementation()
        {
            std::list<mc::ThreadData *>::iterator j(_thread_data.begin());
            for(std::list<Thread *>::iterator i(_threads.begin()), i_end(_threads.end()) ; i != i_end ; ++i, ++j)
            {
                {
                    Lock l(*_pool_sync->mutex);
                    (*j)->terminate = true;
                    _pool_sync->barrier->broadcast();
                }

                delete (*i);
                delete (*j);
            }

            delete _pool_sync;
        }

        virtual Ticket<tags::CPU::MultiCore> enqueue(const function<void ()> & task, mc::DispatchPolicy p) = 0;
        virtual Ticket<tags::CPU::MultiCore> enqueue(const function<void ()> & task) = 0;
    };

    /* StandardImplementation assumes affinity to be disabled. */

    template <typename ListType> struct StandardImplementation;

    template <> struct StandardImplementation<mc::CASDeque<mc::ThreadTask *> > :
        public Implementation<mc::ThreadPool>
    {
        /// Waiting list of worker tasks to be executed (otherwise)
        mc::CASDeque<mc::ThreadTask *> _tasks;

        StandardImplementation() :
            Implementation<mc::ThreadPool>()
        {
            CONTEXT("When initializing the standard implementation of thread pool:");

            int inst_ctr(0);

            for (int i(_num_threads - 1) ; i >= 0 ; --i)
            {
                mc::ThreadData * td = new mc::ThreadData;
                mc::SimpleThreadFunction<mc::CASDeque<mc::ThreadTask *> > tobj(_pool_sync, td, &_tasks, inst_ctr);

                Thread * t;
                {
                    Lock l(*_pool_sync->mutex);
                    t = new Thread(tobj); // tobj will be copied here!
                    _pool_sync->barrier->wait(*_pool_sync->mutex); // Wait until the thread is really setup / got cpu time for the first time
                }

                _threads.push_back(t);
                _thread_data.push_back(td);

                ++inst_ctr;
            }
        }

        ~StandardImplementation()
        {
        }

        virtual Ticket<tags::CPU::MultiCore> enqueue(const function<void ()> & task, mc::DispatchPolicy)
        {
            return enqueue(task);
        }

        virtual Ticket<tags::CPU::MultiCore> enqueue(const function<void ()> & task)
        {
            CONTEXT("When creating a ThreadTask:");

            Ticket<tags::CPU::MultiCore> ticket(policy().apply());
            mc::ThreadTask * t_task(new mc::ThreadTask(task, ticket));

            _tasks.push_back(t_task);

            {
                Lock l(*_pool_sync->mutex);
                _pool_sync->barrier->broadcast();
            }

            return ticket;
        }
    };

    template <> struct StandardImplementation<mc::ConcurrentDeque<mc::ThreadTask *> > :
        public Implementation<mc::ThreadPool>
    {
        /// Waiting list of worker tasks to be executed (otherwise)
        mc::ConcurrentDeque<mc::ThreadTask *> _tasks;

        StandardImplementation() :
            Implementation<mc::ThreadPool>()
        {
            CONTEXT("When initializing the standard implementation of thread pool:");

            int inst_ctr(0);

            for (int i(_num_threads - 1) ; i >= 0 ; --i)
            {
                mc::ThreadData * td = new mc::ThreadData;
                mc::SimpleThreadFunction<mc::ConcurrentDeque<mc::ThreadTask *> > tobj(_pool_sync, td, &_tasks, inst_ctr);

                Thread * t;
                {
                    Lock l(*_pool_sync->mutex);
                    t = new Thread(tobj); // tobj will be copied here!
                    _pool_sync->barrier->wait(*_pool_sync->mutex); // Wait until the thread is really setup / got cpu time for the first time
                }

                _threads.push_back(t);
                _thread_data.push_back(td);

                ++inst_ctr;
            }
        }

        ~StandardImplementation()
        {
        }

        virtual Ticket<tags::CPU::MultiCore> enqueue(const function<void ()> & task, mc::DispatchPolicy)
        {
            return enqueue(task);
        }

        virtual Ticket<tags::CPU::MultiCore> enqueue(const function<void ()> & task)
        {
            CONTEXT("When creating a ThreadTask:");

            Ticket<tags::CPU::MultiCore> ticket(policy().apply());
            mc::ThreadTask * t_task(new mc::ThreadTask(task, ticket));

            _tasks.push_back(t_task);

            {
                Lock l(*_pool_sync->mutex);
                _pool_sync->barrier->broadcast();
            }

            return ticket;
        }
    };

    template <> struct StandardImplementation<std::deque<mc::ThreadTask *> > :
        public Implementation<mc::ThreadPool>
    {
        /// Waiting list of worker tasks to be executed (otherwise)
        std::deque<mc::ThreadTask *> _tasks;

        StandardImplementation() :
            Implementation<mc::ThreadPool>()
        {
            CONTEXT("When initializing the standard implementation of thread pool:");

            int inst_ctr(0);

            for (int i(_num_threads - 1) ; i >= 0 ; --i)
            {
                mc::ThreadData * td = new mc::ThreadData;
                mc::SimpleThreadFunction<std::deque<mc::ThreadTask *> > tobj(_pool_sync, td, &_tasks, inst_ctr);
                Thread * t;
                {
                    Lock l(*_pool_sync->mutex);
                    t = new Thread(tobj); // tobj will be copied here!
                    _pool_sync->barrier->wait(*_pool_sync->mutex); // Wait until the thread is really setup / got cpu time for the first time
                }

                _threads.push_back(t);
                _thread_data.push_back(td);

                ++inst_ctr;
            }
        }

        ~StandardImplementation()
        {
        }

        virtual Ticket<tags::CPU::MultiCore> enqueue(const function<void ()> & task, mc::DispatchPolicy)
        {
            return enqueue(task);
        }

        virtual Ticket<tags::CPU::MultiCore> enqueue(const function<void ()> & task)
        {
            CONTEXT("When creating a ThreadTask:");

            Ticket<tags::CPU::MultiCore> ticket(policy().apply());
            mc::ThreadTask * t_task(new mc::ThreadTask(task, ticket));

            {
                Lock l(*_pool_sync->mutex);
                _tasks.push_back(t_task);
                _pool_sync->barrier->broadcast();
            }

            return ticket;
        }
    };

    struct AffinityImplementation :
        public Implementation<mc::ThreadPool>
    {
        /// Waiting list of worker tasks to be executed (if affinity is enabled)
        std::deque<mc::ThreadTask *> _tasks;

        /// Mapping of threads to the scheduler ids of the cores they run on
        std::vector<unsigned> _sched_ids;

        /// Array of affinity masks for main process and all controlled threads
        cpu_set_t * _affinity_mask;

        AffinityImplementation() :
            Implementation<mc::ThreadPool>()
        {
            CONTEXT("When initializing the affinity implementation of thread pool:");

            _affinity_mask = new cpu_set_t[_num_threads + 1];

            // set own affinity first
            CPU_ZERO(&_affinity_mask[_num_threads]);
            CPU_SET(_topology->num_lpus() - 1, &_affinity_mask[_num_threads]);
            if(sched_setaffinity(syscall(__NR_gettid), sizeof(cpu_set_t), &_affinity_mask[_num_threads]) != 0)
                throw ExternalError("Unix: sched_setaffinity()", "could not set affinity! errno: " + stringify(errno));

#ifdef DEBUG
            std::string msg = "THREAD \t\t POOL_ID \t LPU \t NODE \n";
            msg += "MAIN \t\t - \t\t" + stringify(_topology->num_lpus() - 1) + "\t\t" + stringify(_topology->get_node(_topology->num_lpus() - 1)) + " \n";
#endif

            int inst_ctr(0);

            for (int i(_num_threads - 1) ; i >= 0 ; --i)
            {
                unsigned sched_id(i % (_topology->num_lpus()));
                mc::ThreadData * td = new mc::ThreadData;

                mc::AffinityThreadFunction tobj(_pool_sync, td, &_tasks, inst_ctr, sched_id);

                Thread * t;
                {
                    Lock l(*_pool_sync->mutex);
                    t = new Thread(tobj); // tobj will be copied here!
                    _pool_sync->barrier->wait(*_pool_sync->mutex); // Wait until the thread is really setup / got cpu time for the first time
                }

                _threads.push_back(t);
                _thread_data.push_back(td);
                _sched_ids.push_back(sched_id);
                CPU_ZERO(&_affinity_mask[i]);
                CPU_SET(sched_id, &_affinity_mask[i]);
                if(sched_setaffinity(tobj.tid(), sizeof(cpu_set_t), &_affinity_mask[i]) != 0)
                    throw ExternalError("Unix: sched_setaffinity()", "could not set affinity! errno: " + stringify(errno));
#ifdef DEBUG
                msg += stringify(tobj.tid()) + "\t\t" + stringify(inst_ctr) + "\t\t" + stringify(sched_id) + "\t\t" + stringify(_topology->get_node(sched_id)) + " \n";
#endif

                ++inst_ctr;
            }

#ifdef DEBUG
            LOGMESSAGE(lc_backend, msg);
#endif
        }

        ~AffinityImplementation()
        {
            delete[] _affinity_mask;
        }

        virtual Ticket<tags::CPU::MultiCore> enqueue(const function<void ()> & task, mc::DispatchPolicy p)
        {
            CONTEXT("When creating a ThreadTask:");

            Ticket<tags::CPU::MultiCore> ticket(p.apply());
            mc::ThreadTask * t_task(new mc::ThreadTask(task, ticket));

            {
                Lock l(*_pool_sync->mutex);
                _tasks.push_back(t_task);
                _pool_sync->barrier->broadcast();
            }

            return ticket;
        }

        virtual Ticket<tags::CPU::MultiCore> enqueue(const function<void ()> & task)
        {
            CONTEXT("When creating a ThreadTask:");

            Ticket<tags::CPU::MultiCore> ticket(policy().apply());
            mc::ThreadTask * t_task(new mc::ThreadTask(task, ticket));

            {
                Lock l(*_pool_sync->mutex);
                _tasks.push_back(t_task);
                _pool_sync->barrier->broadcast();
            }

            return ticket;
        }
    };

    template <typename ListType> struct WorkStealingImplementation;

    template <> struct WorkStealingImplementation<std::deque<mc::ThreadTask *> > :
        public Implementation<mc::ThreadPool>
    {
        /// List of user POSIX threads
        std::vector<mc::WorkStealingThreadFunction<std::deque<mc::ThreadTask *> > *> _thread_fn;

        /// Mapping of threads to the scheduler ids of the cores they run on
        std::vector<unsigned> _sched_ids;

        /// Array of affinity masks for main process and all controlled threads
        cpu_set_t * _affinity_mask;

        volatile bool global_terminate;

        WorkStealingImplementation() :
            Implementation<mc::ThreadPool>(),
            _affinity_mask(NULL),
            global_terminate(false)
        {
            CONTEXT("When initializing the work stealing implementation thread pool:");

            bool affinity = Configuration::instance()->get_value("mc::affinity", true);

#ifdef DEBUG
            std::string msg;
#endif

            if (affinity)
            {
                _affinity_mask = new cpu_set_t[_num_threads + 1];

                // set own affinity first
                CPU_ZERO(&_affinity_mask[_num_threads]);
                CPU_SET(_topology->num_lpus() - 1, &_affinity_mask[_num_threads]);
                if(sched_setaffinity(syscall(__NR_gettid), sizeof(cpu_set_t), &_affinity_mask[_num_threads]) != 0)
                    throw ExternalError("Unix: sched_setaffinity()", "could not set affinity! errno: " + stringify(errno));

#ifdef DEBUG
                msg += "THREAD \t\t POOL_ID \t LPU \t NODE \n";
                msg += "MAIN \t\t - \t\t" + stringify(_topology->num_lpus() - 1) + "\t\t" + stringify(_topology->get_node(_topology->num_lpus() - 1)) + " \n";
#endif
            }

            int inst_ctr(0);

            Lock l(*_pool_sync->steal_mutex); // Prevent threads from stealing before all threads are alive
            for (int i(_num_threads - 1) ; i >= 0 ; --i)
            {
                unsigned sched_id(affinity ? (i % (_topology->num_lpus())) : 0xFFFF );
                mc::ThreadData * td = new mc::ThreadData;

                mc::WorkStealingThreadFunction<std::deque<mc::ThreadTask *> > * tobj =
                    new mc::WorkStealingThreadFunction<std::deque<mc::ThreadTask *> >(_pool_sync, td, inst_ctr, sched_id,
                            _thread_fn, _num_threads, global_terminate);

                Thread * t;
                {
                    Lock l(*_pool_sync->mutex);
                    t = new Thread(*tobj); // tobj will be copied here!
                    _pool_sync->barrier->wait(*_pool_sync->mutex); // Wait until the thread is really setup / got cpu time for the first time
                }

                _threads.push_back(t);
                _thread_data.push_back(td);
                _thread_fn.push_back(tobj);

                if (affinity)
                {
                    _sched_ids.push_back(sched_id);
                    CPU_ZERO(&_affinity_mask[i]);
                    CPU_SET(sched_id, &_affinity_mask[i]);
                    if(sched_setaffinity(tobj->tid(), sizeof(cpu_set_t), &_affinity_mask[i]) != 0)
                        throw ExternalError("Unix: sched_setaffinity()", "could not set affinity! errno: " + stringify(errno));
#ifdef DEBUG
                    msg += stringify(tobj->tid()) + "\t\t" + stringify(inst_ctr) + "\t\t" + stringify(sched_id) + "\t\t" + stringify(_topology->get_node(sched_id)) + " \n";
#endif
                }

                ++inst_ctr;
            }
            // Now let all threads go and run

#ifdef DEBUG
            LOGMESSAGE(lc_backend, msg);
#endif
        }

        ~WorkStealingImplementation()
        {
            {
                Lock l(*_pool_sync->steal_mutex);
                global_terminate = true;
            }

            // This will delete our copies - Thread has an own one!
            for (std::vector<mc::WorkStealingThreadFunction<std::deque<mc::ThreadTask *> > *>::iterator i(_thread_fn.begin()),
                    i_end(_thread_fn.end()) ; i != i_end ; ++i)
            {
                delete (*i);
            }

            delete[] _affinity_mask;
        }

        virtual Ticket<tags::CPU::MultiCore> enqueue(const function<void ()> & task, mc::DispatchPolicy p)
        {
            CONTEXT("When creating a ThreadTask:");

            Ticket<tags::CPU::MultiCore> ticket(p.apply());
            mc::ThreadTask * t_task(new mc::ThreadTask(task, ticket));

            int idx((ticket.sid_min() == 0xFFFF) ? rand() % _num_threads : ticket.sid_min());

            mc::WorkStealingThreadFunction<std::deque<mc::ThreadTask *> > * wfunc(_thread_fn[idx]);
            wfunc->enqueue(t_task);

            {
                Lock l(*_pool_sync->mutex);
                _pool_sync->barrier->broadcast();
            }

            return ticket;
        }

        virtual Ticket<tags::CPU::MultiCore> enqueue(const function<void ()> & task)
        {
            CONTEXT("When creating a ThreadTask:");

            Ticket<tags::CPU::MultiCore> ticket(policy().apply());
            mc::ThreadTask * t_task(new mc::ThreadTask(task, ticket));

            int idx((ticket.sid_min() == 0xFFFF) ? rand() % _num_threads : ticket.sid_min());

            mc::WorkStealingThreadFunction<std::deque<mc::ThreadTask *> > * wfunc(_thread_fn[idx]);
            wfunc->enqueue(t_task);

            {
                Lock l(*_pool_sync->mutex);
                _pool_sync->barrier->broadcast();
            }

            return ticket;
        }
    };

    template <> struct WorkStealingImplementation<mc::ConcurrentDeque<mc::ThreadTask *> > :
        public Implementation<mc::ThreadPool>
    {
        /// List of user POSIX threads
        std::vector<mc::WorkStealingThreadFunction<mc::ConcurrentDeque<mc::ThreadTask *> > *> _thread_fn;

        /// Mapping of threads to the scheduler ids of the cores they run on
        std::vector<unsigned> _sched_ids;

        /// Array of affinity masks for main process and all controlled threads
        cpu_set_t * _affinity_mask;

        volatile bool global_terminate;

        WorkStealingImplementation() :
            Implementation<mc::ThreadPool>(),
            _affinity_mask(NULL),
            global_terminate(false)
        {
            CONTEXT("When initializing the work stealing implementation thread pool:");

            bool affinity = Configuration::instance()->get_value("mc::affinity", true);

#ifdef DEBUG
            std::string msg;
#endif

            if (affinity)
            {
                _affinity_mask = new cpu_set_t[_num_threads + 1];

                // set own affinity first
                CPU_ZERO(&_affinity_mask[_num_threads]);
                CPU_SET(_topology->num_lpus() - 1, &_affinity_mask[_num_threads]);
                if(sched_setaffinity(syscall(__NR_gettid), sizeof(cpu_set_t), &_affinity_mask[_num_threads]) != 0)
                    throw ExternalError("Unix: sched_setaffinity()", "could not set affinity! errno: " + stringify(errno));

#ifdef DEBUG
                msg += "THREAD \t\t POOL_ID \t LPU \t NODE \n";
                msg += "MAIN \t\t - \t\t" + stringify(_topology->num_lpus() - 1) + "\t\t" + stringify(_topology->get_node(_topology->num_lpus() - 1)) + " \n";
#endif
            }

            int inst_ctr(0);
            Lock l(*_pool_sync->steal_mutex); // Prevent threads from stealing before all threads are alive
            for (int i(_num_threads - 1) ; i >= 0 ; --i)
            {
                unsigned sched_id(affinity ? (i % (_topology->num_lpus())) : 0xFFFF);

                mc::ThreadData * td = new mc::ThreadData;

                mc::WorkStealingThreadFunction<mc::ConcurrentDeque<mc::ThreadTask *> > * tobj =
                    new mc::WorkStealingThreadFunction<mc::ConcurrentDeque<mc::ThreadTask *> >(_pool_sync, td, inst_ctr, sched_id,
                            _thread_fn, _num_threads, global_terminate);

                Thread * t;
                {
                    Lock l(*_pool_sync->mutex);
                    t = new Thread(*tobj); // tobj will be copied here!
                    _pool_sync->barrier->wait(*_pool_sync->mutex); // Wait until the thread is really setup / got cpu time for the first time
                }

                _threads.push_back(t);
                _thread_data.push_back(td);
                _thread_fn.push_back(tobj);

                if (affinity)
                {
                    _sched_ids.push_back(sched_id);
                    CPU_ZERO(&_affinity_mask[i]);
                    CPU_SET(sched_id, &_affinity_mask[i]);
                    if(sched_setaffinity(tobj->tid(), sizeof(cpu_set_t), &_affinity_mask[i]) != 0)
                        throw ExternalError("Unix: sched_setaffinity()", "could not set affinity! errno: " + stringify(errno));
#ifdef DEBUG
                    msg += stringify(tobj->tid()) + "\t\t" + stringify(inst_ctr) + "\t\t" + stringify(sched_id) + "\t\t" + stringify(_topology->get_node(sched_id)) + " \n";
#endif
                }

                ++inst_ctr;
            }

            // Now let all threads go and run

#ifdef DEBUG
            LOGMESSAGE(lc_backend, msg);
#endif
        }

        ~WorkStealingImplementation()
        {
            {
                Lock l(*_pool_sync->steal_mutex);
                global_terminate = true;
            }

            // This will delete our copies - Thread has an own one!
            for (std::vector<mc::WorkStealingThreadFunction<mc::ConcurrentDeque<mc::ThreadTask *> > *>::iterator i(_thread_fn.begin()),
                    i_end(_thread_fn.end()) ; i != i_end ; ++i)
            {
                delete (*i);
            }

            delete[] _affinity_mask;
        }

        virtual Ticket<tags::CPU::MultiCore> enqueue(const function<void ()> & task, mc::DispatchPolicy p)
        {
            CONTEXT("When creating a ThreadTask:");

            Ticket<tags::CPU::MultiCore> ticket(p.apply());
            mc::ThreadTask * t_task(new mc::ThreadTask(task, ticket));

            int idx((ticket.sid_min() == 0xFFFF) ? rand() % _num_threads : ticket.sid_min());

            mc::WorkStealingThreadFunction<mc::ConcurrentDeque<mc::ThreadTask *> > * wfunc(_thread_fn[idx]);
            wfunc->enqueue(t_task);

            {
                Lock l(*_pool_sync->mutex);
                _pool_sync->barrier->broadcast();
            }

            return ticket;
        }

        virtual Ticket<tags::CPU::MultiCore> enqueue(const function<void ()> & task)
        {
            CONTEXT("When creating a ThreadTask:");

            Ticket<tags::CPU::MultiCore> ticket(policy().apply());
            mc::ThreadTask * t_task(new mc::ThreadTask(task, ticket));

            int idx((ticket.sid_min() == 0xFFFF) ? rand() % _num_threads : ticket.sid_min());

            mc::WorkStealingThreadFunction<mc::ConcurrentDeque<mc::ThreadTask *> > * wfunc(_thread_fn[idx]);
            wfunc->enqueue(t_task);

            {
                Lock l(*_pool_sync->mutex);
                _pool_sync->barrier->broadcast();
            }

            return ticket;
        }
    };

    template <> struct WorkStealingImplementation<mc::CASDeque<mc::ThreadTask *> > :
        public Implementation<mc::ThreadPool>
    {
        /// List of user POSIX threads
        std::vector<mc::WorkStealingThreadFunction<mc::CASDeque<mc::ThreadTask *> > *> _thread_fn;

        /// Mapping of threads to the scheduler ids of the cores they run on
        std::vector<unsigned> _sched_ids;

        /// Array of affinity masks for main process and all controlled threads
        cpu_set_t * _affinity_mask;

        volatile bool global_terminate;

        WorkStealingImplementation() :
            Implementation<mc::ThreadPool>(),
            _affinity_mask(NULL),
            global_terminate(false)
        {
            CONTEXT("When initializing the work stealing implementation thread pool:");

            bool affinity = Configuration::instance()->get_value("mc::affinity", true);

#ifdef DEBUG
            std::string msg;
#endif

            if (affinity)
            {
                _affinity_mask = new cpu_set_t[_num_threads + 1];

                // set own affinity first
                CPU_ZERO(&_affinity_mask[_num_threads]);
                CPU_SET(_topology->num_lpus() - 1, &_affinity_mask[_num_threads]);
                if(sched_setaffinity(syscall(__NR_gettid), sizeof(cpu_set_t), &_affinity_mask[_num_threads]) != 0)
                    throw ExternalError("Unix: sched_setaffinity()", "could not set affinity! errno: " + stringify(errno));

#ifdef DEBUG
                msg += "THREAD \t\t POOL_ID \t LPU \t NODE \n";
                msg += "MAIN \t\t - \t\t" + stringify(_topology->num_lpus() - 1) + "\t\t" + stringify(_topology->get_node(_topology->num_lpus() - 1)) + " \n";
#endif
            }

            int inst_ctr(0);
            Lock l(*_pool_sync->steal_mutex); // Prevent threads from stealing before all threads are alive
            for (int i(_num_threads - 1) ; i >= 0 ; --i)
            {
                unsigned sched_id(affinity ? (i % (_topology->num_lpus())) : 0xFFFF);

                mc::ThreadData * td = new mc::ThreadData;

                mc::WorkStealingThreadFunction<mc::CASDeque<mc::ThreadTask *> > * tobj =
                    new mc::WorkStealingThreadFunction<mc::CASDeque<mc::ThreadTask *> >(_pool_sync, td, inst_ctr, sched_id,
                            _thread_fn, _num_threads, global_terminate);

                Thread * t;
                {
                    Lock l(*_pool_sync->mutex);
                    t = new Thread(*tobj); // tobj will be copied here!
                    _pool_sync->barrier->wait(*_pool_sync->mutex); // Wait until the thread is really setup / got cpu time for the first time
                }

                _threads.push_back(t);
                _thread_data.push_back(td);
                _thread_fn.push_back(tobj);

                if (affinity)
                {
                    _sched_ids.push_back(sched_id);
                    CPU_ZERO(&_affinity_mask[i]);
                    CPU_SET(sched_id, &_affinity_mask[i]);
                    if(sched_setaffinity(tobj->tid(), sizeof(cpu_set_t), &_affinity_mask[i]) != 0)
                        throw ExternalError("Unix: sched_setaffinity()", "could not set affinity! errno: " + stringify(errno));
#ifdef DEBUG
                    msg += stringify(tobj->tid()) + "\t\t" + stringify(inst_ctr) + "\t\t" + stringify(sched_id) + "\t\t" + stringify(_topology->get_node(sched_id)) + " \n";
#endif
                }

                ++inst_ctr;
            }

            // Now let all threads go and run

#ifdef DEBUG
            LOGMESSAGE(lc_backend, msg);
#endif
        }

        ~WorkStealingImplementation()
        {
            {
                Lock l(*_pool_sync->steal_mutex);
                global_terminate = true;
            }

            // This will delete our copies - Thread has an own one!
            for (std::vector<mc::WorkStealingThreadFunction<mc::CASDeque<mc::ThreadTask *> > *>::iterator i(_thread_fn.begin()),
                    i_end(_thread_fn.end()) ; i != i_end ; ++i)
            {
                delete (*i);
            }

            delete[] _affinity_mask;
        }

        virtual Ticket<tags::CPU::MultiCore> enqueue(const function<void ()> & task, mc::DispatchPolicy p)
        {
            CONTEXT("When creating a ThreadTask:");

            Ticket<tags::CPU::MultiCore> ticket(p.apply());
            mc::ThreadTask * t_task(new mc::ThreadTask(task, ticket));

            int idx((ticket.sid_min() == 0xFFFF) ? rand() % _num_threads : ticket.sid_min());

            mc::WorkStealingThreadFunction<mc::CASDeque<mc::ThreadTask *> > * wfunc(_thread_fn[idx]);
            wfunc->enqueue(t_task);

            {
                Lock l(*_pool_sync->mutex);
                _pool_sync->barrier->broadcast();
            }

            return ticket;
        }

        virtual Ticket<tags::CPU::MultiCore> enqueue(const function<void ()> & task)
        {
            CONTEXT("When creating a ThreadTask:");

            Ticket<tags::CPU::MultiCore> ticket(policy().apply());
            mc::ThreadTask * t_task(new mc::ThreadTask(task, ticket));

            int idx((ticket.sid_min() == 0xFFFF) ? rand() % _num_threads : ticket.sid_min());

            mc::WorkStealingThreadFunction<mc::CASDeque<mc::ThreadTask *> > * wfunc(_thread_fn[idx]);
            wfunc->enqueue(t_task);

            {
                Lock l(*_pool_sync->mutex);
                _pool_sync->barrier->broadcast();
            }

            return ticket;
        }
    };
}

using namespace honei;
using namespace honei::mc;

Ticket<tags::CPU::MultiCore> DispatchPolicy::last;

template class InstantiationPolicy<ThreadPool, Singleton>;

ThreadPool::ThreadPool() :
    PrivateImplementationPattern<ThreadPool, Shared>(select_impl())
{
}

ThreadPool::~ThreadPool()
{
}

Implementation<ThreadPool> * ThreadPool::select_impl()
{
    bool works = Configuration::instance()->get_value("mc::work_stealing", false);

    bool affinity = Configuration::instance()->get_value("mc::affinity", true);

    int listtype = Configuration::instance()->get_value("mc::listtype", 0);

#ifndef linux
    affinity = false;
#endif

    if (works)
    {
        if (listtype == 2)
        {
            return new WorkStealingImplementation<mc::CASDeque<mc::ThreadTask *> >;
        }
        else if (listtype == 1)
        {
            return new WorkStealingImplementation<mc::ConcurrentDeque<mc::ThreadTask *> >;
        }
        else
        {
            return new WorkStealingImplementation<std::deque<mc::ThreadTask *> >;
        }

    }
    else if (affinity)
    {
        return new AffinityImplementation;
    }
    else
    {
        if (listtype == 2)
        {
            return new StandardImplementation<mc::CASDeque<mc::ThreadTask *> >;
        }
        else if (listtype == 1)
        {
            return new StandardImplementation<mc::ConcurrentDeque<mc::ThreadTask *> >;
        }
        else
        {
            return new StandardImplementation<std::deque<mc::ThreadTask *> >;
        }
    }
}

unsigned ThreadPool::num_threads() const
{
    return _imp->_num_threads;
}

Ticket<tags::CPU::MultiCore> ThreadPool::enqueue(const function<void ()> & task, DispatchPolicy p)
{
    return _imp->enqueue(task, p);
}

/// Use default policy
Ticket<tags::CPU::MultiCore> ThreadPool::enqueue(const function<void ()> & task)
{
    return _imp->enqueue(task);
}
