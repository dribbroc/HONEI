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

#pragma once
#ifndef MULTICORE_GUARD_THREAD_POOL_HH
#define MULTICORE_GUARD_THREAD_POOL_HH 1

#include <honei/backends/multicore/dispatch_policy.hh>
#include <honei/backends/multicore/ticket.hh>
#include <honei/backends/multicore/thread_function.hh>
#include <honei/backends/multicore/topology.hh>
#include <honei/util/attributes.hh>
#include <honei/util/instantiation_policy.hh>
#include <honei/util/thread.hh>

#include <vector>

namespace honei
{
    namespace mc
    {
        class ThreadPool :
            public InstantiationPolicy<ThreadPool, Singleton>
        {
            private:
                /// \name Private members
                /// \{

                /// Information about processor topology (such as number of processing units)
                Topology * const _topology;

                /// Number of currently pooled threads
                unsigned _num_threads;

                /// A thread instantiation counter
                unsigned _inst_ctr;

                /// List of user POSIX threads
                std::list<std::pair<Thread *, ThreadFunction *> > _threads;

                /// Waiting list of worker tasks to be executed
                std::list<ThreadTask *> _tasks;

                /// Our Mutex
                Mutex * const _mutex;

                /// Condition Variable used to synchronize all threads
                ConditionVariable * const _global_barrier;

                /// Flag whether to use thread affinity
                const bool _affinity;

                /// Function pointer to any of the default dispatch strategies
                DispatchPolicy (* policy) ();

#ifdef linux
                /// Mapping of threads to the scheduler ids of the cores they run on
                std::vector<unsigned> _sched_ids;

                /// Array of affinity masks for main process and all controlled threads
                cpu_set_t * _affinity_mask;
#endif
                /// \}

            protected:

                friend class InstantiationPolicy<ThreadPool, Singleton>;

                /// \name Basic Operations
                /// \{

                /// Constructor

                ThreadPool();

                /// \}

            public:

                /// \name Basic Operations
                /// \{

                /// Destructor

                ~ThreadPool();

                /// \}

                /// \name Public members
                /// \{

                /// Add threads to the pool
                void add_threads(const unsigned num);

                /// Remove threads from the pool
                void delete_threads(const unsigned num);

                /// Retrieve the number of NUMA nodes
                unsigned num_nodes() const;

                /// Retrieve the node on which the main thread runs (use only when affinity enabled)
                unsigned main_node() const;

                /// Retrieve the number of created threads
                unsigned num_threads() const;

                /// Dispatch a task using a specified dispatch strategy
                Ticket<tags::CPU::MultiCore> * enqueue(const function<void ()> & task, DispatchPolicy p);

                /// Dispatch a task using the default dispatch strategy setup through honeirc
                Ticket<tags::CPU::MultiCore> * enqueue(const function<void ()> & task);

                /// \}
        };
    }
}
#endif
