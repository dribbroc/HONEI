/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008, 2009 Sven Mallach <sven.mallach@cs.tu-dortmund.de>
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

#ifndef MULTICORE_GUARD_THREAD_POOL_HH
#define MULTICORE_GUARD_THREAD_POOL_HH 1

#include <honei/backends/multicore/dispatch_policy.hh>
#include <honei/backends/multicore/ticket.hh>
#include <honei/backends/multicore/thread_function.hh>
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

                // Number of available logical processing units
                const unsigned num_lpus;

                // Number of threads use
                unsigned num_threads;

                // List of user POSIX threads
                std::list<std::pair<Thread *, ThreadFunction *> > threads HONEI_ALIGNED(128);

                // Waiting list of worker tasks to be executed
                std::list<ThreadTask *> tasks HONEI_ALIGNED(128);

                // Our Mutex
                Mutex * const mutex;

                // Condition Variable used to synchronize all threads
                ConditionVariable * const global_barrier;

                std::vector<unsigned> thread_ids;

                // Flag whether to use thread affinity
                const bool affinity;

#ifdef linux
                // Array of affinity masks for main process and all controlled threads
                cpu_set_t * affinity_mask;
#endif

            public:
                ThreadPool();
                ~ThreadPool();


                void add_threads(const unsigned num);

                void delete_threads(const unsigned num);

                unsigned get_num_threads() const;

                std::tr1::shared_ptr<Ticket<tags::CPU::MultiCore> > & enqueue(const std::tr1::function<void ()> & task, DispatchPolicy p = DispatchPolicy::any_core());
        };
    }
}
#endif
