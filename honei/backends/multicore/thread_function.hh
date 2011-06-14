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

#pragma once
#ifndef MULTICORE_GUARD_THREAD_FUNCTION_HH
#define MULTICORE_GUARD_THREAD_FUNCTION_HH 1

#include <honei/backends/multicore/concurrent_list.hh>
#include <honei/backends/multicore/thread_task.hh>
#include <honei/util/condition_variable.hh>
#include <honei/util/private_implementation_pattern.hh>
#include <honei/util/thread.hh>

#include <list>
#include <vector>
#include <deque>

namespace honei
{
    namespace mc
    {
        class ThreadFunctionBase
        {
            private:

            public:

                virtual ~ThreadFunctionBase();

                virtual void operator() () = 0;

                virtual void stop() = 0;

                virtual unsigned tid() const = 0;
        };

        class AffinityThreadFunction :
            public ThreadFunctionBase,
            public PrivateImplementationPattern<AffinityThreadFunction, Shared>
        {
            private:

            public:

                AffinityThreadFunction(PoolSyncData * const psync, std::deque<ThreadTask *> * const list, unsigned pool_id, unsigned sched_id);

                virtual ~AffinityThreadFunction();

                /// The threads' main function
                virtual void operator() ();

                virtual void stop();

                virtual unsigned tid() const;

                unsigned pool_id() const;
        };

        template <typename ListType> class SimpleThreadFunction;

        template <> class SimpleThreadFunction<ConcurrentList<ThreadTask *> > :
            public ThreadFunctionBase,
            public PrivateImplementationPattern<SimpleThreadFunction<ConcurrentList<ThreadTask *> >, Shared>
        {
            private:

            public:

                SimpleThreadFunction(PoolSyncData * const psync, ConcurrentList<ThreadTask *> * const list, unsigned pool_id);

                virtual ~SimpleThreadFunction();

                /// The threads' main function
                virtual void operator() ();

                virtual void stop();

                virtual unsigned tid() const;
        };

        template <> class SimpleThreadFunction<std::deque<ThreadTask *> > :
            public ThreadFunctionBase,
            public PrivateImplementationPattern<SimpleThreadFunction<std::deque<ThreadTask *> >, Shared>
        {
            private:

            public:

                SimpleThreadFunction(PoolSyncData * const psync, std::deque<ThreadTask *> * const list, unsigned pool_id);

                virtual ~SimpleThreadFunction();

                /// The threads' main function
                virtual void operator() ();

                virtual void stop();

                virtual unsigned tid() const;
        };

        template <typename ListType> class WorkStealingThreadFunction;

        template <> class WorkStealingThreadFunction<std::deque<mc::ThreadTask *> > :
            public ThreadFunctionBase,
            public PrivateImplementationPattern<WorkStealingThreadFunction<std::deque<mc::ThreadTask *> >, Shared>
        {
            private:

            public:

                WorkStealingThreadFunction(PoolSyncData * const psync, unsigned pool_id, unsigned sched_id,
                        const std::vector<std::pair<Thread *, mc::WorkStealingThreadFunction<std::deque<mc::ThreadTask *> > *> > & threads,
                        unsigned num_thr, volatile bool & terminate, Mutex * const steal_mutex);

                virtual ~WorkStealingThreadFunction();

                /// The threads' main function
                virtual void operator() ();

                virtual void stop();

                void enqueue(mc::ThreadTask * task);

                bool steal(std::deque<mc::ThreadTask *> & thief_list);

                virtual unsigned tid() const;

                unsigned pool_id() const;
        };

        template <> class WorkStealingThreadFunction<ConcurrentList<mc::ThreadTask *> > :
            public ThreadFunctionBase,
            public PrivateImplementationPattern<WorkStealingThreadFunction<ConcurrentList<mc::ThreadTask *> >, Shared>
        {
            private:

            public:

                WorkStealingThreadFunction(PoolSyncData * const psync, unsigned pool_id, unsigned sched_id,
                        const std::vector<std::pair<Thread *, mc::WorkStealingThreadFunction<ConcurrentList<mc::ThreadTask *> > *> > & threads,
                        unsigned num_thr, volatile bool & terminate, Mutex * const steal_mutex);

                virtual ~WorkStealingThreadFunction();

                /// The threads' main function
                virtual void operator() ();

                virtual void stop();

                void enqueue(mc::ThreadTask * task);

                bool steal(mc::ConcurrentList<mc::ThreadTask *> & thief_list);

                virtual unsigned tid() const;

                unsigned pool_id() const;
        };
    }
}
#endif
