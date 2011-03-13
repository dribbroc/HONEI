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

#include <honei/backends/multicore/atomic_slist.hh>
#include <honei/backends/multicore/thread_task.hh>
#include <honei/util/condition_variable.hh>
#include <honei/util/private_implementation_pattern.hh>
#include <honei/util/thread.hh>

#include <list>
#include <vector>

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

                AffinityThreadFunction(PoolSyncData * const psync, std::list<ThreadTask *> * const list, unsigned pool_id, unsigned sched_id);

                virtual ~AffinityThreadFunction();

                /// The threads' main function
                virtual void operator() ();

                virtual void stop();

                virtual unsigned tid() const;

                unsigned pool_id() const;
        };

        class SimpleThreadFunction :
            public ThreadFunctionBase,
            public PrivateImplementationPattern<SimpleThreadFunction, Shared>
        {
            private:

            public:

                SimpleThreadFunction(PoolSyncData * const psync, AtomicSList<ThreadTask *> * const list, unsigned pool_id);

                virtual ~SimpleThreadFunction();

                /// The threads' main function
                virtual void operator() ();

                virtual void stop();

                virtual unsigned tid() const;
        };

        class WorkStealingThreadFunction :
            public ThreadFunctionBase,
            public PrivateImplementationPattern<WorkStealingThreadFunction, Shared>
        {
            private:

            public:

                WorkStealingThreadFunction(PoolSyncData * const psync, unsigned pool_id, unsigned sched_id, const std::vector<std::pair<Thread *, mc::WorkStealingThreadFunction *> > & threads, unsigned num_thr, volatile bool & terminate);

                virtual ~WorkStealingThreadFunction();

                /// The threads' main function
                virtual void operator() ();

                virtual void stop();

                void enqueue(mc::ThreadTask * task);

                virtual unsigned tid() const;

                unsigned pool_id() const;
        };
    }
}
#endif
