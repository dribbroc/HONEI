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
#include <honei/util/private_implementation_pattern.hh>
#include <honei/util/thread.hh>

#include <vector>

namespace honei
{
    namespace mc
    {
        class ThreadPool :
            public PrivateImplementationPattern<ThreadPool, Single>,
            public InstantiationPolicy<ThreadPool, Singleton>
        {
            private:

                /// \name Private members
                /// \{
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
