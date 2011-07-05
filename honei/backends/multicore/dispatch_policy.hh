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
#ifndef MULTICORE_GUARD_POLICY_HH
#define MULTICORE_GUARD_POLICY_HH 1

#include <honei/backends/multicore/ticket.hh>
#include <honei/backends/multicore/topology.hh>
#include <honei/util/attributes.hh>
#include <honei/util/tr1_boost.hh>
#include <vector>
#include <algorithm>

namespace honei
{
    namespace mc
    {
        /* Dispatch on an arbitrary processing unit.
         * It will usually be the one that is able to
         * response at fastest */

        class AnyCorePolicy
        {
            private:

            public:

                Ticket<tags::CPU::MultiCore> operator() ()
                {
                    Ticket<tags::CPU::MultiCore> ticket;
                    return ticket;
                }
        };

        /* Dispatch on the same processing unit as
         * with a previous ticket gained by OnCorePolicy,
         * which is to provide here */

        class SameCorePolicy
        {
            private:

                Ticket<tags::CPU::MultiCore> other;

            public:

                SameCorePolicy(Ticket<tags::CPU::MultiCore> & ticket) :
                    other(ticket)
                {
                }

                Ticket<tags::CPU::MultiCore> operator() ()
                {
                    unsigned sched_id(other.sid());

                    Ticket<tags::CPU::MultiCore> ticket(sched_id, sched_id);

                    return ticket;
                }
        };

        /* Dispatch on an explicit processing unit */

        class OnCorePolicy
        {
            private:

                unsigned core_id;

            public:

                OnCorePolicy(const unsigned core_nr) :
                    core_id(core_nr)
                {
                }

                Ticket<tags::CPU::MultiCore> operator() ()
                {
                    Topology * top = Topology::instance();

                    if (core_id > top->num_lpus() - 1)
                        core_id = top->num_lpus() - 1;

                    Ticket<tags::CPU::MultiCore> ticket(core_id, core_id);

                    return ticket;
                }
        };

        /* Dispatch on an arbitrary core on the same node as
         * with a previous ticket gained by OnNodePolicy,
         * which is to provide here */

        class SameNodePolicy
        {
            private:

                Ticket<tags::CPU::MultiCore> other;

            public:

                SameNodePolicy(Ticket<tags::CPU::MultiCore> & ticket) :
                    other(ticket)
                {
                }

                Ticket<tags::CPU::MultiCore> operator() ()
                {
                    unsigned sched_min(other.sid_min());
                    unsigned sched_max(other.sid_max());

                    Ticket<tags::CPU::MultiCore> ticket(sched_min, sched_max);

                    return ticket;
                }
        };

        /* Dispatch on an arbitrary core on an explicit node. */

        class OnNodePolicy
        {
            private:

                unsigned node_id;

            public:

                OnNodePolicy(const unsigned node_nr) :
                    node_id(node_nr)
                {
                }

                Ticket<tags::CPU::MultiCore> operator() ()
                {
                    Topology * top = Topology::instance();

                    Ticket<tags::CPU::MultiCore> ticket(top->node_min(node_id), top->node_max(node_id));

                    return ticket;
                }
        };

        /* Fill the avaiable nodes processing unit by processing unit */

        class LinearNodePolicy
        {
            private:

                Ticket<tags::CPU::MultiCore> last;

            public:

                LinearNodePolicy(Ticket<tags::CPU::MultiCore> & l) :
                    last(l)
                {
                }

                Ticket<tags::CPU::MultiCore> operator() ()
                {
                    Topology * top = Topology::instance();

                    if (last.uid() == 0)
                    {
                        Ticket<tags::CPU::MultiCore> ticket(0, 0);
                        return ticket;
                    }
                    else
                    {
                        unsigned last_core = last.sid_min();

                        if (last_core == top->num_lpus() - 1)
                            last_core = 0;
                        else
                            ++last_core;

                        Ticket<tags::CPU::MultiCore> ticket(last_core, last_core);
                        return ticket;
                    }
                }
        };

        /* Fill the avaiable nodes in an alternaing manner concerning
         * their processing units */

        class AlternatingNodePolicy
        {
            private:

                Ticket<tags::CPU::MultiCore> last;

            public:

                AlternatingNodePolicy(Ticket<tags::CPU::MultiCore> & l) :
                    last(l)
                {
                }

                Ticket<tags::CPU::MultiCore> operator() ()
                {
                    Topology * top = Topology::instance();
                    unsigned num_nodes = top->num_nodes();

                    if (last.uid() == 0)
                    {
                        Ticket<tags::CPU::MultiCore> ticket(0, 0);
                        return ticket;
                    }
                    else
                    {
                        unsigned last_node = top->get_node(last.sid_min());
                        unsigned core_pos = last.sid_min() % top->lpus_per_node();

                        unsigned next_node(1 + last_node);

                        if (last_node == num_nodes - 1)
                        {
                            next_node = 0; // overwrite increase and reset to 0

                            if (core_pos == top->lpus_per_node() - 1)
                                core_pos = 0;
                            else
                                ++core_pos;
                        }
                        // else leave increase as set above

                        unsigned next_core = next_node * top->lpus_per_node() + core_pos;

                        Ticket<tags::CPU::MultiCore> ticket(next_core, next_core);
                        return ticket;
                    }
                }
        };

        /**
         * DispatchPolicy realizes different dispatch strategies
         * for multicore tasks.
         */

        class DispatchPolicy
        {
            private:

                const function<Ticket<tags::CPU::MultiCore> ()> policy;

                /// \name Basic Operations
                /// \{

                /// Constructor
                DispatchPolicy(const function<Ticket<tags::CPU::MultiCore> ()> & p) :
                    policy(p)
                {
                }

                /// \}

            public:

                static Ticket<tags::CPU::MultiCore> last; // Store history

                /// Creates a new Ticket and assures the given way of dispatching
                Ticket<tags::CPU::MultiCore> apply()
                {
                    last = policy();
                    return last;
                }

                /// \name Named constructors
                /// \{

                /// \name Policies that are able to be used as default
                /// \{

                /// Dispatch on any core available
                static DispatchPolicy any_core()
                {
                    return DispatchPolicy(AnyCorePolicy());
                }

                /// Dispatch alternatingly on next NUMA node
                static DispatchPolicy alternating_node()
                {
                    return DispatchPolicy(AlternatingNodePolicy(last));
                }

                /// Dispatch linearly on NUMA nodes
                static DispatchPolicy linear_node()
                {
                    return DispatchPolicy(LinearNodePolicy(last));
                }

                /// \}

                /// \name Policies that are designed to be used
                /// explicitly as a parameter to ThreadPool::enqueue()
                /// \{

                /// Dispatch on same core as earlier task
                static DispatchPolicy same_core_as(Ticket<tags::CPU::MultiCore> & ticket)
                {
                    return DispatchPolicy(SameCorePolicy(ticket));
                }

                /// Dispatch on explicit core
                static DispatchPolicy on_core(unsigned core_id)
                {
                    return DispatchPolicy(OnCorePolicy(core_id));
                }

                /// Dispatch on same node as earlier task
                static DispatchPolicy same_node_as(Ticket<tags::CPU::MultiCore> & ticket)
                {
                    return DispatchPolicy(SameNodePolicy(ticket));
                }

                /// Dispatch on explicit NUMA node
                static DispatchPolicy on_node(unsigned node_id)
                {
                    return DispatchPolicy(OnNodePolicy(node_id));
                }

                /// \}
        };
    }
}
#endif
