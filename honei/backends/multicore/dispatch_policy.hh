/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008, 2009, 2010 Sven Mallach <mallach@honei.org>
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
        class AnyCorePolicy
        {
            private:

            public:

                Ticket<tags::CPU::MultiCore> * operator() (HONEI_UNUSED std::vector<unsigned> & sids)
                {
                    return new Ticket<tags::CPU::MultiCore>();
                }
        };

        class SameCorePolicy
        {
            private:

                Ticket<tags::CPU::MultiCore> * other;

            public:

                SameCorePolicy(Ticket<tags::CPU::MultiCore> * ticket) :
                    other(ticket)
                {
                }

                Ticket<tags::CPU::MultiCore> * operator() (std::vector<unsigned> & sids)
                {
                    unsigned sched_id(other->sid());

                    // Make sure that there is STILL a thread running on that core.
                    if (sids.end() == std::find(sids.begin(), sids.end(), sched_id))
                        sched_id = 0xFFFF;

                    Ticket<tags::CPU::MultiCore> * ticket = new Ticket<tags::CPU::MultiCore>(sched_id, sched_id);

                    return ticket;
                }
        };

        class OnCorePolicy
        {
            private:

                unsigned core_id;

            public:

                OnCorePolicy(const unsigned core_nr) :
                    core_id(core_nr)
                {
                }

                Ticket<tags::CPU::MultiCore> * operator() (std::vector<unsigned> & sids)
                {
                    // Use core_nr as equal to sched_id and make sure that there is a thread on it
                    if (sids.end() == std::find(sids.begin(), sids.end(), core_id))
                        core_id = 0xFFFF;

                    Ticket<tags::CPU::MultiCore> * ticket =
                        new Ticket<tags::CPU::MultiCore>(core_id, core_id);

                    return ticket;
                }
        };

        class SameNodePolicy
        {
            private:

                Ticket<tags::CPU::MultiCore> * other;

            public:

                SameNodePolicy(Ticket<tags::CPU::MultiCore> * ticket) :
                    other(ticket)
                {
                }

                Ticket<tags::CPU::MultiCore> * operator() (HONEI_UNUSED std::vector<unsigned> & sids)
                {
                    unsigned sched_min(other->sid_min());
                    unsigned sched_max(other->sid_max());

                    Ticket<tags::CPU::MultiCore> * ticket = new Ticket<tags::CPU::MultiCore>(sched_min, sched_max);

                    return ticket;
                }
        };

        class OnNodePolicy
        {
            private:

                unsigned node_id;

            public:

                OnNodePolicy(const unsigned node_nr) :
                    node_id(node_nr)
                {
                }

                Ticket<tags::CPU::MultiCore> * operator() (HONEI_UNUSED std::vector<unsigned> & sids)
                {
                    Topology * top = Topology::instance();

                    Ticket<tags::CPU::MultiCore> * ticket =
                        new Ticket<tags::CPU::MultiCore>(top->node_min(node_id), top->node_max(node_id));

                    return ticket;
                }
        };

        /**
         * DispatchPolicy realizes different dispatch strategies
         * for multicore tasks.
         */

        class DispatchPolicy
        {
            private:

                const function<Ticket<tags::CPU::MultiCore> * (std::vector<unsigned> & sids)> policy;

                /// \name Basic Operations
                /// \{

                /// Constructor
                DispatchPolicy(const function<Ticket<tags::CPU::MultiCore> * (std::vector<unsigned> & sids)> p) :
                    policy(p)
                {
                }

                /// \}

            public:

                /// Creates a new Ticket and assures the given way of dispatching
                Ticket<tags::CPU::MultiCore> * apply(std::vector<unsigned> & sids)
                {
                    return policy(sids);
                }

                /// Named constructors

                /// Dispatch on any core available
                static DispatchPolicy any_core()
                {
                    return DispatchPolicy(AnyCorePolicy());
                }

                /// Dispatch on same core as earlier task
                static DispatchPolicy same_core_as(Ticket<tags::CPU::MultiCore> * ticket)
                {
                    return DispatchPolicy(SameCorePolicy(ticket));
                }

                /// Dispatch on explicit core
                static DispatchPolicy on_core(unsigned core_id)
                {
                    return DispatchPolicy(OnCorePolicy(core_id));
                }

                /// Dispatch on same node as earlier task
                static DispatchPolicy same_node_as(Ticket<tags::CPU::MultiCore> * ticket)
                {
                    return DispatchPolicy(SameNodePolicy(ticket));
                }

                /// Dispatch on explicit NUMA node
                static DispatchPolicy on_node(unsigned node_id)
                {
                    return DispatchPolicy(OnNodePolicy(node_id));
                }
        };
    }
}
#endif
