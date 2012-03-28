/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 - 2012 Sven Mallach <mallach@honei.org>
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
#ifndef MULTICORE_GUARD_DISPATCH_POLICY_HH
#define MULTICORE_GUARD_DISPATCH_POLICY_HH 1

#include <honei/backends/multicore/ticket.hh>
#include <honei/backends/multicore/topology.hh>
#include <honei/util/attributes.hh>
#include <honei/util/log.hh>
#include <honei/util/stringify.hh>
#include <honei/util/tr1_boost.hh>
#include <vector>
#include <algorithm>
#include <iostream>

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
                    CONTEXT("When using the AnyCorePolicy to dispatch a task:");

                    Ticket<tags::CPU::MultiCore> ticket;

#ifdef DEBUG
                    std::string msg = "Dispatching ticket " + stringify(ticket.uid()) + " with the AnyCore Policy \n";
                    LOGMESSAGE(lc_backend, msg);
#endif

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
                    CONTEXT("When using the OnCorePolicy to dispatch a task:");

                    Topology * top = Topology::instance();

                    if (core_id >= top->num_lpus() || (! top->lpu(core_id)->has_thread))
                        core_id = 0;

                    Ticket<tags::CPU::MultiCore> ticket(0xFFFF, core_id);

#ifdef DEBUG
                    std::string msg = "Dispatching ticket " + stringify(ticket.uid()) + " on LPU " + stringify(core_id) + "\n";
                    LOGMESSAGE(lc_backend, msg);
#endif

                    return ticket;
                }
        };

        /* Dispatch on the same processing unit as with a
         * previous ticket which is to provide here */

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
                    CONTEXT("When using the SameCorePolicy to dispatch a task:");

                    Ticket<tags::CPU::MultiCore> ticket(0xFFFF, other.req_sched());

#ifdef DEBUG
                    std::string msg = "Dispatching ticket " + stringify(ticket.uid()) + " on LPU " + stringify(other.req_sched()) + "\n";
                    LOGMESSAGE(lc_backend, msg);
#endif

                    return ticket;
                }
        };

        /* Dispatch on an arbitrary core on an explicit socket. */

        class OnSocketPolicy
        {
            private:

                unsigned socket_id;

            public:

                OnSocketPolicy(const unsigned socket_nr) :
                    socket_id(socket_nr)
                {
                }

                Ticket<tags::CPU::MultiCore> operator() ()
                {
                    CONTEXT("When using the OnSocketPolicy to dispatch a task:");

                    Topology * top = Topology::instance();

                    if (socket_id >= top->num_cpus() || (top->sockets()[socket_id]->_num_threads == 0))
                        socket_id = 0;

                    Ticket<tags::CPU::MultiCore> ticket(socket_id, 0xFFFF);

#ifdef DEBUG
                    std::string msg = "Dispatching ticket " + stringify(ticket.uid()) + " on socket " + stringify(socket_id) + "\n";
                    LOGMESSAGE(lc_backend, msg);
#endif

                    return ticket;
                }
        };

        /* Dispatch on an arbitrary core on the same node as
         * with a previous ticket which is to provide here */

        class SameSocketPolicy
        {
            private:

                Ticket<tags::CPU::MultiCore> other;

            public:

                SameSocketPolicy(Ticket<tags::CPU::MultiCore> & ticket) :
                    other(ticket)
                {
                }

                Ticket<tags::CPU::MultiCore> operator() ()
                {
                    CONTEXT("When using the SameSocketPolicy to dispatch a task:");

                    unsigned socket;

                    if (other.req_sched() == 0xFFFF)
                    {
                        socket = 0xFFFF;
                    }
                    else
                    {
                       LPU * lpu = Topology::instance()->lpu(other.req_sched());
                       socket = lpu->socket_id;
                    }

                    Ticket<tags::CPU::MultiCore> ticket(socket, 0xFFFF);

#ifdef DEBUG
                    std::string msg = "Dispatching ticket " + stringify(ticket.uid()) + " on socket " + stringify(socket) + "\n";
                    LOGMESSAGE(lc_backend, msg);
#endif

                    return ticket;
                }
        };


        /* Fill the available sockets processing unit by processing unit */

        class LinearSocketPolicy
        {
            private:

                Ticket<tags::CPU::MultiCore> last;

            public:

                LinearSocketPolicy(Ticket<tags::CPU::MultiCore> & l) :
                    last(l)
                {
                }

                Ticket<tags::CPU::MultiCore> operator() ()
                {
                    CONTEXT("When using the LinearSocketPolicy to dispatch a task:");

                    // Only happens at the very beginning...
                    if (last.req_sched() >= Topology::instance()->num_lpus())
                    {
                        Ticket<tags::CPU::MultiCore> ticket(0, 0);
#ifdef DEBUG
                        std::string msg = "Dispatching ticket " + stringify(ticket.uid()) + " on LPU 0" + "\n";
                        LOGMESSAGE(lc_backend, msg);
#endif
                        return ticket;
                    }

                    LPU * lpu = Topology::instance()->lpu(last.req_sched());
                    Ticket<tags::CPU::MultiCore> ticket(lpu->socket_id, lpu->linear_enqueue_succ->sched_id);

#ifdef DEBUG
                    std::string msg = "Dispatching ticket " + stringify(ticket.uid()) + " on LPU " + stringify(lpu->linear_enqueue_succ->sched_id) + "\n";
                    LOGMESSAGE(lc_backend, msg);
#endif
                    return ticket;
                }
        };

        /* Fill the available sockets in an alternating manner concerning
         * their processing units */

        class AlternatingSocketPolicy
        {
            private:

                Ticket<tags::CPU::MultiCore> last;

            public:

                AlternatingSocketPolicy(Ticket<tags::CPU::MultiCore> & l) :
                    last(l)
                {
                }

                Ticket<tags::CPU::MultiCore> operator() ()
                {
                    CONTEXT("When using the AlternatingSocketPolicy to dispatch a task:");

                    // Only happens at the very beginning....
                    if (last.req_sched() >= Topology::instance()->num_lpus())
                    {
                        Ticket<tags::CPU::MultiCore> ticket(0, 0);
#ifdef DEBUG
                        std::string msg = "Dispatching ticket " + stringify(ticket.uid()) + " on LPU 0" + "\n";
                        LOGMESSAGE(lc_backend, msg);
#endif
                        return ticket;
                    }

                    LPU * lpu = Topology::instance()->lpu(last.req_sched());
                    Ticket<tags::CPU::MultiCore> ticket(lpu->alternating_enqueue_succ->socket_id, lpu->alternating_enqueue_succ->sched_id);

#ifdef DEBUG
                    std::string msg = "Dispatching ticket " + stringify(ticket.uid()) + " on LPU " + stringify(lpu->alternating_enqueue_succ->sched_id) + " (last ticket was " + stringify(last.uid()) + " with requested sched " + stringify(last.req_sched()) + ")\n";
                    LOGMESSAGE(lc_backend, msg);
#endif
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

                /// Dispatch alternatingly on next (NUMA) socket
                static DispatchPolicy alternating_socket()
                {
                    return DispatchPolicy(AlternatingSocketPolicy(last));
                }

                /// Dispatch linearly on (NUMA) socket
                static DispatchPolicy linear_socket()
                {
                    return DispatchPolicy(LinearSocketPolicy(last));
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
                static DispatchPolicy same_socket_as(Ticket<tags::CPU::MultiCore> & ticket)
                {
                    return DispatchPolicy(SameSocketPolicy(ticket));
                }

                /// Dispatch on explicit NUMA node
                static DispatchPolicy on_socket(unsigned socket_id)
                {
                    return DispatchPolicy(OnSocketPolicy(socket_id));
                }

                /// \}
        };
    }
}
#endif
