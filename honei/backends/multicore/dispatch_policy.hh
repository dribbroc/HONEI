/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Sven Mallach <sven.mallach@cs.tu-dortmund.de>
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

#ifndef MULTICORE_GUARD_POLICY_HH
#define MULTICORE_GUARD_POLICY_HH 1

#include <honei/backends/multicore/multicore_ticket.hh>

#include <map>
#include <tr1/functional>

namespace honei
{
    namespace mc
    {
        class AnyCorePolicy
        {
            private:

            public:

                MultiCoreTicket * operator() (unsigned * tid_array)
                {
                    return new MultiCoreTicket;
                }
        };

        class SameCorePolicy
        {
            private:

                MultiCoreTicket & other;

            public:

                SameCorePolicy(MultiCoreTicket & ticket) :
                    other(ticket)
                {
                }

                MultiCoreTicket * operator() (unsigned * tid_array)
                {
                    MultiCoreTicket * ticket = new MultiCoreTicket(other.tid());
                    return ticket;
                }
        };

        class OnCorePolicy
        {
            private:

                const unsigned core_id;

            public:

                OnCorePolicy(const unsigned core_nr) :
                    core_id(core_nr)
                {
                }

                MultiCoreTicket * operator() (unsigned * tid_array)
                {
                    MultiCoreTicket * ticket = new MultiCoreTicket(tid_array[core_id]);
                    return ticket;
                }
        };

        class DispatchPolicy
        {
            private:

                const std::tr1::function<MultiCoreTicket * (unsigned *)> policy;

                DispatchPolicy(const std::tr1::function<MultiCoreTicket * (unsigned *)> p) :
                    policy(p)
                {
                }

            public:

                MultiCoreTicket * apply(unsigned * tid_array)
                {
                    return policy(tid_array);
                }

                /// Named constructors

                /// Dispatch on any core available
                static DispatchPolicy any_core()
                {
                    return DispatchPolicy(AnyCorePolicy());
                }

                /// Dispatch on same core as earlier task
                static DispatchPolicy same_core_as(MultiCoreTicket & ticket)
                {
                    return DispatchPolicy(SameCorePolicy(ticket));
                }

                /// Dispatch on explicit core
                static DispatchPolicy on_core(unsigned core_id)
                {
                    return DispatchPolicy(OnCorePolicy(core_id));
                }
        };
    }
}
#endif
