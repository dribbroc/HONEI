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

#include <honei/util/ticket.hh>

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

                unsigned & operator() (std::map<unsigned, unsigned> & mapping, Ticket * ticket, unsigned * tid_array)
                {
                    return mapping[ticket->id()] = 0;
                }
        };

        class SameCorePolicy
        {
            private:

                const unsigned ticket_id;

            public:

                SameCorePolicy(const unsigned ticket_nr) :
                    ticket_id(ticket_nr)
                {
                }

                unsigned & operator() (std::map<unsigned, unsigned> & mapping, Ticket * ticket, unsigned * tid_array)
                {
                    unsigned & thread_id = mapping[ticket_id];
                    return mapping[ticket->id()] = thread_id;
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

                unsigned & operator() (std::map<unsigned, unsigned> & mapping, Ticket * ticket, unsigned * tid_array)
                {
                    return mapping[ticket->id()] = tid_array[core_id];
                }
        };

        class DispatchPolicy
        {
            private:

                // Mapping of ticket IDs to thread IDs
              static std::map<unsigned, unsigned> mapping;

                const std::tr1::function<unsigned & (std::map<unsigned, unsigned> &, Ticket *, unsigned *)> policy;

                DispatchPolicy(const std::tr1::function<unsigned & (std::map<unsigned, unsigned> &, Ticket *, unsigned *)> p) :
                    policy(p)
                {
                }

            public:

                unsigned & apply(Ticket * ticket, unsigned * tid_array)
                {
                    return policy(DispatchPolicy::mapping, ticket, tid_array);
                }

                /// Named constructors

                /// Dispatch on any core available
                static DispatchPolicy any_core()
                {
                    return DispatchPolicy(AnyCorePolicy());
                }

                /// Dispatch on same core as earlier task
                static DispatchPolicy same_core_as(unsigned ticket_id)
                {
                    return DispatchPolicy(SameCorePolicy(ticket_id));
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
