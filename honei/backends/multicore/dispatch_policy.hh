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

#ifndef MULTICORE_GUARD_POLICY_HH
#define MULTICORE_GUARD_POLICY_HH 1

#include <honei/backends/multicore/ticket.hh>
#include <honei/util/attributes.hh>

#include <tr1/functional>
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

                Ticket<tags::CPU::MultiCore> * operator() ( HONEI_UNUSED std::vector<unsigned> & tids)
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

                Ticket<tags::CPU::MultiCore> * operator() (std::vector<unsigned> & tids)
                {
                    unsigned thread_id(other->tid());
                    if (tids.end() == std::find(tids.begin(), tids.end(), thread_id)) // Thread could already have been deleted
                        thread_id = 0;

                    Ticket<tags::CPU::MultiCore> * ticket = new Ticket<tags::CPU::MultiCore>(thread_id);

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

                Ticket<tags::CPU::MultiCore> * operator() (std::vector<unsigned> & tids)
                {
                    Ticket<tags::CPU::MultiCore> * ticket =
                        new Ticket<tags::CPU::MultiCore>(tids[core_id % tids.size()]);

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

                const std::tr1::function<Ticket<tags::CPU::MultiCore> * (std::vector<unsigned> & tids)> policy;

                /// \name Basic Operations
                /// \{

                /// Constructor
                DispatchPolicy(const std::tr1::function<Ticket<tags::CPU::MultiCore> * (std::vector<unsigned> & tids)> p) :
                    policy(p)
                {
                }

                /// \}

            public:

                /// Creates a new Ticket and assures the given way of dispatching
                Ticket<tags::CPU::MultiCore> * apply(std::vector<unsigned> & tids)
                {
                    return policy(tids);
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
        };
    }
}
#endif
