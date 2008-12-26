/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#ifndef HONEI_GUARD_UTIL_TICKET_HH
#define HONEI_GUARD_UTIL_TICKET_HH 1

#include <honei/util/condition_variable.hh>
#include <honei/util/private_implementation_pattern.hh>

namespace honei
{
    // Forward declaration for Ticket.
    class TicketList;

    /**
     * Ticket is used by all asynchronous function calls to relay/query
     * information on a function's completion status.
     */
    class Ticket :
        public PrivateImplementationPattern<Ticket, Shared>
    {
        private:
            /// \name Private members
            /// \{
            /// \}

        public:
            /// \name Friends of Ticket
            /// \{

            friend class TicketList;

            /// \}

            /// \name Basic Operations
            /// \{

            /// Constructor.
            Ticket();

            /// Destructor.
            ~Ticket();

            /// \}

            /// Mark ticket as completed.
            void mark();

            /// Wait for ticket completion.
            void wait() const;
   };

    /**
     * TicketList is used when multiple tickets are needed.
     */
    class TicketList :
        public InstantiationPolicy<TicketList, NonCopyable>,
        public PrivateImplementationPattern<TicketList, Single>
    {
        public:
            /// \name Basic Operations
            /// \{

            /// Constructor.
            TicketList();

            /// Destructor.
            ~TicketList();

            /// \}

            /// Push a ticket to the back of the list.
            void push_back(const Ticket & ticket);

            /// Wait for ticket completion.
            void wait() const;

            bool empty();

            const unsigned remove_front();
    };
}

#endif
