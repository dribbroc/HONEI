/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#ifndef HONEI_GUARD_UTIL_TICKET_HH
#define HONEI_GUARD_UTIL_TICKET_HH 1

#include <honei/util/condition_variable.hh>
#include <honei/util/private_implementation_pattern.hh>
#include <honei/util/tags.hh>

namespace honei
{
    // Forward declaration for Ticket.
    class TicketList;

    /**
     * TicketBase is the virtual base class for all types
     * of Tickets implemented for HONEIs backends.
     */

    class TicketBase
    {
        public:

            /// \name Basic Operations
            /// \{

            /// Destructor.
            virtual ~TicketBase()
            {
            }
            /// \}

            /// Mark ticket as completed.
            virtual void mark() = 0;

            /// Wait for ticket completion.
            virtual void wait() const = 0;

            /// \name Friends of TicketBase
            /// \{
            friend class TicketList;
            /// \}
    };

    /**
     * Ticket is used by all asynchronous function calls to relay/query
     * information on a function's completion status.
     */

    template <typename Tag_> class Ticket;

    template <> class Ticket<tags::CPU> :
        public TicketBase,
        public PrivateImplementationPattern<Ticket<tags::CPU>, Shared>
    {
        /// \name Basic Operations
        /// \{

        /// Constructor.
        Ticket();

        /// \}

        /// Mark ticket as completed.
        virtual void mark();

        /// Wait for ticket completion.
        virtual void wait() const;
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
            void push_back(const std::tr1::shared_ptr<TicketBase> & ticket);

            /// Wait for ticket completion.
            void wait() const;
    };
}

#endif
