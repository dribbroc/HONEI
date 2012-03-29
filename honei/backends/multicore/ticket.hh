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
#ifndef HONEI_GUARD_MULTICORE_TICKET_HH
#define HONEI_GUARD_MULTICORE_TICKET_HH 1

#include <honei/util/ticket.hh>

namespace honei
{
        template <> class Ticket<tags::CPU::MultiCore> :
            public TicketBase,
            public PrivateImplementationPattern<Ticket<tags::CPU::MultiCore>, Shared>
        {
            private:
                /// \name Private members
                /// \{

                /// \}

            public:
                /// \name Basic Operations
                /// \{

                /// Constructor
                Ticket(const unsigned req_socket = 0xFFFF, const unsigned req_sched = 0xFFFF);

                /// Destructor
                virtual ~Ticket();

                /// \}

                /// Derive mark and wait from TicketBase

                Ticket<tags::CPU::MultiCore> & operator= (const Ticket<tags::CPU::MultiCore> & other);

                /// Retrieve unique ticket ID
                unsigned uid() const;

                /// Retrieve the socket the associated task shall be assigned to
                unsigned & req_socket() const;

                /// Retrieve the LPU the associated task shall be assigned to
                unsigned & req_sched() const;

                /// Retrieve sched ID of the thread executing the task
                unsigned & sid();

                /// Retrieve thread ID
                unsigned & tid();
       };
}
#endif
