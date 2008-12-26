/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Sven Mallach <sven.mallach@cs.uni-dortmund.de>
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

#ifndef HONEI_GUARD_MULTICORE_TICKET_HH
#define HONEI_GUARD_MULTICORE_TICKET_HH 1

#include <honei/util/ticket.hh>

namespace honei
{
    namespace mc
    {
        class MultiCoreTicket :
            public Ticket
        {
            private:
            /// \name Private members
            /// \{

            /// Counter for unique global ticket IDs
            static unsigned counter;

            /// Unique ID
            unsigned id;

            /// ID of the executing thread
            unsigned thread_id;

            /// \}

        public:
            /// \name Basic Operations
            /// \{

            /// Constructor.
            MultiCoreTicket(const unsigned tid = 0);

            /// Destructor.
            ~MultiCoreTicket();

            /// \}

            /// Retrieve unique ticket ID
            unsigned uid() const;

            /// Retrieve thread ID
            unsigned & tid();
       };
    }
}

#endif
