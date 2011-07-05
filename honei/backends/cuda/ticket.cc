/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@tu-dortmund.de>
 * Copyright (c) 2011 Sven Mallach <mallach@honei.org>
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

#include <honei/backends/cuda/ticket.hh>
#include <honei/util/assertion.hh>
#include <honei/util/condition_variable.hh>
#include <honei/util/exception.hh>
#include <honei/util/lock.hh>
#include <honei/util/private_implementation_pattern-impl.hh>

namespace honei
{
        template <> struct Implementation<Ticket<tags::GPU::MultiCore> > :
            public TicketBaseImpl
        {
            /// Counter for unique global ticket IDs
            static unsigned counter;

            /// Unique ID
            unsigned id;

            /// ID of the executing thread
            unsigned thread_id;

            Implementation() :
                TicketBaseImpl(),
                id(counter)
            {
                ++counter;
            }

            virtual ~Implementation()
            {
            }
        };

        Ticket<tags::GPU::MultiCore>::Ticket() :
            TicketBase(),
            PrivateImplementationPattern<Ticket<tags::GPU::MultiCore>, Shared>(new Implementation<Ticket<tags::GPU::MultiCore> >())
        {
        }

        Ticket<tags::GPU::MultiCore>::~Ticket()
        {
        }

        unsigned Ticket<tags::GPU::MultiCore>::uid() const
        {
            return _imp->id;
        }

        unsigned Implementation<Ticket<tags::GPU::MultiCore> >::counter(0);
}
