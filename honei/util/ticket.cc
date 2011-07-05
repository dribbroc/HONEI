/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2008, 2009, 2011 Sven Mallach <mallach@honei.org>
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

#include <honei/util/condition_variable.hh>
#include <honei/util/lock.hh>
#include <honei/util/mutex.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <honei/util/ticket.hh>

#include <vector>

namespace honei
{
    /* TicketBase member implementation */

    TicketBase::TicketBase() :
        _base(new TicketBaseImpl)
    {
    }

    TicketBase::~TicketBase()
    {
    }

    void TicketBase::mark() const
    {
        _base->mark();
    }

    void TicketBase::wait() const
    {
        _base->wait();
    }



    template <> struct Implementation<Ticket<tags::CPU> > :
        public TicketBaseImpl
    {
        Implementation() :
            TicketBaseImpl()
        {
        }

        virtual ~Implementation()
        {
        }
    };


    /* Ticket<tags::CPU>'s implementation - just adopting
     * the basic implementation from TicketBase */

    Ticket<tags::CPU>::Ticket() :
        TicketBase(),
        PrivateImplementationPattern<Ticket<tags::CPU>, Shared>(new Implementation<Ticket<tags::CPU> >)
    {
    }

    Ticket<tags::CPU>::~Ticket()
    {
    }

    /* Implementation of TicketVector */

    template <> struct Implementation<TicketVector>
    {
        std::vector<TicketBase> tickets;
    };

    TicketVector::TicketVector() :
        PrivateImplementationPattern<TicketVector, Single>(new Implementation<TicketVector>)
    {
    }

    TicketVector::~TicketVector()
    {
    }

    void
    TicketVector::push_back(TicketBase ticket)
    {
        CONTEXT("When pushing a Ticket to the back of a TicketVector:");

        _imp->tickets.push_back(ticket);
    }

    void
    TicketVector::wait() const
    {
        CONTEXT("When waiting for the Tickets in a TicketVector:");

        while (! _imp->tickets.empty())
        {
            _imp->tickets.back().wait();
            _imp->tickets.pop_back();
        }
    }
}
