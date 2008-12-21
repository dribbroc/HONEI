/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#include <honei/util/assertion.hh>
#include <honei/util/condition_variable.hh>
#include <honei/util/exception.hh>
#include <honei/util/lock.hh>
#include <honei/util/mutex.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <honei/util/ticket.hh>

#include <list>
#include <tr1/memory>

namespace honei
{
    Ticket::Ticket() :
        PrivateImplementationPattern<Ticket, Shared>(new Implementation<Ticket>(counter++))
    {
    }

    Ticket::~Ticket()
    {
    }

    void
    Ticket::mark()
    {
        CONTEXT("When marking a ticket for completion:");
        Lock l(_imp->mutex);

        ASSERT(! _imp->completed, "ticket marked more than once!");

        _imp->completed = true;
        _imp->completion.signal();
    }

    void
    Ticket::wait() const
    {
        CONTEXT("When waiting for ticket completion:");
        Lock l(_imp->mutex);

        while (! _imp->completed)
        {
            _imp->completion.wait(_imp->mutex);
        }
    }

    unsigned
    Ticket::id() const
    {
        return _imp->uid;
    }

    unsigned Ticket::counter(0);

    template <> struct Implementation<TicketList>
    {
        std::list<std::tr1::shared_ptr<Implementation<Ticket> > > tickets;
    };

    TicketList::TicketList() :
        PrivateImplementationPattern<TicketList, Single>(new Implementation<TicketList>)
    {
    }

    TicketList::~TicketList()
    {
    }

    void
    TicketList::push_back(const Ticket & ticket)
    {
        CONTEXT("When pushing a Ticket to the back of a TicketList:");

        _imp->tickets.push_back(ticket._imp);
    }

    bool
    TicketList::empty()
    {
        return _imp->tickets.empty();
    }

    const unsigned
    TicketList::remove_front()
    {
        std::tr1::shared_ptr<Implementation<Ticket> > ticket(_imp->tickets.front());
        _imp->tickets.pop_front();
        return ticket->uid;
    }

    void
    TicketList::wait() const
    {
        for (std::list<std::tr1::shared_ptr<Implementation<Ticket> > >::iterator i(_imp->tickets.begin()) , i_end(_imp->tickets.end()) ; i != i_end ; ++i)
        {
            std::tr1::shared_ptr<Implementation<Ticket> > ticket(*i);
            Lock l(ticket->mutex);

            while (! ticket->completed)
            {
                ticket->completion.wait(ticket->mutex);
            }
        }
    }
}
