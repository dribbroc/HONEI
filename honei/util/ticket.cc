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
    template <> struct Implementation<Ticket>
    {
        Mutex mutex;

        ConditionVariable completion;

        bool completed;

        Implementation() :
            completed(false)
        {
        }
    };

    Ticket::Ticket() :
        PrivateImplementationPattern<Ticket, Shared>(new Implementation<Ticket>)
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

    void
    TicketList::wait() const
    {
        CONTEXT("When waiting for the Tickets in a TicketList:");

        while (! _imp->tickets.empty())
        {
            std::tr1::shared_ptr<Implementation<Ticket> > ticket(_imp->tickets.front());
            Lock l(ticket->mutex);

            while (! ticket->completed)
            {
                ticket->completion.wait(ticket->mutex);
            }

            _imp->tickets.pop_front();
        }
    }
}
