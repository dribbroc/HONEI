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

#include <honei/backends/multicore/multicore_ticket.hh>

using namespace honei::mc;

MultiCoreTicket::MultiCoreTicket(const unsigned tid) :
        id(counter),
        thread_id(tid)
{
    ++counter;
}

MultiCoreTicket::~MultiCoreTicket()
{
}

unsigned
MultiCoreTicket::uid() const
{
    return id;
}

unsigned &
MultiCoreTicket::tid() 
{
    return thread_id;
}

unsigned MultiCoreTicket::counter(0);
