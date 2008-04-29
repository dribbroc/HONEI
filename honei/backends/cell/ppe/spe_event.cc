/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * This file is part of the Utility C++ library. LibUtil is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibUtil is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/backends/cell/ppe/spe_error.hh>
#include <honei/backends/cell/ppe/spe_event.hh>
#include <honei/util/assertion.hh>
#include <honei/util/exception.hh>
#include <honei/util/instantiation_policy-impl.hh>
#include <honei/util/stringify.hh>

#include <cerrno>
#include <cstring>

namespace honei
{
    SPEEvent::SPEEvent(const SPE & spe, const unsigned events) :
        _handler(spe_event_handler_create())
    {
        CONTEXT("When creating SPEEvent for SPE #" + stringify(spe.id()) + ":");

        if (0 == _handler)
        {
            throw SPEError("spe_event_handler_create", errno);
        }

        _events.events = events;
        _events.spe = spe.context();

        unsigned retval;
        if (0 != (retval = spe_event_handler_register(_handler, &_events)))
        {
            throw SPEError("spe_event_handler_register", retval);
        }
    }

    SPEEvent::~SPEEvent()
    {
        int result(0);

        unsigned retval;
        if (0 != (retval = spe_event_handler_deregister(_handler, &_events)))
        {
            throw SPEError("spe_event_handler_deregister", retval);
        }

        do
        {
            result = spe_event_handler_destroy(_handler);
        }
        while ((-1 == result) && (EAGAIN == errno));

        if (-1 == result)
            throw SPEError("spe_event_handler_destroy", errno);
    }

    spe_event_unit_t *
    SPEEvent::wait(int timeout)
    {
        memset(&_recent_event, sizeof(spe_event_unit_t), 0);

        spe_event_wait(_handler, &_recent_event, 1, timeout);

        return &_recent_event;
    }
}
