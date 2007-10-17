/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#include <libutil/assertion.hh>
#include <libutil/spe_event.hh>

#include <cstring>
#include <cerrno>

namespace honei
{
    SPEEvent::SPEEvent(const SPE & spe, const unsigned events) :
        _handler(spe_event_handler_create())
    {
        CONTEXT("When creating SPEEvent for SPE #" + stringify(spe.id()) + ":");

        _events.events = events;
        _events.spe = spe.context();

        if (0 != spe_event_handler_register(_handler, &_events))
        {
            /// \todo check for ENOTSUP
        }
    }

    SPEEvent::~SPEEvent()
    {
        int result(0);

        if (0 != spe_event_handler_deregister(_handler, &_events))
        {
            /// \todo check for ENOTSUP
        }

        do
        {
            result = spe_event_handler_destroy(_handler);
        }
        while ((-1 == result) && (EAGAIN == errno));

        ASSERT(0 == result, "spe_event_handler_destroy returned '" + stringify(strerror(errno)) + "'!");
    }

    spe_event_unit_t *
    SPEEvent::wait(int timeout)
    {
        memset(&_recent_event, sizeof(spe_event_unit_t), 0);

        spe_event_wait(_handler, &_recent_event, 1, timeout);

        return &_recent_event;
    }
}
