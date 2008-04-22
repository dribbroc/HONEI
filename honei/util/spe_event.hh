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

#ifndef LIBUTIL_GUARD_SPE_EVENT_HH
#define LIBUTIL_GUARD_SPE_EVENT_HH 1

#include <honei/util/instantiation_policy.hh>
#include <honei/util/spe_manager.hh>

namespace honei
{
    class SPEEvent :
        public InstantiationPolicy<SPEEvent, NonCopyable>
    {
        private:
            /// Our event handler.
            spe_event_handler_ptr_t _handler;

            /// Our event interests.
            spe_event_unit_t _events;

            /// Our most recent event.
            spe_event_unit_t _recent_event;

        public:
            /// Constructor.
            SPEEvent(const SPE & spe, const unsigned int events);

            /// Destructor.
            ~SPEEvent();

            /// Wait until event happened.
            spe_event_unit_t * wait(int timeout = -1);
    };
}

#endif
