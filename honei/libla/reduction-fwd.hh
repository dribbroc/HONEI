/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * This file is part of the LA C++ library. LibLa is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibLa is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef LIBLA_GUARD_REDUCTION_FWD_HH
#define LIBLA_GUARD_REDUCTION_FWD_HH 1

namespace honei
{
    /**
     * ReductionType is a template tag parameter for the Reduction class
     * template. It governs the type of the (mathematical) reduction that shall
     * be computed.
     *
     * \ingroup grplaoperations
     */
    enum ReductionType
    {
        rt_sum = 0,
        rt_max,
        rt_min
    };
}

#endif
