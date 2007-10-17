/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#ifndef LIBLA_CELL_GUARD_CELL_UTILITY_HH
#define LIBLA_CELL_GUARD_CELL_UTILITY_HH 1

inline unsigned multiple_of_sixteen(unsigned u) __attribute__((always_inline));

unsigned multiple_of_sixteen(unsigned u)
{
    unsigned r(u & 0xF);

    return (r > 0) ? u + 16 - r : u;
}

#endif
