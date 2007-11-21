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

#ifndef CELL_GUARD_UTIL_TRANSFER_HH
#define CELL_GUARD_UTIL_TRANSFER_HH 1

#include <cell/cell.hh>

#include <spu_intrinsics.h>

namespace intern
{
    extern const vector unsigned char extract_patterns[4];
}

inline unsigned multiple_of_sixteen(unsigned u) __attribute__((always_inline));

unsigned multiple_of_sixteen(unsigned u)
{
    unsigned r(u & 0xF);

    return (r > 0) ? u + 16 - r : u;
}

template <typename DT_> inline void extract(DT_ & first, const DT_ & second, unsigned offset) __attribute__((always_inline));

template <> inline void extract<vector float>(vector float & first, const vector float & second, unsigned offset)
{
    first = spu_shuffle(first, second, intern::extract_patterns[offset]);
}

template <> inline void extract<vector double>(vector double & first, const vector double & second, unsigned offset)
{
    /// \todo Should we use an own lookup pattern for double version to avoid multiplication?
    first = spu_shuffle(first, second, intern::extract_patterns[offset * 2]);
}

template <typename DT_> inline void fill(void * address, unsigned long size, DT_ value) __attribute__((always_inline));

template <> inline void fill<float>(void * address, unsigned long size, float value)
{
    honei::cell::Pointer<float> p = { address };
    vector float v(spu_splats(value));

    unsigned i(0);
    for ( ; i < size / sizeof(float) ; ++i)
    {
        p.vectorised[i] = v;
    }

    for (unsigned j(0) ; j < size % sizeof(float) ; j++)
    {
        p.typed[i * 4 + j] = value;
    }
}

template <> inline void fill<double>(void * address, unsigned long size, double value)
{
    honei::cell::Pointer<double> p = { address };
    vector double v(spu_splats(value));

    unsigned i(0);
    for ( ; i < size / sizeof(double) ; ++i)
    {
        p.vectorised[i] = v;
    }

    for (unsigned j(0) ; j < size % sizeof(double) ; j++)
    {
        p.typed[i * 2 + j] = value;
    }
}
#endif

