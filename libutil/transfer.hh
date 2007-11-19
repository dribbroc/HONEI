/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#ifndef LIBUTIL_GUARD_TRANSFER_HH
#define LIBUTIL_GUARD_TRANSFER_HH 1

namespace intern
{
    extern const int sse_double_shuffle_patterns[2];
    extern const int sse_float_shuffle_patterns[4];
}

namespace honei
{
    template <typename DT_> inline int sse_shuffle_mask (unsigned long offset) __attribute__((always_inline));

    template <> inline int sse_shuffle_mask<double> (unsigned long offset)
    {
        return intern::sse_double_shuffle_patterns[offset];
    }
}

#endif
