/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2011 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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


#pragma once
#ifndef WOOLB3_GUARD_COLLIDE_STREAM_HH
#define WOOLB3_GUARD_COLLIDE_STREAM_HH 1



#include <honei/la/dense_vector.hh>
#include <honei/woolb3/packed_grid3.hh>
#include <honei/util/profiler.hh>
#include <cmath>

namespace honei
{
    template <typename Tag_>
    struct CollideStream
    {
    };

    template <>
    struct CollideStream<tags::CPU>
    {
        template <typename DT_, unsigned long directions>
        static void value(PackedGrid3<DT_, directions> & pgrid, DT_ tau, bool inner)
        {
            PROFILER_START("CollStream");

            for (unsigned long direction(0) ; direction < directions ; ++direction)
            {
                const DT_ * const f(pgrid.f[direction]->elements());
                const DT_ * const f_eq(pgrid.f_eq[direction]->elements());
                DT_ * f_temp2(pgrid.f_temp2[direction]->elements());
                /*unsigned long * dir_index;
                unsigned long * dir;
                unsigned long size;
                if (inner)
                {
                    dir_index = pgrid.dir_index_inner[direction]->elements();
                    dir = pgrid.dir_inner[direction]->elements();
                    size = pgrid.dir_index_inner[direction]->size() - 1;
                }
                else
                {
                    dir_index = pgrid.dir_index_outer[direction]->elements();
                    dir = pgrid.dir_outer[direction]->elements();
                    size = pgrid.dir_index_outer[direction]->size() - 1;
                }
                for (unsigned long begin(0), half(0) ; begin < size ; begin+=2, ++half)
                {
                    const unsigned long end_dir(dir_index[begin + 1]);
                    for (unsigned long i(dir_index[begin]), offset(0) ; i < end_dir ; ++i, ++offset)
                    {
                        f_temp2[dir[half] + offset] = f[i] - (f[i] - f_eq[i])/tau;
                    }
                }
*/
                unsigned long start(inner == true ? 0 : pgrid.grid.size() - pgrid.grid.inner_halo_size());
                unsigned long end(inner == true ? pgrid.grid.local_size() : pgrid.grid.size());
                const unsigned long * neighbours(pgrid.neighbours[direction]->elements());
                for (unsigned long i(start) ; i < end ; ++i)
                {
                    if (neighbours[i] < pgrid.h->size())
                        f_temp2[neighbours[i]] = f[i] - (f[i] - f_eq[i])/tau;
                }
            }
            PROFILER_STOP("CollStream");
        }
    };
}
#endif
