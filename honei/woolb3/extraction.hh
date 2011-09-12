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
#ifndef WOOLB3_GUARD_EXTRACTION_HH
#define WOOLB3_GUARD_EXTRACTION_HH 1



#include <honei/lbm/tags.hh>
#include <honei/la/dense_vector.hh>
#include <honei/woolb3/packed_grid.hh>
#include <cmath>

namespace honei
{
    template <typename Tag_>
    struct Extraction
    {
    };

    template <>
    struct Extraction<tags::CPU>
    {
        template <typename DT_, unsigned long directions>
        static void value(PackedGrid<DT_, directions> & pgrid)
        {
            for (unsigned long direction(0) ; direction < directions ; ++direction)
            {
                DenseVector<DT_> swap = *pgrid.f[direction];
                pgrid.f[direction].reset(new DenseVector<DT_>(*pgrid.f_temp[direction]));
                pgrid.f_temp[direction].reset(new DenseVector<DT_>(swap));
            }

            for (unsigned long i(0) ; i < pgrid.h->size() ; ++i)
            {
                (*pgrid.h)[i] = (*pgrid.f[0])[i] +
                    (*pgrid.f[1])[i] +
                    (*pgrid.f[2])[i] +
                    (*pgrid.f[3])[i] +
                    (*pgrid.f[4])[i] +
                    (*pgrid.f[5])[i] +
                    (*pgrid.f[6])[i] +
                    (*pgrid.f[7])[i] +
                    (*pgrid.f[8])[i];

            (*pgrid.u)[i] = ((*pgrid.distribution_x)[0] * (*pgrid.f[0])[i] +
                    (*pgrid.distribution_x)[1] * (*pgrid.f[1])[i] +
                    (*pgrid.distribution_x)[2] * (*pgrid.f[2])[i] +
                    (*pgrid.distribution_x)[3] * (*pgrid.f[3])[i] +
                    (*pgrid.distribution_x)[4] * (*pgrid.f[4])[i] +
                    (*pgrid.distribution_x)[5] * (*pgrid.f[5])[i] +
                    (*pgrid.distribution_x)[6] * (*pgrid.f[6])[i] +
                    (*pgrid.distribution_x)[7] * (*pgrid.f[7])[i] +
                    (*pgrid.distribution_x)[8] * (*pgrid.f[8])[i]) / (*pgrid.h)[i];

            (*pgrid.v)[i] = ((*pgrid.distribution_y)[0] * (*pgrid.f[0])[i] +
                    (*pgrid.distribution_y)[1] * (*pgrid.f[1])[i] +
                    (*pgrid.distribution_y)[2] * (*pgrid.f[2])[i] +
                    (*pgrid.distribution_y)[3] * (*pgrid.f[3])[i] +
                    (*pgrid.distribution_y)[4] * (*pgrid.f[4])[i] +
                    (*pgrid.distribution_y)[5] * (*pgrid.f[5])[i] +
                    (*pgrid.distribution_y)[6] * (*pgrid.f[6])[i] +
                    (*pgrid.distribution_y)[7] * (*pgrid.f[7])[i] +
                    (*pgrid.distribution_y)[8] * (*pgrid.f[8])[i]) / (*pgrid.h)[i];
            }
        }
    };
}
#endif
