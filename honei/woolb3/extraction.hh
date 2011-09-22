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



#include <honei/la/dense_vector.hh>
#include <honei/woolb3/packed_grid3.hh>
#include <honei/lbm/lbm_limiter.hh>
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
        static void value(PackedGrid3<DT_, directions> & pgrid, unsigned long start = 0, unsigned long end = 0)
        {
            if (end == 0)
                end = pgrid.h->size();

            DT_ * f[directions];
            for (unsigned long dir(0) ; dir < directions ; ++dir)
                f[dir] = pgrid.f[dir]->elements();
            DT_ * f_temp[directions];
            for (unsigned long dir(0) ; dir < directions ; ++dir)
                f_temp[dir] = pgrid.f_temp[dir]->elements();

            for (unsigned long direction(0) ; direction < directions ; ++direction)
            {
                for (unsigned long i(start) ; i < end ; ++i)
                {
                    DT_ t = (f[direction])[i];
                    (f[direction])[i] = (f_temp[direction])[i];
                    (f_temp[direction])[i] = t;
                }
            }

            DT_ * h(pgrid.h->elements());
            DT_ * u(pgrid.u->elements());
            DT_ * v(pgrid.v->elements());

            DT_ lax_upper(10e-5);
            DT_ lax_lower(-lax_upper);

            for (unsigned long i(start) ; i < end ; ++i)
            {
                (h)[i] = (f[0])[i] +
                    (f[1])[i] +
                    (f[2])[i] +
                    (f[3])[i] +
                    (f[4])[i] +
                    (f[5])[i] +
                    (f[6])[i] +
                    (f[7])[i] +
                    (f[8])[i];

                //if (fabs(h[i]) > 10e-5)
                if(h[i] < lax_lower || h[i] > lax_upper)
                {
                    (u)[i] = ((*pgrid.distribution_x)[0] * (f[0])[i] +
                            (*pgrid.distribution_x)[1] * (f[1])[i] +
                            (*pgrid.distribution_x)[2] * (f[2])[i] +
                            (*pgrid.distribution_x)[3] * (f[3])[i] +
                            (*pgrid.distribution_x)[4] * (f[4])[i] +
                            (*pgrid.distribution_x)[5] * (f[5])[i] +
                            (*pgrid.distribution_x)[6] * (f[6])[i] +
                            (*pgrid.distribution_x)[7] * (f[7])[i] +
                            (*pgrid.distribution_x)[8] * (f[8])[i]) / h[i];

                    (v)[i] = ((*pgrid.distribution_y)[0] * (f[0])[i] +
                            (*pgrid.distribution_y)[1] * (f[1])[i] +
                            (*pgrid.distribution_y)[2] * (f[2])[i] +
                            (*pgrid.distribution_y)[3] * (f[3])[i] +
                            (*pgrid.distribution_y)[4] * (f[4])[i] +
                            (*pgrid.distribution_y)[5] * (f[5])[i] +
                            (*pgrid.distribution_y)[6] * (f[6])[i] +
                            (*pgrid.distribution_y)[7] * (f[7])[i] +
                            (*pgrid.distribution_y)[8] * (f[8])[i]) / h[i];
                }
                else
                {
                    h[i] = DT_(0);
                    u[i] = DT_(0);
                    v[i] = DT_(0);
                }
                h[i] = lbm::MinModLimiter<tags::CPU>::value(h[i]);
            }
        }
    };
}
#endif
