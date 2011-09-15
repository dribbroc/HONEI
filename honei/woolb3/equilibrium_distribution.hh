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
#ifndef WOOLB3_GUARD_EQUILIBRIUM_DISTRIBUTION_HH
#define WOOLB3_GUARD_EQUILIBRIUM_DISTRIBUTION_HH 1



#include <honei/lbm/tags.hh>
#include <honei/la/dense_vector.hh>
#include <honei/woolb3/packed_grid3.hh>
#include <cmath>

namespace honei
{
    template <typename Tag_>
    struct EquilibriumDistribution
    {
    };

    template <>
    struct EquilibriumDistribution<tags::CPU>
    {
        template <typename DT_, unsigned long directions>
        static void value(PackedGrid3<DT_, directions> & pgrid, DT_ g, DT_ e, unsigned long start = 0, unsigned long end = 0)
        {
            if (end == 0)
                end = pgrid.h->size();
            unsigned long direction(0);

            DT_ * f_eq(pgrid.f_eq[direction]->elements());
            DT_ * const h(pgrid.h->elements());
            DT_ * const u(pgrid.u->elements());
            DT_ * const v(pgrid.v->elements());
            DT_ d_x;
            DT_ d_y;

            for (unsigned long i(start) ; i < end ; ++i)
            {
                f_eq[i] = h[i] - ((DT_(5.) * g * h[i] * h[i]) / (DT_(6.) * e * e)) -
                    ((DT_(2.) * h[i]) /(DT_(3.) * e * e) * (u[i] * u[i] + v[i] * v[i]));
            }

            for (unsigned long direction(1) ; direction < directions ; direction+=2)
            {
                f_eq = pgrid.f_eq[direction]->elements();
                d_x = (*pgrid.distribution_x)[direction];
                d_y = (*pgrid.distribution_y)[direction];
                for (unsigned long i(start) ; i < end ; ++i)
                {
                    f_eq[i] = ((g * h[i] * h[i]) /(DT_(6.) * e * e)) +
                        ((h[i] / (DT_(3.) * e * e)) * (d_x * u[i] + d_y * v[i])) +
                        ((h[i] / (DT_(2.) * e * e * e * e)) * (d_x * u[i] * d_x * u[i] + DT_(2.) * d_x * u[i] * d_y * v[i] + d_y * v[i] * d_y * v[i])) -
                        ((h[i] / (DT_(6.) * e * e)) * (u[i] * u[i] + v[i] * v[i]));
                }
            }

            for (unsigned long direction(2) ; direction < directions ; direction+=2)
            {
                f_eq = pgrid.f_eq[direction]->elements();
                d_x = (*pgrid.distribution_x)[direction];
                d_y = (*pgrid.distribution_y)[direction];
                for (unsigned long i(start) ; i < end ; ++i)
                {
                            f_eq[i] = ((g * h[i] * h[i]) /(DT_(24.) * e * e)) +
                                          ((h[i] / (DT_(12.) * e * e)) * (d_x * u[i] + d_y * v[i])) +
                                          ((h[i] / (DT_(8.) * e * e * e * e)) * (d_x * u[i] * d_x * u[i] + DT_(2.) * d_x * u[i] * d_y * v[i] + d_y * v[i] * d_y * v[i])) -
                                          ((h[i] / (DT_(24.) * e * e)) * (u[i] * u[i] + v[i] * v[i]));
                }
            }
        }
    };
}
#endif
