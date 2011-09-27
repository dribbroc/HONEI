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
#ifndef WOOLB3_GUARD_FORCE_HH
#define WOOLB3_GUARD_FORCE_HH 1



#include <honei/lbm/tags.hh>
#include <honei/la/dense_vector.hh>
#include <honei/woolb3/packed_grid3.hh>
#include <honei/util/profiler.hh>
#include <cmath>
#include <limits>

namespace honei
{
    template <typename Tag_, typename SourceScheme_>
    struct Force
    {
    };

    template <>
    struct Force<tags::CPU, lbm::lbm_source_schemes::BED_SLOPE>
    {
        template <typename DT_, unsigned long directions>
        static void value(PackedGrid3<DT_, directions> & pgrid, DT_ g, DT_ d_x, DT_ d_y, DT_ d_t, DT_ /*manning*/,
                bool inner)
        {
            unsigned long start(inner == true ? 0 : pgrid.grid.size() - pgrid.grid.inner_halo_size());
            unsigned long end(inner == true ? pgrid.grid.local_size() : pgrid.grid.size());

            DT_ force_multiplier(d_t / (6 * d_x * d_x / (d_t * d_t)));
            DT_ gravity_multiplier(-g);
            DT_ force_times_gravity(force_multiplier * gravity_multiplier);

            DT_ temp(0);
            DT_ * h(pgrid.h->elements());
            DT_ * b(pgrid.b->elements());

            unsigned long direction(1);
            DT_ distribution_x((*pgrid.distribution_x)[direction]);
            DT_ distribution_y((*pgrid.distribution_y)[direction]);
            unsigned long * neighbours(pgrid.neighbours[direction]->elements());
            unsigned long * neighbours_other(pgrid.neighbours[direction]->elements());
            DT_ * f_temp(pgrid.f_temp[direction]->elements());

            for (unsigned long i(start) ; i < end ; ++i)
            {
                if (neighbours[i] < pgrid.h->size())
                {
                    temp = force_times_gravity * distribution_x;
                    temp *= (h[i] + h[neighbours[i]]) / (DT_(2));
                    if (neighbours_other[i] < pgrid.h->size())
                        temp *= (b[neighbours_other[i]] - b[i]) / (d_x);
                    f_temp[i] += temp;
                }
            }

            direction = 2;
            distribution_x = (*pgrid.distribution_x)[direction];
            distribution_y = (*pgrid.distribution_y)[direction];
            neighbours = pgrid.neighbours[direction]->elements();
            neighbours_other = pgrid.neighbours[1]->elements();
            f_temp = pgrid.f_temp[direction]->elements();
            for (unsigned long i(start) ; i < end ; ++i)
            {
                if (neighbours[i] < pgrid.h->size())
                {
                    temp = force_times_gravity * distribution_x;
                    temp *= (h[i] + h[neighbours[i]]) / (DT_(2));
                    if (neighbours_other[i] < pgrid.h->size())
                        temp *= (b[neighbours_other[i]] - b[i]) / (d_x);
                    f_temp[i] += temp;
                }
            }

            neighbours_other = pgrid.neighbours[3]->elements();
            for (unsigned long i(start) ; i < end ; ++i)
            {
                if (neighbours[i] < pgrid.h->size())
                {
                    temp = force_times_gravity * distribution_y;
                    temp *= (h[i] + h[(*pgrid.neighbours[direction])[i]]) / (DT_(2));
                    if (neighbours_other[i] < pgrid.h->size())
                        temp *= (b[neighbours_other[i]] - b[i]) / (d_y);
                    f_temp[i] += temp;
                }
            }

            direction = 3;
            distribution_x = (*pgrid.distribution_x)[direction];
            distribution_y = (*pgrid.distribution_y)[direction];
            neighbours = pgrid.neighbours[direction]->elements();
            neighbours_other = pgrid.neighbours[direction]->elements();
            f_temp = pgrid.f_temp[direction]->elements();
            for (unsigned long i(start) ; i < end ; ++i)
            {
                if (neighbours[i] < pgrid.h->size())
                {
                    temp = force_times_gravity * distribution_y;
                    temp *= (h[i] + h[neighbours[i]]) / (DT_(2));
                    if (neighbours_other[i] < pgrid.h->size())
                        temp *= (b[neighbours_other[i]] - b[i]) / (d_y);
                    f_temp[i] += temp;
                }
            }

            direction = 4;
            distribution_x = (*pgrid.distribution_x)[direction];
            distribution_y = (*pgrid.distribution_y)[direction];
            neighbours = pgrid.neighbours[direction]->elements();
            neighbours_other = pgrid.neighbours[1]->elements();
            f_temp = pgrid.f_temp[direction]->elements();
            for (unsigned long i(start) ; i < end ; ++i)
            {
                if (neighbours[i] < pgrid.h->size())
                {
                    temp = force_times_gravity * distribution_x;
                    temp *= (h[i] + h[neighbours[i]]) / (DT_(2));
                    if (neighbours_other[i] < pgrid.h->size())
                        temp *= (b[neighbours_other[i]] - b[i]) / (d_x);
                    f_temp[i] += temp;
                }
            }

            neighbours_other = pgrid.neighbours[3]->elements();
            for (unsigned long i(start) ; i < end ; ++i)
            {
                if (neighbours[i] < pgrid.h->size())
                {
                    temp = force_times_gravity * distribution_y;
                    temp *= (h[i] + h[neighbours[i]]) / (DT_(2));
                    if (neighbours_other[i] < pgrid.h->size())
                        temp *= (b[neighbours_other[i]] - b[i]) / (d_y);
                    f_temp[i] += temp;
                }
            }

            direction = 5;
            distribution_x = (*pgrid.distribution_x)[direction];
            distribution_y = (*pgrid.distribution_y)[direction];
            neighbours = pgrid.neighbours[direction]->elements();
            neighbours_other = pgrid.neighbours[1]->elements();
            f_temp = pgrid.f_temp[direction]->elements();
            for (unsigned long i(start) ; i < end ; ++i)
            {
                if (neighbours[i] < pgrid.h->size())
                {
                    temp = force_times_gravity * distribution_x;
                    temp *= (h[i] + h[neighbours[i]]) / (DT_(2));
                    if (neighbours_other[i] < pgrid.h->size())
                        temp *= (b[neighbours_other[i]] - b[i]) / (d_x);
                    f_temp[i] += temp;
                }
            }

            direction = 6;
            distribution_x = (*pgrid.distribution_x)[direction];
            distribution_y = (*pgrid.distribution_y)[direction];
            neighbours = pgrid.neighbours[direction]->elements();
            neighbours_other = pgrid.neighbours[1]->elements();
            f_temp = pgrid.f_temp[direction]->elements();
            for (unsigned long i(start) ; i < end ; ++i)
            {
                if (neighbours[i] < pgrid.h->size())
                {
                    temp = force_times_gravity * distribution_x;
                    temp *= (h[i] + h[neighbours[i]]) / (DT_(2));
                    if (neighbours_other[i] < pgrid.h->size())
                        temp *= (b[neighbours_other[i]] - b[i]) / (d_x);
                    f_temp[i] += temp;
                }
            }

            neighbours_other = pgrid.neighbours[3]->elements();
            for (unsigned long i(start) ; i < end ; ++i)
            {
                if (neighbours[i] < pgrid.h->size())
                {
                    temp = force_times_gravity * distribution_y;
                    temp *= (h[i] + h[neighbours[i]]) / (DT_(2));
                    if (neighbours_other[i] < pgrid.h->size())
                        temp *= (b[neighbours_other[i]] - b[i]) / (d_y);
                    f_temp[i] += temp;
                }
            }

            direction = 7;
            distribution_x = (*pgrid.distribution_x)[direction];
            distribution_y = (*pgrid.distribution_y)[direction];
            neighbours = pgrid.neighbours[direction]->elements();
            neighbours_other = pgrid.neighbours[3]->elements();
            f_temp = pgrid.f_temp[direction]->elements();
            for (unsigned long i(start) ; i < end ; ++i)
            {
                if (neighbours[i] < pgrid.h->size())
                {
                    temp = force_times_gravity * distribution_y;
                    temp *= (h[i] + h[neighbours[i]]) / (DT_(2));
                    if (neighbours_other[i] < pgrid.h->size())
                        temp *= (b[neighbours_other[i]] - b[i]) / (d_y);
                    f_temp[i] += temp;
                }
            }

            direction = 8;
            distribution_x = (*pgrid.distribution_x)[direction];
            distribution_y = (*pgrid.distribution_y)[direction];
            neighbours = pgrid.neighbours[direction]->elements();
            neighbours_other = pgrid.neighbours[1]->elements();
            f_temp = pgrid.f_temp[direction]->elements();
            for (unsigned long i(start) ; i < end ; ++i)
            {
                if (neighbours[i] < pgrid.h->size())
                {
                    temp = force_times_gravity * distribution_x;
                    temp *= (h[i] + h[neighbours[i]]) / (DT_(2));
                    if (neighbours_other[i] < pgrid.h->size())
                        temp *= (b[neighbours_other[i]] - b[i]) / (d_x);
                    f_temp[i] += temp;
                }
            }

            neighbours_other = pgrid.neighbours[3]->elements();
            for (unsigned long i(start) ; i < end ; ++i)
            {
                if (neighbours[i] < pgrid.h->size())
                {
                    temp = force_times_gravity * distribution_y;
                    temp *= (h[i] + h[neighbours[i]]) / (DT_(2));
                    if (neighbours_other[i] < pgrid.h->size())
                        temp *= (b[neighbours_other[i]] - b[i]) / (d_y);
                    f_temp[i] += temp;
                }
            }
        }
    };

    template <>
    struct Force<tags::CPU, lbm::lbm_source_schemes::BED_FRICTION>
    {
        template <typename DT_, unsigned long directions>
        static void value(PackedGrid3<DT_, directions> & pgrid, DT_ g, DT_ d_x, DT_ /*d_y*/, DT_ d_t, DT_ manning,
                bool inner)
        {
            unsigned long start(inner == true ? 0 : pgrid.grid.size() - pgrid.grid.inner_halo_size());
            unsigned long end(inner == true ? pgrid.grid.local_size() : pgrid.grid.size());

            DT_ force_multiplier(d_t / (DT_(6) * d_x * d_x / (d_t * d_t)));
            force_multiplier *= g;
            DT_ manning2(manning * manning);
            DT_ * h(pgrid.h->elements());
            DT_ * u(pgrid.u->elements());
            DT_ * v(pgrid.v->elements());

            unsigned long direction(1);
            DT_ * f_temp((pgrid.f_temp[direction])->elements());
            DT_ distribution_x((*pgrid.distribution_x)[direction]);
            DT_ distribution_y((*pgrid.distribution_y)[direction]);
            for (unsigned long i(start) ; i < end ; ++i)
            {
                if( pow(h[i], DT_(1./3.)) > std::numeric_limits<DT_>::epsilon() || pow(h[i], DT_(1./3.)) < -std::numeric_limits<DT_>::epsilon())
                    f_temp[i] -= force_multiplier * distribution_x * manning2 *
                        u[i] * sqrt(u[i] * u[i] + v[i] * v[i]) / pow(h[i], DT_(1./3.));
            }

            direction = 2;
            f_temp = (pgrid.f_temp[direction])->elements();
            distribution_x = (*pgrid.distribution_x)[direction];
            distribution_y = (*pgrid.distribution_y)[direction];
            for (unsigned long i(start) ; i < end ; ++i)
            {
                if( pow(h[i], DT_(1./3.)) > std::numeric_limits<DT_>::epsilon() || pow(h[i], DT_(1./3.)) < -std::numeric_limits<DT_>::epsilon())
                    f_temp[i] -= force_multiplier * distribution_x * manning2 *
                        u[i] * sqrt(u[i] * u[i] + v[i] * v[i]) / pow(h[i], DT_(1./3.));
            }
            for (unsigned long i(start) ; i < end ; ++i)
            {
                if( pow(h[i], DT_(1./3.)) > std::numeric_limits<DT_>::epsilon() || pow(h[i], DT_(1./3.)) < -std::numeric_limits<DT_>::epsilon())
                    f_temp[i] -= force_multiplier * distribution_y * manning2 *
                        v[i] * sqrt(u[i] * u[i] + v[i] * v[i]) / pow(h[i], DT_(1./3.));
            }

            direction = 3;
            f_temp = (pgrid.f_temp[direction])->elements();
            distribution_x = (*pgrid.distribution_x)[direction];
            distribution_y = (*pgrid.distribution_y)[direction];
            for (unsigned long i(start) ; i < end ; ++i)
            {
                if( pow(h[i], DT_(1./3.)) > std::numeric_limits<DT_>::epsilon() || pow(h[i], DT_(1./3.)) < -std::numeric_limits<DT_>::epsilon())
                    f_temp[i] -= force_multiplier * distribution_y * manning2 *
                        v[i] * sqrt(u[i] * u[i] + v[i] * v[i]) / pow(h[i], DT_(1./3.));
            }

            direction = 4;
            f_temp = (pgrid.f_temp[direction])->elements();
            distribution_x = (*pgrid.distribution_x)[direction];
            distribution_y = (*pgrid.distribution_y)[direction];
            for (unsigned long i(start) ; i < end ; ++i)
            {
                if( pow(h[i], DT_(1./3.)) > std::numeric_limits<DT_>::epsilon() || pow(h[i], DT_(1./3.)) < -std::numeric_limits<DT_>::epsilon())
                    f_temp[i] -= force_multiplier * distribution_x * manning2 *
                        u[i] * sqrt(u[i] * u[i] + v[i] * v[i]) / pow(h[i], DT_(1./3.));
            }
            for (unsigned long i(start) ; i < end ; ++i)
            {
                if( pow(h[i], DT_(1./3.)) > std::numeric_limits<DT_>::epsilon() || pow(h[i], DT_(1./3.)) < -std::numeric_limits<DT_>::epsilon())
                    f_temp[i] -= force_multiplier * distribution_y * manning2 *
                        v[i] * sqrt(u[i] * u[i] + v[i] * v[i]) / pow(h[i], DT_(1./3.));
            }

            direction = 5;
            f_temp = (pgrid.f_temp[direction])->elements();
            distribution_x = (*pgrid.distribution_x)[direction];
            distribution_y = (*pgrid.distribution_y)[direction];
            for (unsigned long i(start) ; i < end ; ++i)
            {
                if( pow(h[i], DT_(1./3.)) > std::numeric_limits<DT_>::epsilon() || pow(h[i], DT_(1./3.)) < -std::numeric_limits<DT_>::epsilon())
                    f_temp[i] -= force_multiplier * distribution_x * manning2 *
                        u[i] * sqrt(u[i] * u[i] + v[i] * v[i]) / pow(h[i], DT_(1./3.));
            }

            direction = 6;
            f_temp = (pgrid.f_temp[direction])->elements();
            distribution_x = (*pgrid.distribution_x)[direction];
            distribution_y = (*pgrid.distribution_y)[direction];
            for (unsigned long i(start) ; i < end ; ++i)
            {
                if( pow(h[i], DT_(1./3.)) > std::numeric_limits<DT_>::epsilon() || pow(h[i], DT_(1./3.)) < -std::numeric_limits<DT_>::epsilon())
                    f_temp[i] -= force_multiplier * distribution_x * manning2 *
                        u[i] * sqrt(u[i] * u[i] + v[i] * v[i]) / pow(h[i], DT_(1./3.));
            }
            for (unsigned long i(start) ; i < end ; ++i)
            {
                if( pow(h[i], DT_(1./3.)) > std::numeric_limits<DT_>::epsilon() || pow(h[i], DT_(1./3.)) < -std::numeric_limits<DT_>::epsilon())
                    f_temp[i] -= force_multiplier * distribution_y * manning2 *
                        v[i] * sqrt(u[i] * u[i] + v[i] * v[i]) / pow(h[i], DT_(1./3.));
            }

            direction = 7;
            f_temp = (pgrid.f_temp[direction])->elements();
            distribution_x = (*pgrid.distribution_x)[direction];
            distribution_y = (*pgrid.distribution_y)[direction];
            for (unsigned long i(start) ; i < end ; ++i)
            {
                if( pow(h[i], DT_(1./3.)) > std::numeric_limits<DT_>::epsilon() || pow(h[i], DT_(1./3.)) < -std::numeric_limits<DT_>::epsilon())
                    f_temp[i] -= force_multiplier * distribution_y * manning2 *
                        v[i] * sqrt(u[i] * u[i] + v[i] * v[i]) / pow(h[i], DT_(1./3.));
            }

            direction = 8;
            f_temp = (pgrid.f_temp[direction])->elements();
            distribution_x = (*pgrid.distribution_x)[direction];
            distribution_y = (*pgrid.distribution_y)[direction];
            for (unsigned long i(start) ; i < end ; ++i)
            {
                if( pow(h[i], DT_(1./3.)) > std::numeric_limits<DT_>::epsilon() || pow(h[i], DT_(1./3.)) < -std::numeric_limits<DT_>::epsilon())
                    f_temp[i] -= force_multiplier * distribution_x * manning2 *
                        u[i] * sqrt(u[i] * u[i] + v[i] * v[i]) / pow(h[i], DT_(1./3.));
            }
            for (unsigned long i(start) ; i < end ; ++i)
            {
                if( pow(h[i], DT_(1./3.)) > std::numeric_limits<DT_>::epsilon() || pow(h[i], DT_(1./3.)) < -std::numeric_limits<DT_>::epsilon())
                    f_temp[i] -= force_multiplier * distribution_y * manning2 *
                        v[i] * sqrt(u[i] * u[i] + v[i] * v[i]) / pow(h[i], DT_(1./3.));
            }
        }
    };

    template <typename Tag_>
        struct Force<Tag_, lbm::lbm_source_schemes::BED_FULL>
        {
            template <typename DT_, unsigned long directions>
                static void value(PackedGrid3<DT_, directions> & pgrid, DT_ g, DT_ d_x, DT_ d_y, DT_ d_t, DT_ manning,
                        bool inner)
                {
                    PROFILER_START("Force BED_FULL");
                    Force<Tag_, lbm::lbm_source_schemes::BED_SLOPE>::value(pgrid, g, d_x, d_y, d_t, manning, inner);
                    Force<Tag_, lbm::lbm_source_schemes::BED_FRICTION>::value(pgrid, g, d_x, d_y, d_t, manning, inner);
                    PROFILER_STOP("Force BED_FULL");
                }
        };
}
#endif
