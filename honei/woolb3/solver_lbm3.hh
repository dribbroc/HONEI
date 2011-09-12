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
#ifndef WOOLB3_GUARD_SOLVER_LBM3_HH
#define WOOLB3_GUARD_SOLVER_LBM3_HH 1



#include <honei/lbm/tags.hh>
#include <honei/la/dense_vector.hh>
#include <honei/woolb3/grid3.hh>
#include <honei/woolb3/packed_grid3.hh>
#include <honei/woolb3/equilibrium_distribution.hh>
#include <honei/woolb3/collide_stream.hh>
#include <honei/woolb3/extraction.hh>
#include <honei/woolb3/update_velocity_directions.hh>
#include <cmath>

namespace honei
{
    template <typename Tag_, typename DT_, unsigned long directions>
    struct SolverLBM3
    {
        private:
                DT_ _relaxation_time, _delta_x, _delta_y, _delta_t;

                Grid3<DT_, directions> _grid;
                PackedGrid3<DT_, directions> _pgrid;

                DT_ _n_alpha, _e, _gravity, _pi, _e_squared;

        public:
                SolverLBM3(Grid3<DT_, directions> & grid, PackedGrid3<DT_, directions> & pgrid, DT_ dx, DT_ dy, DT_ dt, DT_ rel_time) :
                    _relaxation_time(rel_time),
                    _delta_x(dx),
                    _delta_y(dy),
                    _delta_t(dt),
                    _grid(grid),
                    _pgrid(pgrid),
                    _n_alpha(DT_(6.)),
                    _gravity(9.80665),
                    _pi(3.14159265)
                {
                    _e = _delta_x / _delta_t;
                    _e_squared = _e * _e;
                }

                void do_preprocessing()
                {
                    (*_pgrid.distribution_x)[0] = DT_(0.);
                    (*_pgrid.distribution_x)[1] = DT_(_e * cos(DT_(0.)));
                    (*_pgrid.distribution_x)[2] = DT_(sqrt(DT_(2.)) * _e * cos(_pi / DT_(4.)));
                    (*_pgrid.distribution_x)[3] = DT_(_e * cos(_pi / DT_(2.)));
                    (*_pgrid.distribution_x)[4] = DT_(sqrt(DT_(2.)) * _e * cos(DT_(3.) * _pi / DT_(4.)));
                    (*_pgrid.distribution_x)[5] = DT_(_e * cos(_pi));
                    (*_pgrid.distribution_x)[6] = DT_(sqrt(DT_(2.)) * _e * cos(DT_(5.) * _pi / DT_(4.)));
                    (*_pgrid.distribution_x)[7] = DT_(_e * cos(DT_(3.) * _pi / DT_(2.)));
                    (*_pgrid.distribution_x)[8] = DT_(sqrt(DT_(2.)) * _e * cos(DT_(7.) * _pi / DT_(4.)));
                    (*_pgrid.distribution_y)[0] = DT_(0.);
                    (*_pgrid.distribution_y)[1] = DT_(_e * sin(DT_(0.)));
                    (*_pgrid.distribution_y)[2] = DT_(sqrt(DT_(2.)) * _e * sin(_pi / DT_(4.)));
                    (*_pgrid.distribution_y)[3] = DT_(_e * sin(_pi / DT_(2.)));
                    (*_pgrid.distribution_y)[4] = DT_(sqrt(DT_(2.)) * _e * sin(DT_(3.) * _pi / DT_(4.)));
                    (*_pgrid.distribution_y)[5] = DT_(_e * sin(_pi));
                    (*_pgrid.distribution_y)[6] = DT_(sqrt(DT_(2.)) * _e * sin(DT_(5.) * _pi / DT_(4.)));
                    (*_pgrid.distribution_y)[7] = DT_(_e * sin(DT_(3.) * _pi / DT_(2.)));
                    (*_pgrid.distribution_y)[8] = DT_(sqrt(DT_(2.)) * _e * sin(DT_(7.) * _pi / DT_(4.)));

                    ///Compute initial equilibrium distribution:
                    EquilibriumDistribution<Tag_>::value(_pgrid, _gravity, _e);

                    for (unsigned long direction(0) ; direction < directions ; ++direction)
                        _pgrid.f[direction].reset(new DenseVector<DT_>(_pgrid.f_eq[direction]->copy()));

                    CollideStream<Tag_>::value(_pgrid, _relaxation_time);
                }

                void solve()
                {
                    UpdateVelocityDirections<Tag_>::value(_grid, _pgrid);

                    Extraction<Tag_>::value(_pgrid);

                    EquilibriumDistribution<Tag_>::value(_pgrid, _gravity, _e);

                    CollideStream<Tag_>::value(_pgrid, _relaxation_time);
                }
    };
}
#endif
