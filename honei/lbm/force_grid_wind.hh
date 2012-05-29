/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2012 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the LBM C++ library. LBM is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LBM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef LBM_GUARD_FORCE_GRID_WIND_HH
#define LBM_GUARD_FORCE_GRID_WIND_HH 1

#include<honei/lbm/force_grid.hh>

using namespace honei::lbm;

namespace honei
{
    template <>
    struct ForceGrid<tags::CPU::Generic, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::WIND_FRICTION>
    {
        template<typename DT1_>
            static void value(PackedGridInfo<D2Q9> & info,
                              PackedGridData<D2Q9, DT1_> & data,
                              DT1_ d_x,
                              HONEI_UNUSED DT1_ d_y,
                              DT1_ d_t,
                              HONEI_UNUSED DT1_ manning,
                              DT1_ w_x = 0.,
                              DT1_ w_y = 0.,
                              DT1_ absolute_pressure = 101325.,     // Pa (standard for dry air at 20 degrees celsius)
                              DT1_ specific_gas_constant = 287.058, // J/(kg . K) (dry air)
                              DT1_ gas_temperature = 20)            // degrees celsius
            {
                CONTEXT("When computing wind friction force term:");

                info.limits->lock(lm_read_only);

                data.distribution_x->lock(lm_read_only);
                data.distribution_y->lock(lm_read_only);

                data.f_temp_1->lock(lm_read_and_write);
                data.f_temp_2->lock(lm_read_and_write);
                data.f_temp_3->lock(lm_read_and_write);
                data.f_temp_4->lock(lm_read_and_write);
                data.f_temp_5->lock(lm_read_and_write);
                data.f_temp_6->lock(lm_read_and_write);
                data.f_temp_7->lock(lm_read_and_write);
                data.f_temp_8->lock(lm_read_and_write);

                const unsigned long * const limits(info.limits->elements());
                const DT1_ * const distribution_x(data.distribution_x->elements());
                const DT1_ * const distribution_y(data.distribution_y->elements());

                DT1_ * f_temp_1(data.f_temp_1->elements());
                DT1_ * f_temp_2(data.f_temp_2->elements());
                DT1_ * f_temp_3(data.f_temp_3->elements());
                DT1_ * f_temp_4(data.f_temp_4->elements());
                DT1_ * f_temp_5(data.f_temp_5->elements());
                DT1_ * f_temp_6(data.f_temp_6->elements());
                DT1_ * f_temp_7(data.f_temp_7->elements());
                DT1_ * f_temp_8(data.f_temp_8->elements());

                //precompute constants
                DT1_ force_multiplier(d_t / (6 * d_x * d_x / (d_t * d_t)));
                DT1_ p_alpha(absolute_pressure / (specific_gas_constant * gas_temperature));
                DT1_ veloc_abs(sqrt(w_x * w_x + w_y * w_y));
                DT1_ wind_veloc_contribution(DT1_(0.75) + DT1_(0.067) * veloc_abs);
                DT1_ total(force_multiplier * p_alpha * wind_veloc_contribution * DT1_(0.001));

                //std::cout << "Total: " << total << std::endl;

                for(unsigned long i((limits)[0]); i < (limits)[info.limits->size() - 1]; ++i)
                {
                        (f_temp_1)[i] += distribution_x[1] * total * w_x * veloc_abs + distribution_y[1] * total * w_y * veloc_abs;
                        (f_temp_2)[i] += distribution_x[2] * total * w_x * veloc_abs + distribution_y[2] * total * w_y * veloc_abs;
                        (f_temp_3)[i] += distribution_x[3] * total * w_x * veloc_abs + distribution_y[3] * total * w_y * veloc_abs;
                        (f_temp_4)[i] += distribution_x[4] * total * w_x * veloc_abs + distribution_y[4] * total * w_y * veloc_abs;
                        (f_temp_5)[i] += distribution_x[5] * total * w_x * veloc_abs + distribution_y[5] * total * w_y * veloc_abs;
                        (f_temp_6)[i] += distribution_x[6] * total * w_x * veloc_abs + distribution_y[6] * total * w_y * veloc_abs;
                        (f_temp_7)[i] += distribution_x[7] * total * w_x * veloc_abs + distribution_y[7] * total * w_y * veloc_abs;
                        (f_temp_8)[i] += distribution_x[8] * total * w_x * veloc_abs + distribution_y[8] * total * w_y * veloc_abs;
                }

                info.limits->unlock(lm_read_only);

                data.distribution_x->unlock(lm_read_only);
                data.distribution_y->unlock(lm_read_only);

                data.f_temp_1->unlock(lm_read_and_write);
                data.f_temp_2->unlock(lm_read_and_write);
                data.f_temp_3->unlock(lm_read_and_write);
                data.f_temp_4->unlock(lm_read_and_write);
                data.f_temp_5->unlock(lm_read_and_write);
                data.f_temp_6->unlock(lm_read_and_write);
                data.f_temp_7->unlock(lm_read_and_write);
                data.f_temp_8->unlock(lm_read_and_write);
            }
    };

    ///todo unify interface WIND VELOC IS HARD WIRED
    template <typename Tag_>
        struct ForceGrid<Tag_, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::FULL>
        {
            template<typename DT1_, typename DT2_>
                static void value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT1_> & data, DT2_ g, DT2_ d_x, DT2_ d_y, DT2_ d_t, DT2_ manning, DT2_ w_x = 38.875, DT2_ w_y = -38.875)
                {
                    ForceGrid<Tag_, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_SLOPE>::value(info, data, g, d_x, d_y, d_t, manning);
                    ForceGrid<Tag_, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_FRICTION>::value(info, data, g, d_x, d_y, d_t, manning);
                    ForceGrid<Tag_, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::WIND_FRICTION>::value(info, data, g, d_x, d_y, d_t, manning, w_x, w_y);
                }
        };
}
#endif
