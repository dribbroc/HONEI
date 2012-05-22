/* vim: set number sw=4 sts=4 et nofoldenable : */

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


#pragma once

#ifndef LBM_GUARD_FORCE_GRID_POLLUTANT_HH
#define LBM_GUARD_FORCE_GRID_POLLUTANT_HH 1

#include "config.h"
namespace honei
{
   template <typename Tag_,
              typename App_>
    struct ForceGridPollutant
    {
    };

    template <>
    struct ForceGridPollutant<tags::CPU::Generic, lbm_applications::LABSWE>
    {
        template<typename DT1_, typename DT2_>
            static void value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT1_> & data_flow, PackedGridData<D2Q9, DT1_> & data_poll, DT1_ dt, DT1_ k, DT1_ s_0)
            {
                CONTEXT("When computing force term for pollutant transport:");

                info.limits->lock(lm_read_only);

                data_flow.h->lock(lm_read_only);
                data_flow.u->lock(lm_read_only);
                data_flow.v->lock(lm_read_only);
                data_poll.h->lock(lm_read_only);

                data_poll.f_temp_0->lock(lm_write_only);
                data_poll.f_temp_1->lock(lm_write_only);
                data_poll.f_temp_2->lock(lm_write_only);
                data_poll.f_temp_3->lock(lm_write_only);
                data_poll.f_temp_4->lock(lm_write_only);
                data_poll.f_temp_5->lock(lm_write_only);
                data_poll.f_temp_6->lock(lm_write_only);
                data_poll.f_temp_7->lock(lm_write_only);
                data_poll.f_temp_8->lock(lm_write_only);

                const unsigned long * const limits(info.limits->elements());
                const DT1_ * const h(data_flow.h->elements());
                const DT1_ * const c(data_poll.h->elements());
                const DT1_ * const u(data_flow.u->elements());
                const DT1_ * const v(data_flow.v->elements());

                DT1_ * f_temp_0(data_poll.f_temp_0->elements());
                DT1_ * f_temp_1(data_poll.f_temp_1->elements());
                DT1_ * f_temp_2(data_poll.f_temp_2->elements());
                DT1_ * f_temp_3(data_poll.f_temp_3->elements());
                DT1_ * f_temp_4(data_poll.f_temp_4->elements());
                DT1_ * f_temp_5(data_poll.f_temp_5->elements());
                DT1_ * f_temp_6(data_poll.f_temp_6->elements());
                DT1_ * f_temp_7(data_poll.f_temp_7->elements());
                DT1_ * f_temp_8(data_poll.f_temp_8->elements());

                DT1_ t1, t2;

                DT1_ w_0(4./9.);
                DT1_ w_odd(1./9.);
                DT1_ w_even(1./36.);
                DT1_ const_multiplier_0(-dt * w_0);
                DT1_ const_multiplier_odd(-dt * w_odd);
                DT1_ const_multiplier_even(-dt * w_even);

                for(unsigned long i((limits)[0]); i < (limits)[info.limits->size() - 1]; ++i)
                {
                    t1 = k * h[i] * c[i];
                    t2 = s_0 * h[i];
                    (f_temp_0)[i] += const_multiplier_0 * (t1 + t2);
                    (f_temp_1)[i] += const_multiplier_odd * (t1 + t2);
                    (f_temp_3)[i] += const_multiplier_odd * (t1 + t2);
                    (f_temp_5)[i] += const_multiplier_odd * (t1 + t2);
                    (f_temp_7)[i] += const_multiplier_odd * (t1 + t2);
                    (f_temp_2)[i] += const_multiplier_even * (t1 + t2);
                    (f_temp_4)[i] += const_multiplier_even * (t1 + t2);
                    (f_temp_6)[i] += const_multiplier_even * (t1 + t2);
                    (f_temp_8)[i] += const_multiplier_even * (t1 + t2);
                }

                info.limits->unlock(lm_read_only);

                data_flow.h->unlock(lm_read_only);
                data_flow.u->unlock(lm_read_only);
                data_flow.v->unlock(lm_read_only);
                data_flow.distribution_x->unlock(lm_read_only);
                data_flow.distribution_y->unlock(lm_read_only);

                data_flow.f_temp_1->unlock(lm_write_only);
                data_flow.f_temp_2->unlock(lm_write_only);
                data_flow.f_temp_3->unlock(lm_write_only);
                data_flow.f_temp_4->unlock(lm_write_only);
                data_flow.f_temp_5->unlock(lm_write_only);
                data_flow.f_temp_6->unlock(lm_write_only);
                data_flow.f_temp_7->unlock(lm_write_only);
                data_flow.f_temp_8->unlock(lm_write_only);
            }
    };
}

#endif
