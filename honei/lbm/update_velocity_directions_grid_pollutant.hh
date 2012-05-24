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


#pragma once
#ifndef LBM_GUARD_UPDATE_VELOCITY_DIRECTIONS_GRID_POLLUTANT_HH
#define LBM_GUARD_UPDATE_VELOCITY_DIRECTIONS_GRID_POLLUTANT_HH 1

#include <honei/lbm/tags.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/lbm/grid.hh>
#include <cmath>
using namespace honei::lbm;

namespace honei
{
   template <typename Tag_,
              typename BoundaryType_>
    struct UpdateVelocityDirectionsGridPollutant
    {
    };

   template <>
       struct UpdateVelocityDirectionsGridPollutant<tags::CPU::Generic, lbm_boundary_types::DIRICHLET_SLIP>
       {
           template<typename DT1_>
               static void value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT1_> & data_flow, PackedGridData<D2Q9, DT1_> & data_poll, DT1_ e, DT1_ value)
               {
                   CONTEXT("When updating velocity directions (Dirichlet pollutant conditions):");

                   info.limits->lock(lm_read_only);
                   info.types->lock(lm_read_only);

                   data_poll.f_temp_1->lock(lm_read_and_write);
                   data_poll.f_temp_2->lock(lm_read_and_write);
                   data_poll.f_temp_3->lock(lm_read_and_write);
                   data_poll.f_temp_4->lock(lm_read_and_write);
                   data_poll.f_temp_5->lock(lm_read_and_write);
                   data_poll.f_temp_6->lock(lm_read_and_write);
                   data_poll.f_temp_7->lock(lm_read_and_write);
                   data_poll.f_temp_8->lock(lm_read_and_write);

                   data_flow.u->lock(lm_read_only);
                   data_flow.v->lock(lm_read_only);
                   data_flow.h->lock(lm_read_only);

                   data_flow.distribution_x->lock(lm_read_only);
                   data_flow.distribution_y->lock(lm_read_only);
                   const unsigned long * const limits(info.limits->elements());
                   const unsigned long * const types(info.types->elements());

                   DT1_ * f_temp_1(data_poll.f_temp_1->elements());
                   DT1_ * f_temp_2(data_poll.f_temp_2->elements());
                   DT1_ * f_temp_3(data_poll.f_temp_3->elements());
                   DT1_ * f_temp_4(data_poll.f_temp_4->elements());
                   DT1_ * f_temp_5(data_poll.f_temp_5->elements());
                   DT1_ * f_temp_6(data_poll.f_temp_6->elements());
                   DT1_ * f_temp_7(data_poll.f_temp_7->elements());
                   DT1_ * f_temp_8(data_poll.f_temp_8->elements());

                   DT1_ * tu(data_flow.u->elements());
                   DT1_ * tv(data_flow.v->elements());
                   DT1_ * th(data_flow.h->elements());

                   DT1_ * distribution_x(data_flow.distribution_x->elements());
                   DT1_ * distribution_y(data_flow.distribution_y->elements());

                   DT1_ e2(e);
                   DT1_ e2_by_3(DT1_(3.) * e2);
                   DT1_ e2_by_12(DT1_(12.) * e2);
                   DT1_ one_over_9(1./9.);
                   DT1_ one_over_36(1./36.);

                   DT1_ t1, t2;

                   const unsigned long end(info.limits->size() - 1);
                   for (unsigned long begin(0) ; begin < end ; ++begin)
                   {
                       if(((types)[begin] & 1<<0) == 1<<0)
                           for (unsigned long i((limits)[begin]) ; i != (limits)[begin + 1] ; ++i)
                           {
                               const DT1_ u((tu)[i]);
                               const DT1_ v((tv)[i]);
                               const DT1_ h((th)[i]);

                               const DT1_ dxu = (distribution_x)[5] * u;
                               const DT1_ dyv = (distribution_y)[5] * v;

                               t1 = (h / e2_by_3) * dxu;
                               t2 = (h / e2_by_3) * dyv;
                               (f_temp_5)[i] = value * (one_over_9 + t1 + t2);
                           }
                       if(((types)[begin] & 1<<1) == 1<<1)
                           for (unsigned long i((limits)[begin]) ; i != (limits)[begin + 1] ; ++i)
                           {
                               const DT1_ u((tu)[i]);
                               const DT1_ v((tv)[i]);
                               const DT1_ h((th)[i]);

                               const DT1_ dxu = (distribution_x)[6] * u;
                               const DT1_ dyv = (distribution_y)[6] * v;

                               t1 = (h / e2_by_12) * dxu;
                               t2 = (h / e2_by_12) * dyv;
                               (f_temp_6)[i] = value * (one_over_36 + t1 + t2);
                           }
                       if(((types)[begin] & 1<<2) == 1<<2)
                           for (unsigned long i((limits)[begin]) ; i != (limits)[begin + 1] ; ++i)
                           {
                               const DT1_ u((tu)[i]);
                               const DT1_ v((tv)[i]);
                               const DT1_ h((th)[i]);

                               const DT1_ dxu = (distribution_x)[7] * u;
                               const DT1_ dyv = (distribution_y)[7] * v;

                               t1 = (h / e2_by_3) * dxu;
                               t2 = (h / e2_by_3) * dyv;
                               (f_temp_7)[i] =  value * (one_over_9 + t1 + t2);
                           }
                       if(((types)[begin] & 1<<3) == 1<<3)
                           for (unsigned long i((limits)[begin]) ; i != (limits)[begin + 1] ; ++i)
                           {
                               const DT1_ u((tu)[i]);
                               const DT1_ v((tv)[i]);
                               const DT1_ h((th)[i]);

                               const DT1_ dxu = (distribution_x)[8] * u;
                               const DT1_ dyv = (distribution_y)[8] * v;

                               t1 = (h / e2_by_12) * dxu;
                               t2 = (h / e2_by_12) * dyv;
                               (f_temp_8)[i] =  value * (one_over_36 + t1 + t2);
                           }
                       if(((types)[begin] & 1<<4) == 1<<4)
                           for (unsigned long i((limits)[begin]) ; i != (limits)[begin + 1] ; ++i)
                           {
                               const DT1_ u((tu)[i]);
                               const DT1_ v((tv)[i]);
                               const DT1_ h((th)[i]);

                               const DT1_ dxu = (distribution_x)[1] * u;
                               const DT1_ dyv = (distribution_y)[1] * v;

                               t1 = (h / e2_by_3) * dxu;
                               t2 = (h / e2_by_3) * dyv;
                               (f_temp_1)[i] = value * (one_over_9 + t1 + t2);
                           }
                       if(((types)[begin] & 1<<5) == 1<<5)
                           for (unsigned long i((limits)[begin]) ; i != (limits)[begin + 1] ; ++i)
                           {
                               const DT1_ u((tu)[i]);
                               const DT1_ v((tv)[i]);
                               const DT1_ h((th)[i]);

                               const DT1_ dxu = (distribution_x)[2] * u;
                               const DT1_ dyv = (distribution_y)[2] * v;

                               t1 = (h / e2_by_12) * dxu;
                               t2 = (h / e2_by_12) * dyv;
                               (f_temp_2)[i] = value * (one_over_36+  t1 + t2);
                           }
                       if(((types)[begin] & 1<<6) == 1<<6)
                           for (unsigned long i((limits)[begin]) ; i != (limits)[begin + 1] ; ++i)
                           {
                               const DT1_ u((tu)[i]);
                               const DT1_ v((tv)[i]);
                               const DT1_ h((th)[i]);

                               const DT1_ dxu = (distribution_x)[3] * u;
                               const DT1_ dyv = (distribution_y)[3] * v;

                               t1 = (h / e2_by_3) * dxu;
                               t2 = (h / e2_by_3) * dyv;
                               (f_temp_3)[i] = value * (one_over_9 + t1 + t2);
                           }
                       if(((types)[begin] & 1<<7) == 1<<7)
                           for (unsigned long i((limits)[begin]) ; i != (limits)[begin + 1] ; ++i)
                           {
                               const DT1_ u((tu)[i]);
                               const DT1_ v((tv)[i]);
                               const DT1_ h((th)[i]);

                               const DT1_ dxu = (distribution_x)[4] * u;
                               const DT1_ dyv = (distribution_y)[4] * v;

                               t1 = (h / e2_by_12) * dxu;
                               t2 = (h / e2_by_12) * dyv;
                               (f_temp_4)[i] = value * (one_over_36 + t1 + t2);
                           }
                   }

                   info.limits->unlock(lm_read_only);
                   info.types->unlock(lm_read_only);

                   data_poll.f_temp_1->unlock(lm_read_and_write);
                   data_poll.f_temp_2->unlock(lm_read_and_write);
                   data_poll.f_temp_3->unlock(lm_read_and_write);
                   data_poll.f_temp_4->unlock(lm_read_and_write);
                   data_poll.f_temp_5->unlock(lm_read_and_write);
                   data_poll.f_temp_6->unlock(lm_read_and_write);
                   data_poll.f_temp_7->unlock(lm_read_and_write);
                   data_poll.f_temp_8->unlock(lm_read_and_write);

                   data_flow.u->unlock(lm_read_only);
                   data_flow.v->unlock(lm_read_only);
                   data_flow.h->unlock(lm_read_only);

                   data_flow.distribution_x->unlock(lm_read_only);
                   data_flow.distribution_y->unlock(lm_read_only);
               }
       };
}

#endif
