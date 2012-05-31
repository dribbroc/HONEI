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
#ifndef LBM_GUARD_EQUILIBRIUM_DISTRIBUTION_GRID_POLLUTANT_HH
#define LBM_GUARD_EQUILIBRIUM_DISTRIBUTION_GRID_POLLUTANT_HH 1

/**
 * \file
 * Implementation of local equilibrium distribution functions for pollutant transport (convection diffusion) used by LBM - (SWE) solvers using a PackedGrid.
 *
 * \ingroup grpliblbm
 **/

#include <honei/lbm/tags.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/lbm/grid.hh>
#include <honei/util/benchmark_info.hh>

using namespace honei;
using namespace lbm;

namespace honei
{
    template<typename Tag_, typename App_>
        struct EquilibriumDistributionGridPollutant
        {
        };
    /**
    * \brief Equilibrium distribution for direction 0.
    *
    * \ingroup grplbmoperations
    */
    template<>
        struct EquilibriumDistributionGridPollutant<tags::CPU::Generic, lbm_applications::LABSWE>
        {
            /**
             * \name Equilibrium distribution.
             *
             */
            template<typename DT1_, typename DT2_>
                static void value(DT2_ e, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT1_> & data_flow, PackedGridData<D2Q9, DT1_> & data_poll)
                {
                    CONTEXT("When computing local equilibrium distribution function for pollutant transport (convection diffusion equation):");

                    info.limits->lock(lm_read_only);

                    data_flow.u->lock(lm_read_only);
                    data_flow.v->lock(lm_read_only);
                    data_flow.h->lock(lm_read_only);
                    data_poll.h->lock(lm_read_only);

                    data_flow.distribution_x->lock(lm_read_only);
                    data_flow.distribution_y->lock(lm_read_only);

                    data_poll.f_eq_0->lock(lm_write_only);
                    data_poll.f_eq_1->lock(lm_write_only);
                    data_poll.f_eq_2->lock(lm_write_only);
                    data_poll.f_eq_3->lock(lm_write_only);
                    data_poll.f_eq_4->lock(lm_write_only);
                    data_poll.f_eq_5->lock(lm_write_only);
                    data_poll.f_eq_6->lock(lm_write_only);
                    data_poll.f_eq_7->lock(lm_write_only);
                    data_poll.f_eq_8->lock(lm_write_only);

                    const unsigned long * const limits(info.limits->elements());
                    const DT1_ * const tu(data_flow.u->elements());
                    const DT1_ * const tv(data_flow.v->elements());
                    const DT1_ * const th(data_flow.h->elements());
                    const DT1_ * const tc(data_poll.h->elements());
                    const DT1_ * const distribution_x(data_flow.distribution_x->elements());
                    const DT1_ * const distribution_y(data_flow.distribution_y->elements());

                    DT1_ * f_eq_0(data_poll.f_eq_0->elements());
                    DT1_ * f_eq_1(data_poll.f_eq_1->elements());
                    DT1_ * f_eq_2(data_poll.f_eq_2->elements());
                    DT1_ * f_eq_3(data_poll.f_eq_3->elements());
                    DT1_ * f_eq_4(data_poll.f_eq_4->elements());
                    DT1_ * f_eq_5(data_poll.f_eq_5->elements());
                    DT1_ * f_eq_6(data_poll.f_eq_6->elements());
                    DT1_ * f_eq_7(data_poll.f_eq_7->elements());
                    DT1_ * f_eq_8(data_poll.f_eq_8->elements());

                    DT1_ e2(e);
                    DT1_ e2_by_3(DT1_(3.) * e2);
                    DT1_ e2_by_12(DT1_(12.) * e2);
                    DT1_ one_over_9(1./9.);
                    DT1_ five_over_9(5./9.);
                    DT1_ one_over_36(1./36.);

                    const unsigned long start(limits[0]);
                    const unsigned long end(limits[info.limits->size() - 1]);

                    DT1_ t1, t2;


                    for(unsigned long i(start); i < end; ++i)
                    {
                        const DT1_ h((th)[i]);
                        const DT1_ c((tc)[i]);

                        (f_eq_0)[i] = c * DT1_(h - five_over_9);
                    }

                    for(unsigned long i(start); i < end; ++i)
                    {
                        const DT1_ u((tu)[i]);
                        const DT1_ v((tv)[i]);
                        const DT1_ h((th)[i]);
                        const DT1_ c((tc)[i]);

                        const DT1_ dxu = (distribution_x)[1] * u;
                        const DT1_ dyv = (distribution_y)[1] * v;

                        t1 = (h / e2_by_3) * dxu;
                        t2 = (h / e2_by_3) * dyv;
                        (f_eq_1)[i] = c * (one_over_9 + t1 + t2);
                    }

                    for(unsigned long i(start); i < end; ++i)
                    {
                        const DT1_ u((tu)[i]);
                        const DT1_ v((tv)[i]);
                        const DT1_ h((th)[i]);
                        const DT1_ c((tc)[i]);

                        const DT1_ dxu = (distribution_x)[3] * u;
                        const DT1_ dyv = (distribution_y)[3] * v;

                        t1 = (h / e2_by_3) * dxu;
                        t2 = (h / e2_by_3) * dyv;
                        (f_eq_3)[i] = c * (one_over_9 + t1 + t2);
                    }

                    for(unsigned long i(start); i < end; ++i)
                    {
                        const DT1_ u((tu)[i]);
                        const DT1_ v((tv)[i]);
                        const DT1_ h((th)[i]);
                        const DT1_ c((tc)[i]);

                        const DT1_ dxu = (distribution_x)[5] * u;
                        const DT1_ dyv = (distribution_y)[5] * v;

                        t1 = (h / e2_by_3) * dxu;
                        t2 = (h / e2_by_3) * dyv;
                        (f_eq_5)[i] = c * (one_over_9 + t1 + t2);
                    }

                    for(unsigned long i(start); i < end; ++i)
                    {
                        const DT1_ u((tu)[i]);
                        const DT1_ v((tv)[i]);
                        const DT1_ h((th)[i]);
                        const DT1_ c((tc)[i]);

                        const DT1_ dxu = (distribution_x)[7] * u;
                        const DT1_ dyv = (distribution_y)[7] * v;

                        t1 = (h / e2_by_3) * dxu;
                        t2 = (h / e2_by_3) * dyv;
                        (f_eq_7)[i] = c * (one_over_9 + t1 + t2);
                    }

                    for(unsigned long i(start); i < end; ++i)
                    {
                        const DT1_ u((tu)[i]);
                        const DT1_ v((tv)[i]);
                        const DT1_ h((th)[i]);
                        const DT1_ c((tc)[i]);

                        const DT1_ dxu = (distribution_x)[2] * u;
                        const DT1_ dyv = (distribution_y)[2] * v;

                        t1 = (h / e2_by_12) * dxu;
                        t2 = (h / e2_by_12) * dyv;
                        (f_eq_2)[i] = c * (one_over_36 + t1 + t2);
                    }

                    for(unsigned long i(start); i < end; ++i)
                    {
                        const DT1_ u((tu)[i]);
                        const DT1_ v((tv)[i]);
                        const DT1_ h((th)[i]);
                        const DT1_ c((tc)[i]);

                        const DT1_ dxu = (distribution_x)[4] * u;
                        const DT1_ dyv = (distribution_y)[4] * v;

                        t1 = (h / e2_by_12) * dxu;
                        t2 = (h / e2_by_12) * dyv;
                        (f_eq_4)[i] = c * (one_over_36 + t1 + t2);
                    }

                    for(unsigned long i(start); i < end; ++i)
                    {
                        const DT1_ u((tu)[i]);
                        const DT1_ v((tv)[i]);
                        const DT1_ h((th)[i]);
                        const DT1_ c((tc)[i]);

                        const DT1_ dxu = (distribution_x)[6] * u;
                        const DT1_ dyv = (distribution_y)[6] * v;

                        t1 = (h / e2_by_12) * dxu;
                        t2 = (h / e2_by_12) * dyv;
                        (f_eq_6)[i] = c * (one_over_36 + t1 + t2);
                    }

                    for(unsigned long i(start); i < end; ++i)
                    {
                        const DT1_ u((tu)[i]);
                        const DT1_ v((tv)[i]);
                        const DT1_ h((th)[i]);
                        const DT1_ c((tc)[i]);

                        const DT1_ dxu = (distribution_x)[8] * u;
                        const DT1_ dyv = (distribution_y)[8] * v;

                        t1 = (h / e2_by_12) * dxu;
                        t2 = (h / e2_by_12) * dyv;
                        (f_eq_8)[i] = c * (one_over_36 + t1 + t2);
                    }

                    info.limits->unlock(lm_read_only);

                    data_flow.u->unlock(lm_read_only);
                    data_flow.v->unlock(lm_read_only);
                    data_flow.h->unlock(lm_read_only);
                    data_poll.h->unlock(lm_read_only);

                    data_flow.distribution_x->unlock(lm_read_only);
                    data_flow.distribution_y->unlock(lm_read_only);

                    data_poll.f_eq_0->unlock(lm_write_only);
                    data_poll.f_eq_1->unlock(lm_write_only);
                    data_poll.f_eq_2->unlock(lm_write_only);
                    data_poll.f_eq_3->unlock(lm_write_only);
                    data_poll.f_eq_4->unlock(lm_write_only);
                    data_poll.f_eq_5->unlock(lm_write_only);
                    data_poll.f_eq_6->unlock(lm_write_only);
                    data_poll.f_eq_7->unlock(lm_write_only);
                    data_poll.f_eq_8->unlock(lm_write_only);
                }
        };

    template<>
        struct EquilibriumDistributionGridPollutant<tags::GPU::CUDA, lbm_applications::LABSWE>
        {
                static void value(float e, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data_flow, PackedGridData<D2Q9, float> & data_poll);
                static void value(double e, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, double> & data_flow, PackedGridData<D2Q9, double> & data_poll);
        };
}


#endif
