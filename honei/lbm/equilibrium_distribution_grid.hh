/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008, 2009 Markus Geveler <apryde@gmx.de>
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
#ifndef LBM_GUARD_EQUILIBRIUM_DISTRIBUTION_GRID_HH
#define LBM_GUARD_EQUILIBRIUM_DISTRIBUTION_GRID_HH 1

/**
 * \file
 * Implementation of local equilibrium distribution functions used by LBM - (SWE) solvers using a PackedGrid.
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
        struct EquilibriumDistributionGrid
        {
        };
    /**
    * \brief Equilibrium distribution for direction 0.
    *
    * \ingroup grplbmoperations
    */
    template<>
        struct EquilibriumDistributionGrid<tags::CPU, lbm_applications::LABSWE>
        {
            /**
             * \name Equilibrium distribution.
             *
             */
            template<typename DT1_, typename DT2_>
                static void value(DT2_ g, DT2_ e, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT1_> & data)
                {
                    CONTEXT("When computing LABSWE local equilibrium distribution function:");

                    info.limits->lock(lm_read_only);

                    data.u->lock(lm_read_only);
                    data.v->lock(lm_read_only);
                    data.h->lock(lm_read_only);

                    data.distribution_x->lock(lm_read_only);
                    data.distribution_y->lock(lm_read_only);

                    data.f_eq_0->lock(lm_write_only);
                    data.f_eq_1->lock(lm_write_only);
                    data.f_eq_2->lock(lm_write_only);
                    data.f_eq_3->lock(lm_write_only);
                    data.f_eq_4->lock(lm_write_only);
                    data.f_eq_5->lock(lm_write_only);
                    data.f_eq_6->lock(lm_write_only);
                    data.f_eq_7->lock(lm_write_only);
                    data.f_eq_8->lock(lm_write_only);

                    DT1_ e2(e);
                    DT1_ e42(DT1_(2.) * e2 * e2);
                    DT1_ e23(DT1_(3.) * e2);
                    DT1_ e26(DT1_(6.) * e2);
                    DT1_ e48(DT1_(8.) * e2 * e2);
                    DT1_ e212(DT1_(12.) * e2);
                    DT1_ e224(DT1_(24.) * e2);
                    for(unsigned long i((*info.limits)[0]); i < (*info.limits)[info.limits->size() - 1]; ++i)
                    {
                        DT1_ u((*data.u)[i]);
                        DT1_ v((*data.v)[i]);
                        DT1_ h((*data.h)[i]);
                        DT1_ u2(u * u);
                        DT1_ v2(v * v);
                        DT1_ gh(g * h);

                        DT1_ dxu, dyv;
                        DT1_ t1, t2, t3, t4;

                        t1 = (DT1_(5.) * gh) / e26;
                        t2 = DT1_(2.) / e23 * (u2 + v2);
                        (*data.f_eq_0)[i] = h * (DT1_(1) - t1 - t2);

                        dxu = (*data.distribution_x)[1] * u;
                        dyv = (*data.distribution_y)[1] * v;
                        t1 = (gh) / e26;
                        t2 = (dxu + dyv) / e23;
                        t3 = (dxu * dxu + DT1_(2.) * dxu * dyv + dyv * dyv) / e42;
                        t4 = (u2 + v2) / e26;
                        (*data.f_eq_1)[i] = h * (t1 + t2 + t3 - t4);

                        dxu = (*data.distribution_x)[3] * u;
                        dyv = (*data.distribution_y)[3] * v;
                        t1 = (gh) / e26;
                        t2 = (dxu + dyv) / e23;
                        t3 = (dxu * dxu + DT1_(2.) * dxu * dyv + dyv * dyv) / e42;
                        t4 = (u2 + v2) / e26;
                        (*data.f_eq_3)[i] = h * (t1 + t2 + t3 - t4);

                        dxu = (*data.distribution_x)[5] * u;
                        dyv = (*data.distribution_y)[5] * v;
                        t1 = (gh) / e26;
                        t2 = (dxu + dyv) / e23;
                        t3 = (dxu * dxu + DT1_(2.) * dxu * dyv + dyv * dyv) / e42;
                        t4 = (u2 + v2) / e26;
                        (*data.f_eq_5)[i] = h * (t1 + t2 + t3 - t4);

                        dxu = (*data.distribution_x)[7] * u;
                        dyv = (*data.distribution_y)[7] * v;
                        t1 = (gh) / e26;
                        t2 = (dxu + dyv) / e23;
                        t3 = (dxu * dxu + DT1_(2.) * dxu * dyv + dyv * dyv) / e42;
                        t4 = (u2 + v2) / e26;
                        (*data.f_eq_7)[i] = h * (t1 + t2 + t3 - t4);

                        dxu = (*data.distribution_x)[2] * u;
                        dyv = (*data.distribution_y)[2] * v;
                        t1 = (gh) / e224;
                        t2 = (dxu + dyv) / e212;
                        t3 = (dxu * dxu + DT1_(2.) * dxu * dyv + dyv * dyv) / e48;
                        t4 = (u2 + v2) / e224;
                        (*data.f_eq_2)[i] =  h * (t1 + t2 + t3 - t4);

                        dxu = (*data.distribution_x)[4] * u;
                        dyv = (*data.distribution_y)[4] * v;
                        t1 = (gh) / e224;
                        t2 = (dxu + dyv) / e212;
                        t3 = (dxu * dxu + DT1_(2.) * dxu * dyv + dyv * dyv) / e48;
                        t4 = (u2 + v2) / e224;
                        (*data.f_eq_4)[i] =  h * (t1 + t2 + t3 - t4);

                        dxu = (*data.distribution_x)[6] * u;
                        dyv = (*data.distribution_y)[6] * v;
                        t1 = (gh) / e224;
                        t2 = (dxu + dyv) / e212;
                        t3 = (dxu * dxu + DT1_(2.) * dxu * dyv + dyv * dyv) / e48;
                        t4 = (u2 + v2) / e224;
                        (*data.f_eq_6)[i] =  h * (t1 + t2 + t3 - t4);

                        dxu = (*data.distribution_x)[8] * u;
                        dyv = (*data.distribution_y)[8] * v;
                        t1 = (gh) / e224;
                        t2 = (dxu + dyv) / e212;
                        t3 = (dxu * dxu + DT1_(2.) * dxu * dyv + dyv * dyv) / e48;
                        t4 = (u2 + v2) / e224;
                        (*data.f_eq_8)[i] =  h * (t1 + t2 + t3 - t4);
                    }

                    info.limits->unlock(lm_read_only);

                    data.u->unlock(lm_read_only);
                    data.v->unlock(lm_read_only);
                    data.h->unlock(lm_read_only);

                    data.distribution_x->unlock(lm_read_only);
                    data.distribution_y->unlock(lm_read_only);

                    data.f_eq_0->unlock(lm_write_only);
                    data.f_eq_1->unlock(lm_write_only);
                    data.f_eq_2->unlock(lm_write_only);
                    data.f_eq_3->unlock(lm_write_only);
                    data.f_eq_4->unlock(lm_write_only);
                    data.f_eq_5->unlock(lm_write_only);
                    data.f_eq_6->unlock(lm_write_only);
                    data.f_eq_7->unlock(lm_write_only);
                    data.f_eq_8->unlock(lm_write_only);
                }

            template<typename DT1_>
            static inline BenchmarkInfo get_benchmark_info(HONEI_UNUSED PackedGridInfo<D2Q9> * info, PackedGridData<D2Q9, DT1_> * data)
                {
                    BenchmarkInfo result;
                    result.flops = data->h->size() * 12 + 8 * data->h->size() * 27;
                    result.load = data->h->size() * 9 * 5 * sizeof(DT1_);
                    result.store = data->h->size() * 9 * sizeof(DT1_);
                    result.size.push_back(data->h->size());
                    return result;
                }
        };

    template<typename Tag_>
        struct EquilibriumDistributionGrid<Tag_, lbm_applications::LABNAVSTO>
        {
            /**
             * \name Equilibrium distribution for NAVSTO equations
             *
             */
            template<typename DT1_, typename DT2_>
                static void value(HONEI_UNUSED DT2_ g, DT2_ e, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT1_> & data)
                {
                    CONTEXT("When computing LABNAVSTO local equilibrium distribution function:");

                    info.limits->lock(lm_read_only);

                    data.u->lock(lm_read_only);
                    data.v->lock(lm_read_only);
                    data.h->lock(lm_read_only);

                    data.distribution_x->lock(lm_read_only);
                    data.distribution_y->lock(lm_read_only);

                    data.f_eq_0->lock(lm_write_only);
                    data.f_eq_1->lock(lm_write_only);
                    data.f_eq_2->lock(lm_write_only);
                    data.f_eq_3->lock(lm_write_only);
                    data.f_eq_4->lock(lm_write_only);
                    data.f_eq_5->lock(lm_write_only);
                    data.f_eq_6->lock(lm_write_only);
                    data.f_eq_7->lock(lm_write_only);
                    data.f_eq_8->lock(lm_write_only);

                    DT1_ e2(e); //e squared is passed!!
                    DT1_ e4(e2 * e2);
                    DT1_ three_by_e2(DT1_(3.) / e2);
                    DT1_ nine_by_2e4(DT1_(9.) / (DT1_(2.) * e4));
                    DT1_ three_by_2e2(DT1_(3.) / (DT1_(2.) * e2));
                    for(unsigned long i((*info.limits)[0]); i < (*info.limits)[info.limits->size() - 1]; ++i)
                    {
                        DT1_ u((*data.u)[i]);
                        DT1_ v((*data.v)[i]);
                        DT1_ h((*data.h)[i]);
                        DT1_ u2(u * u);
                        DT1_ v2(v * v);

                        DT1_ dxu, dyv;
                        DT1_ omega_0(DT1_(4./9.));
                        DT1_ omega_1(DT1_(1./9.));
                        DT1_ omega_2(DT1_(1./36.));
                        DT1_ one(1.);

                        (*data.f_eq_0)[i] = h * omega_0 * (one - (three_by_2e2 * (u2 + v2)));

                        dxu = (*data.distribution_x)[1] * u;
                        dyv = (*data.distribution_y)[1] * v;
                        (*data.f_eq_1)[i] = h * omega_1 * (one + (three_by_e2 * (dxu + dyv)) + (nine_by_2e4 * ((dxu + dyv) * (dxu + dyv))) -  (three_by_2e2 * (u2 + v2)));

                        dxu = (*data.distribution_x)[3] * u;
                        dyv = (*data.distribution_y)[3] * v;
                        (*data.f_eq_3)[i] = h * omega_1 * (one + (three_by_e2 * (dxu + dyv)) + (nine_by_2e4 * ((dxu + dyv) * (dxu + dyv))) -  (three_by_2e2 * (u2 + v2)));

                        dxu = (*data.distribution_x)[5] * u;
                        dyv = (*data.distribution_y)[5] * v;
                        (*data.f_eq_5)[i] = h * omega_1 * (one + (three_by_e2 * (dxu + dyv)) + (nine_by_2e4 * ((dxu + dyv) * (dxu + dyv))) -  (three_by_2e2 * (u2 + v2)));

                        dxu = (*data.distribution_x)[7] * u;
                        dyv = (*data.distribution_y)[7] * v;
                        (*data.f_eq_7)[i] = h * omega_1 * (one + (three_by_e2 * (dxu + dyv)) + (nine_by_2e4 * ((dxu + dyv) * (dxu + dyv))) -  (three_by_2e2 * (u2 + v2)));

                        dxu = (*data.distribution_x)[2] * u;
                        dyv = (*data.distribution_y)[2] * v;
                        (*data.f_eq_2)[i] = h * omega_2 * (one + (three_by_e2 * (dxu + dyv)) + (nine_by_2e4 * ((dxu + dyv) * (dxu + dyv))) -  (three_by_2e2 * (u2 + v2)));

                        dxu = (*data.distribution_x)[4] * u;
                        dyv = (*data.distribution_y)[4] * v;
                        (*data.f_eq_4)[i] = h * omega_2 * (one + (three_by_e2 * (dxu + dyv)) + (nine_by_2e4 * ((dxu + dyv) * (dxu + dyv))) -  (three_by_2e2 * (u2 + v2)));

                        dxu = (*data.distribution_x)[6] * u;
                        dyv = (*data.distribution_y)[6] * v;
                        (*data.f_eq_6)[i] = h * omega_2 * (one + (three_by_e2 * (dxu + dyv)) + (nine_by_2e4 * ((dxu + dyv) * (dxu + dyv))) -  (three_by_2e2 * (u2 + v2)));

                        dxu = (*data.distribution_x)[8] * u;
                        dyv = (*data.distribution_y)[8] * v;
                        (*data.f_eq_8)[i] = h * omega_2 * (one + (three_by_e2 * (dxu + dyv)) + (nine_by_2e4 * ((dxu + dyv) * (dxu + dyv))) -  (three_by_2e2 * (u2 + v2)));
                    }

                    info.limits->unlock(lm_read_only);

                    data.u->unlock(lm_read_only);
                    data.v->unlock(lm_read_only);
                    data.h->unlock(lm_read_only);

                    data.distribution_x->unlock(lm_read_only);
                    data.distribution_y->unlock(lm_read_only);

                    data.f_eq_0->unlock(lm_write_only);
                    data.f_eq_1->unlock(lm_write_only);
                    data.f_eq_2->unlock(lm_write_only);
                    data.f_eq_3->unlock(lm_write_only);
                    data.f_eq_4->unlock(lm_write_only);
                    data.f_eq_5->unlock(lm_write_only);
                    data.f_eq_6->unlock(lm_write_only);
                    data.f_eq_7->unlock(lm_write_only);
                    data.f_eq_8->unlock(lm_write_only);
                }

        };

    template<>
        struct EquilibriumDistributionGrid<tags::GPU::CUDA, lbm_applications::LABSWE>
        {
            static void value(float g, float e, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data);
            static void value(double g, double e, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, double> & data);
        };

    template<>
        struct EquilibriumDistributionGrid<tags::Cell, lbm_applications::LABSWE>
        {
            static void value(float g, float e, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data);
        };

    template<>
        struct EquilibriumDistributionGrid<tags::CPU::SSE, lbm_applications::LABSWE>
        {
            template <typename DT1_>
                static void value(DT1_ g, DT1_ e, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT1_> & data);
        };

    template<>
        struct EquilibriumDistributionGrid<tags::CPU::Itanium, lbm_applications::LABSWE>
        {
            template <typename DT1_>
            static void value(DT1_ g, DT1_ e, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT1_> & data);
        };
}
#endif
