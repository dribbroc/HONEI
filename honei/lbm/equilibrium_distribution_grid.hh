/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
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
             * \brief Computes the equilibrium distribution for the zeroth direction.
             *
             * \param result The destination matrix.
             * \param h The height matrix.
             * \param g The gravitational constant to be used.
             * \param e The ratio of space and time stepping.
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

                    for(unsigned long i((*info.limits)[0]); i < (*info.limits)[info.limits->size() - 1]; ++i)
                    {
                        DT1_ u((*data.u)[i]);
                        DT1_ v((*data.v)[i]);
                        DT1_ h((*data.h)[i]);

                        (*data.f_eq_0)[i] = h -
                            ((DT1_(5.) * g * h * h) / (DT1_(6.) * e * e)) -
                            ((DT1_(2.) * h) /(DT1_(3.) * e * e) * (u * u + v * v));

                        (*data.f_eq_1)[i] = ((g * h * h) /(DT1_(6.) * e * e)) +
                            ((h / (DT1_(3.) * e * e)) * ((*data.distribution_x)[1] * u + (*data.distribution_y)[1] * v)) +
                            ((h / (DT1_(2.) * e * e)) * ((*data.distribution_x)[1] * u * (*data.distribution_x)[1] * u + DT1_(2.) * (*data.distribution_x)[1] * u * (*data.distribution_y)[1] * v + (*data.distribution_y)[1] * v * (*data.distribution_y)[1] * v)) -
                            ((h / (DT1_(6.) * e * e)) * (u * u + v * v));
                        (*data.f_eq_3)[i] = ((g * h * h) /(DT1_(6.) * e * e)) +
                            ((h / (DT1_(3.) * e * e)) * ((*data.distribution_x)[3] * u + (*data.distribution_y)[3] * v)) +
                            ((h / (DT1_(2.) * e * e)) * ((*data.distribution_x)[3] * u * (*data.distribution_x)[3] * u + DT1_(2.) * (*data.distribution_x)[3] * u * (*data.distribution_y)[3] * v + (*data.distribution_y)[3] * v * (*data.distribution_y)[3] * v)) -
                            ((h / (DT1_(6.) * e * e)) * (u * u + v * v));
                        (*data.f_eq_5)[i] = ((g * h * h) /(DT1_(6.) * e * e)) +
                            ((h / (DT1_(3.) * e * e)) * ((*data.distribution_x)[5] * u + (*data.distribution_y)[5] * v)) +
                            ((h / (DT1_(2.) * e * e)) * ((*data.distribution_x)[5] * u * (*data.distribution_x)[5] * u + DT1_(2.) * (*data.distribution_x)[5] * u * (*data.distribution_y)[5] * v + (*data.distribution_y)[5] * v * (*data.distribution_y)[5] * v)) -
                            ((h / (DT1_(6.) * e * e)) * (u * u + v * v));
                        (*data.f_eq_7)[i] = ((g * h * h) /(DT1_(6.) * e * e)) +
                            ((h / (DT1_(3.) * e * e)) * ((*data.distribution_x)[7] * u + (*data.distribution_y)[7] * v)) +
                            ((h / (DT1_(2.) * e * e)) * ((*data.distribution_x)[7] * u * (*data.distribution_x)[7] * u + DT1_(2.) * (*data.distribution_x)[7] * u * (*data.distribution_y)[7] * v + (*data.distribution_y)[7] * v * (*data.distribution_y)[7] * v)) -
                            ((h / (DT1_(6.) * e * e)) * (u * u + v * v));

                        (*data.f_eq_2)[i] = ((g * h * h) /(DT1_(24.) * e * e)) +
                            ((h / (DT1_(12.) * e * e)) * ((*data.distribution_x)[2] * u + (*data.distribution_y)[2] * v)) +
                            ((h / (DT1_(8.) * e * e)) * ((*data.distribution_x)[2] * u * (*data.distribution_x)[2] * u + DT1_(2.) * (*data.distribution_x)[2] * u * (*data.distribution_y)[2] * v + (*data.distribution_y)[2] * v * (*data.distribution_y)[2] * v)) -
                            ((h / (DT1_(24.) * e * e)) * (u * u + v * v));

                        (*data.f_eq_4)[i] = ((g * h * h) /(DT1_(24.) * e * e)) +
                            ((h / (DT1_(12.) * e * e)) * ((*data.distribution_x)[4] * u + (*data.distribution_y)[4] * v)) +
                            ((h / (DT1_(8.) * e * e)) * ((*data.distribution_x)[4] * u * (*data.distribution_x)[4] * u + DT1_(2.) * (*data.distribution_x)[4] * u * (*data.distribution_y)[4] * v + (*data.distribution_y)[4] * v * (*data.distribution_y)[4] * v)) -
                            ((h / (DT1_(24.) * e * e)) * (u * u + v * v));

                        (*data.f_eq_6)[i] = ((g * h * h) /(DT1_(24.) * e * e)) +
                            ((h / (DT1_(12.) * e * e)) * ((*data.distribution_x)[6] * u + (*data.distribution_y)[6] * v)) +
                            ((h / (DT1_(8.) * e * e)) * ((*data.distribution_x)[6] * u * (*data.distribution_x)[6] * u + DT1_(2.) * (*data.distribution_x)[6] * u * (*data.distribution_y)[6] * v + (*data.distribution_y)[6] * v * (*data.distribution_y)[6] * v)) -
                            ((h / (DT1_(24.) * e * e)) * (u * u + v * v));

                        (*data.f_eq_8)[i] = ((g * h * h) /(DT1_(24.) * e * e)) +
                            ((h / (DT1_(12.) * e * e)) * ((*data.distribution_x)[8] * u + (*data.distribution_y)[8] * v)) +
                            ((h / (DT1_(8.) * e * e)) * ((*data.distribution_x)[8] * u * (*data.distribution_x)[8] * u + DT1_(2.) * (*data.distribution_x)[8] * u * (*data.distribution_y)[8] * v + (*data.distribution_y)[8] * v * (*data.distribution_y)[8] * v)) -
                            ((h / (DT1_(24.) * e * e)) * (u * u + v * v));
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
}
#endif
