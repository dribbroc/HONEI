/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
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
#include <grid.hh>

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
    template<typename Tag_>
        struct EquilibriumDistributionGrid<Tag_, lbm_applications::LABSWE>
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
                static void value(DT2_ g, DT2_ e, DenseVector<DT1_> & e_u, DenseVector<DT1_> & e_v, PackedGridData<D2Q9, DT1_> & grid_data)
                {
                    CONTEXT("When computing LABSWE local equilibrium distribution function:");
                    for(unsigned long i(0); i < grid_data.h->size(); ++i)
                    {
                        DT1_ u((*grid_data.u)[i]);
                        DT1_ v((*grid_data.v)[i]);
                        DT1_ h((*grid_data.h)[i]);

                        (*grid_data.f_eq_0)[i] = h -
                            ((DT1_(5.) * g * h * h) / (DT1_(6.) * e * e)) -
                            ((DT1_(2.) * h) /(DT1_(3.) * e * e) * (u * u + v * v));

                        (*grid_data.f_eq_1)[i] = ((g * h * h) /(DT1_(6.) * e * e)) +
                            ((h / (DT1_(3.) * e * e)) * (e_u[1] * u + e_v[1] * v)) +
                            ((h / (DT1_(2.) * e * e)) * (e_u[1] * u * e_u[1] * u + DT1_(2.) * e_u[1] * u * e_v[1] * v + e_v[1] * v * e_v[1] * v)) -
                            ((h / (DT1_(6.) * e * e)) * (u * u + v * v));
                        (*grid_data.f_eq_3)[i] = ((g * h * h) /(DT1_(6.) * e * e)) +
                            ((h / (DT1_(3.) * e * e)) * (e_u[3] * u + e_v[3] * v)) +
                            ((h / (DT1_(2.) * e * e)) * (e_u[3] * u * e_u[3] * u + DT1_(2.) * e_u[3] * u * e_v[3] * v + e_v[3] * v * e_v[3] * v)) -
                            ((h / (DT1_(6.) * e * e)) * (u * u + v * v));
                        (*grid_data.f_eq_5)[i] = ((g * h * h) /(DT1_(6.) * e * e)) +
                            ((h / (DT1_(3.) * e * e)) * (e_u[5] * u + e_v[5] * v)) +
                            ((h / (DT1_(2.) * e * e)) * (e_u[5] * u * e_u[5] * u + DT1_(2.) * e_u[5] * u * e_v[5] * v + e_v[5] * v * e_v[5] * v)) -
                            ((h / (DT1_(6.) * e * e)) * (u * u + v * v));
                        (*grid_data.f_eq_7)[i] = ((g * h * h) /(DT1_(6.) * e * e)) +
                            ((h / (DT1_(3.) * e * e)) * (e_u[7] * u + e_v[7] * v)) +
                            ((h / (DT1_(2.) * e * e)) * (e_u[7] * u * e_u[7] * u + DT1_(2.) * e_u[7] * u * e_v[7] * v + e_v[7] * v * e_v[7] * v)) -
                            ((h / (DT1_(6.) * e * e)) * (u * u + v * v));

                        (*grid_data.f_eq_2)[i] = ((g * h * h) /(DT1_(24.) * e * e)) +
                            ((h / (DT1_(12.) * e * e)) * (e_u[2] * u + e_v[2] * v)) +
                            ((h / (DT1_(8.) * e * e)) * (e_u[2] * u * e_u[2] * u + DT1_(2.) * e_u[2] * u * e_v[2] * v + e_v[2] * v * e_v[2] * v)) -
                            ((h / (DT1_(24.) * e * e)) * (u * u + v * v));

                        (*grid_data.f_eq_4)[i] = ((g * h * h) /(DT1_(24.) * e * e)) +
                            ((h / (DT1_(12.) * e * e)) * (e_u[4] * u + e_v[4] * v)) +
                            ((h / (DT1_(8.) * e * e)) * (e_u[4] * u * e_u[4] * u + DT1_(2.) * e_u[4] * u * e_v[4] * v + e_v[4] * v * e_v[4] * v)) -
                            ((h / (DT1_(24.) * e * e)) * (u * u + v * v));

                        (*grid_data.f_eq_6)[i] = ((g * h * h) /(DT1_(24.) * e * e)) +
                            ((h / (DT1_(12.) * e * e)) * (e_u[6] * u + e_v[6] * v)) +
                            ((h / (DT1_(8.) * e * e)) * (e_u[6] * u * e_u[6] * u + DT1_(2.) * e_u[6] * u * e_v[6] * v + e_v[6] * v * e_v[6] * v)) -
                            ((h / (DT1_(24.) * e * e)) * (u * u + v * v));

                        (*grid_data.f_eq_8)[i] = ((g * h * h) /(DT1_(24.) * e * e)) +
                            ((h / (DT1_(12.) * e * e)) * (e_u[8] * u + e_v[8] * v)) +
                            ((h / (DT1_(8.) * e * e)) * (e_u[8] * u * e_u[8] * u + DT1_(2.) * e_u[8] * u * e_v[8] * v + e_v[8] * v * e_v[8] * v)) -
                            ((h / (DT1_(24.) * e * e)) * (u * u + v * v));
                    }
                }
        };
}
#endif
