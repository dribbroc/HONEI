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


#ifndef LBM_GUARD_EQUILIBRIUM_DISTRIBUTION_HH
#define LBM_GUARD_EQUILIBRIUM_DISTRIBUTION_HH 1

/**
 * \file
 * Implementation of local equilibrium distribution functions used by LBM - (SWE) solvers.
 *
 * \ingroup grpliblbm
 **/

#include <honei/lbm/tags.hh>
#include <honei/la/dense_matrix.hh>

using namespace honei;
using namespace lbm;

namespace honei
{
    template<typename Tag_, typename App_, typename Direction_>
        struct EquilibriumDistribution
        {
        };

    /**
    * \brief Equilibrium distribution for direction 0.
    *
    * \ingroup grplbmoperations
    */
    template<typename Tag_>
        struct EquilibriumDistribution<Tag_, lbm_applications::LABSWE, lbm_lattice_types::D2Q9::DIR_0>
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
                static void value(DenseMatrix<DT1_>& result, DenseMatrix<DT1_>& h, DenseMatrix<DT1_>& u, DenseMatrix<DT1_>& v, DT2_ g, DT2_ e)
                {
                    CONTEXT("When computing LABSWE local equilibrium distribution function (direction 0):");
                    for(unsigned long i(0); i < h.rows(); ++i)
                    {
                        for(unsigned long j(0); j < h.columns(); ++j)
                        {
                            result(i,j) = h(i,j) -
                                          ((DT1_(5.) * g * h(i,j) * h(i,j)) / (DT1_(6.) * e * e)) -
                                          ((DT1_(2.) * h(i,j)) /(DT1_(3.) * e * e) * (u(i,j) * u(i,j) + v(i,j) * v(i,j)));
                        }
                    }
                }
        };
    /**
     * \brief Equilibrium distribution for odd direction.
     *
     * \ingroup grplbmoperations
     */
    template<typename Tag_>
        struct EquilibriumDistribution<Tag_, lbm_applications::LABSWE, lbm_lattice_types::D2Q9::DIR_ODD>
        {
            /**
             * \name Equilibrium distribution.
             *
             * \brief Computes the equilibrium distribution for odd directions.
             *
             * \param result The destination matrix.
             * \param h The height matrix.
             * \param g The gravitational constant to be used.
             * \param u The velocity in x direction.
             * \param v The velocity in y direction.
             * \param e The ratio of space and time stepping.
             * \param e_u The corresponding distribution vector entry.
             * \param e_v The corresponding distribution vector entry.
             */
            template<typename DT1_, typename DT2_>
                static void value(DenseMatrix<DT1_>& result, DenseMatrix<DT1_>& h, DenseMatrix<DT1_>& u, DenseMatrix<DT1_>& v, DT2_ g, DT2_ e, DT2_ e_u, DT2_ e_v)
                {
                    CONTEXT("When computing LABSWE local equilibrium distribution function (odd direction):");
                    for(unsigned long i(0); i < h.rows(); ++i)
                    {
                        for(unsigned long j(0); j < h.columns(); ++j)
                        {
                            result(i,j) = ((g * h(i,j) * h(i,j)) /(DT1_(6.) * e * e)) +
                                          ((h(i,j) / (DT1_(3.) * e * e)) * (e_u * u(i,j) + e_v * v(i,j))) +
                                          ((h(i,j) / (DT1_(2.) * e * e)) * (e_u * u(i,j) * e_u * u(i,j) + DT1_(2.) * e_u * u(i,j) * e_v * v(i,j) + e_v * v(i,j) * e_v * v(i,j))) -
                                          ((h(i,j) / (DT1_(6.) * e * e)) * (u(i,j) * u(i,j) + v(i,j) * v(i,j)));
                        }
                    }
                }
        };
    /**
     * \brief Equilibrium distribution for even direction.
     * \ingroup grplbmoperations
     */
    template<typename Tag_>
        struct EquilibriumDistribution<Tag_, lbm_applications::LABSWE, lbm_lattice_types::D2Q9::DIR_EVEN>
        {
            /**
             * \name Equilibrium distribution.
             *
             * \brief Computes the equilibrium distribution for even directions.
             *
             * \param result The destination matrix.
             * \param h The height matrix.
             * \param g The gravitational constant to be used.
             * \param u The velocity in x direction.
             * \param v The velocity in y direction.
             * \param e The ratio of space and time stepping.
             * \param e_u The corresponding distribution vector entry.
             * \param e_v The corresponding distribution vector entry.
             */
            template<typename DT1_, typename DT2_>
                static void value(DenseMatrix<DT1_>& result, DenseMatrix<DT1_>& h, DenseMatrix<DT1_>& u, DenseMatrix<DT1_>& v, DT2_ g, DT2_ e, DT2_ e_u, DT2_ e_v)
                {
                    CONTEXT("When computing LABSWE local equilibrium distribution function (even direction):");
                    for(unsigned long i(0); i < h.rows(); ++i)
                    {
                        for(unsigned long j(0); j < h.columns(); ++j)
                        {
                            result(i,j) = ((g * h(i,j) * h(i,j)) /(DT1_(24.) * e * e)) +
                                          ((h(i,j) / (DT1_(12.) * e * e)) * (e_u * u(i,j) + e_v * v(i,j))) +
                                          ((h(i,j) / (DT1_(8.) * e * e)) * (e_u * u(i,j) * e_u * u(i,j) + DT1_(2.) * e_u * u(i,j) * e_v * v(i,j) + e_v * v(i,j) * e_v * v(i,j))) -
                                          ((h(i,j) / (DT1_(24.) * e * e)) * (u(i,j) * u(i,j) + v(i,j) * v(i,j)));

                        }
                    }
                }
        };

//---------------------------------------------------NAVIER STOKES-------------------------------------------------------------------------------------
    /**
    * \brief Equilibrium distribution for direction 0.
    *
    * \ingroup grplbmoperations
    */
    template<typename Tag_>
        struct EquilibriumDistribution<Tag_, lbm_applications::LABNAVSTO, lbm_lattice_types::D2Q9::DIR_0>
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
                static void value(DenseMatrix<DT1_>& result, DenseMatrix<DT1_>& h, DenseMatrix<DT1_>& u, DenseMatrix<DT1_>& v, DT2_ g, DT2_ e)
                {
                    CONTEXT("When computing LABNAVSTO local equilibrium distribution function (direction 0):");
                    for(unsigned long i(0); i < h.rows(); ++i)
                    {
                        for(unsigned long j(0); j < h.columns(); ++j)
                        {
                            result(i,j) = DT1_(4./9.) * h(i,j) -
                                          ((DT1_(2.) * h(i,j)) /(DT1_(3.) * e * e) * (u(i,j) * u(i,j) + v(i,j) * v(i,j)));
                        }
                    }
                }
        };
    /**
     * \brief Equilibrium distribution for odd direction.
     *
     * \ingroup grplbmoperations
     */
    template<typename Tag_>
        struct EquilibriumDistribution<Tag_, lbm_applications::LABNAVSTO, lbm_lattice_types::D2Q9::DIR_ODD>
        {
            /**
             * \name Equilibrium distribution.
             *
             * \brief Computes the equilibrium distribution for odd directions.
             *
             * \param result The destination matrix.
             * \param h The height matrix.
             * \param g The gravitational constant to be used.
             * \param u The velocity in x direction.
             * \param v The velocity in y direction.
             * \param e The ratio of space and time stepping.
             * \param e_u The corresponding distribution vector entry.
             * \param e_v The corresponding distribution vector entry.
             */
            template<typename DT1_, typename DT2_>
                static void value(DenseMatrix<DT1_>& result, DenseMatrix<DT1_>& h, DenseMatrix<DT1_>& u, DenseMatrix<DT1_>& v, DT2_ g, DT2_ e, DT2_ e_u, DT2_ e_v)
                {
                    CONTEXT("When computing LABNAVSTO local equilibrium distribution function (odd direction):");
                    for(unsigned long i(0); i < h.rows(); ++i)
                    {
                        for(unsigned long j(0); j < h.columns(); ++j)
                        {
                            result(i,j) = ((g * h(i,j) * h(i,j)) /(DT1_(6.) * e * e)) +
                                          ((h(i,j) / (DT1_(3.) * e * e)) * (e_u * u(i,j) + e_v * v(i,j))) +
                                          ((h(i,j) / (DT1_(2.) * e * e)) * (e_u * u(i,j) * e_u * u(i,j) + DT1_(2.) * e_u * u(i,j) * e_v * v(i,j) + e_v * v(i,j) * e_v * v(i,j))) -
                                          ((h(i,j) / (DT1_(6.) * e * e)) * (u(i,j) * u(i,j) + v(i,j) * v(i,j))) +
                                          DT1_(1./9.) * h(i,j);
                        }
                    }
                }
        };
    /**
     * \brief Equilibrium distribution for even direction.
     * \ingroup grplbmoperations
     */
    template<typename Tag_>
        struct EquilibriumDistribution<Tag_, lbm_applications::LABNAVSTO, lbm_lattice_types::D2Q9::DIR_EVEN>
        {
            /**
             * \name Equilibrium distribution.
             *
             * \brief Computes the equilibrium distribution for even directions.
             *
             * \param result The destination matrix.
             * \param h The height matrix.
             * \param g The gravitational constant to be used.
             * \param u The velocity in x direction.
             * \param v The velocity in y direction.
             * \param e The ratio of space and time stepping.
             * \param e_u The corresponding distribution vector entry.
             * \param e_v The corresponding distribution vector entry.
             */
            template<typename DT1_, typename DT2_>
                static void value(DenseMatrix<DT1_>& result, DenseMatrix<DT1_>& h, DenseMatrix<DT1_>& u, DenseMatrix<DT1_>& v, DT2_ g, DT2_ e, DT2_ e_u, DT2_ e_v)
                {
                    CONTEXT("When computing LABNAVSTO local equilibrium distribution function (even direction):");
                    for(unsigned long i(0); i < h.rows(); ++i)
                    {
                        for(unsigned long j(0); j < h.columns(); ++j)
                        {
                            result(i,j) = ((g * h(i,j) * h(i,j)) /(DT1_(24.) * e * e)) +
                                          ((h(i,j) / (DT1_(12.) * e * e)) * (e_u * u(i,j) + e_v * v(i,j))) +
                                          ((h(i,j) / (DT1_(8.) * e * e)) * (e_u * u(i,j) * e_u * u(i,j) + DT1_(2.) * e_u * u(i,j) * e_v * v(i,j) + e_v * v(i,j) * e_v * v(i,j))) -
                                          ((h(i,j) / (DT1_(24.) * e * e)) * (u(i,j) * u(i,j) + v(i,j) * v(i,j))) +
                                          DT1_(1./36.) * h(i,j);
                        }
                    }
                }
        };


}
#endif
