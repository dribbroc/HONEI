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


#pragma once
#ifndef LBM_GUARD_COLLIDE_STREAM_HH
#define LBM_GUARD_COLLIDE_STREAM_HH 1


/**
 * \file
 * Implementation of collision and streaming modules used by  LBM - (SWE) solvers.
 *
 * \ingroup grpliblbm
 **/

#include <honei/lbm/tags.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/dense_matrix.hh>
#include <cmath>
using namespace honei::lbm;

namespace honei
{
    template <typename Tag_, typename Application_, typename BoundaryType_, typename Direction_>
    struct CollideStream
    {
    };

    /**
     * \brief Collision and streaming module for LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::PERIODIC, lbm_lattice_types::D2Q9::DIR_1>
    {
        /**
         * \name Collision and Streaming for direction 1..
         *
         * \brief Solves the LB equation.
         *
         * \param result The destination matrix.
         * \param dist The temporary distribution matrix.
         * \param eq_dist The equilibrium distribution matrix..
         * \param s_x Source matrix in x direction.
         * \param s_y Source matrix in y direction..
         * \param e_x Corresponding distribution scalar.
         * \param e_y Corresponding distribution scalar.
         * \param tau The relaxation time.
         */
        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_>& result,
                          DenseMatrix<DT1_>& dist,
                          DenseMatrix<DT1_>& eq_dist,
                          DenseMatrix<DT1_>& s_x,
                          DenseMatrix<DT1_>& s_y,
                          DT2_ e_x,
                          DT2_ e_y,
                          DT2_ tau)
        {
            CONTEXT("When performing collision and streaming in DIR 1:");
            long y_max(result.rows());
            long x_max(result.columns());

            for(long i(0); i < y_max; ++i)
            {
                long i_forward(i + 1);
                long i_backward(i - 1);

                for(long j(0); j < x_max; ++j)
                {
                    long j_forward(j + 1);
                    long j_backward(j - 1);

                    ///Respect periodic boundaries:
                    if(j_forward >= x_max)
                        j_forward = j_forward - x_max;
                    if(j_backward < 0)
                        j_backward = j_backward + x_max;
                    if(i_forward >= y_max)
                        i_forward = i_forward - y_max;
                    if(i_backward < 0)
                        i_backward = i_backward + y_max;

                    ///Perform streaming and collision:
                    result(i,j_forward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                }
            }
        }
    };

    /**
     * \brief Collision and streaming module for LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::PERIODIC, lbm_lattice_types::D2Q9::DIR_2>
    {
        /**
         * \name Collision and Streaming for direction 2.
         *
         * \brief Solves the LB equation.
         *
         * \param result The destination matrix.
         * \param dist The temporary distribution matrix.
         * \param eq_dist The equilibrium distribution matrix..
         * \param s_x Source matrix in x direction.
         * \param s_y Source matrix in y direction..
         * \param e_x Corresponding distribution scalar.
         * \param e_y Corresponding distribution scalar.
         * \param tau The relaxation time.
         */
        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_>& result,
                          DenseMatrix<DT1_>& dist,
                          DenseMatrix<DT1_>& eq_dist,
                          DenseMatrix<DT1_>& s_x,
                          DenseMatrix<DT1_>& s_y,
                          DT2_ e_x,
                          DT2_ e_y,
                          DT2_ tau)
        {
            CONTEXT("When performing collision and streaming in DIR 2:");
            long y_max(result.rows());
            long x_max(result.columns());

            for(long i(0); i < y_max; ++i)
            {
                long i_forward(i + 1);
                long i_backward(i - 1);

                for(long j(0); j < x_max; ++j)
                {
                    long j_forward(j + 1);
                    long j_backward(j - 1);

                    ///Respect periodic boundaries:
                    if(j_forward >= x_max)
                        j_forward = j_forward - x_max;
                    if(j_backward < 0)
                        j_backward = j_backward + x_max;
                    if(i_forward >= y_max)
                        i_forward = i_forward - y_max;
                    if(i_backward < 0)
                        i_backward = i_backward + y_max;

                    ///Perform streaming and collision:
                    result(i_forward, j_forward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                }
            }
        }
    };

    /**
     * \brief Collision and streaming module for LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::PERIODIC, lbm_lattice_types::D2Q9::DIR_3>
    {
        /**
         * \name Collision and Streaming for direction 3.
         *
         * \brief Solves the LB equation.
         *
         * \param result The destination matrix.
         * \param dist The temporary distribution matrix.
         * \param eq_dist The equilibrium distribution matrix..
         * \param s_x Source matrix in x direction.
         * \param s_y Source matrix in y direction..
         * \param e_x Corresponding distribution scalar.
         * \param e_y Corresponding distribution scalar.
         * \param tau The relaxation time.
         */
        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_>& result,
                          DenseMatrix<DT1_>& dist,
                          DenseMatrix<DT1_>& eq_dist,
                          DenseMatrix<DT1_>& s_x,
                          DenseMatrix<DT1_>& s_y,
                          DT2_ e_x,
                          DT2_ e_y,
                          DT2_ tau)
        {

            CONTEXT("When performing collision and streaming in DIR 3:");
            long y_max(result.rows());
            long x_max(result.columns());

            for(long i(0); i < y_max; ++i)
            {
                long i_forward(i + 1);
                long i_backward(i - 1);

                for(long j(0); j < x_max; ++j)
                {
                    long j_forward(j + 1);
                    long j_backward(j - 1);

                    ///Respect periodic boundaries:
                    if(j_forward >= x_max)
                        j_forward = j_forward - x_max;
                    if(j_backward < 0)
                        j_backward = j_backward + x_max;
                    if(i_forward >= y_max)
                        i_forward = i_forward - y_max;
                    if(i_backward < 0)
                        i_backward = i_backward + y_max;

                    ///Perform streaming and collision:
                    result(i_forward, j) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                }
            }
        }
    };

    /**
     * \brief Collision and streaming module for LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::PERIODIC, lbm_lattice_types::D2Q9::DIR_4>
    {
        /**
         * \name Collision and Streaming for direction 4.
         *
         * \brief Solves the LB equation.
         *
         * \param result The destination matrix.
         * \param dist The temporary distribution matrix.
         * \param eq_dist The equilibrium distribution matrix..
         * \param s_x Source matrix in x direction.
         * \param s_y Source matrix in y direction..
         * \param e_x Corresponding distribution scalar.
         * \param e_y Corresponding distribution scalar.
         * \param tau The relaxation time.
         */
        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_>& result,
                          DenseMatrix<DT1_>& dist,
                          DenseMatrix<DT1_>& eq_dist,
                          DenseMatrix<DT1_>& s_x,
                          DenseMatrix<DT1_>& s_y,
                          DT2_ e_x,
                          DT2_ e_y,
                          DT2_ tau)
        {

            CONTEXT("When performing collision and streaming in DIR 4:");
            long y_max(result.rows());
            long x_max(result.columns());

            for(long i(0); i < y_max; ++i)
            {
                long i_forward(i + 1);
                long i_backward(i - 1);

                for(long j(0); j < x_max; ++j)
                {
                    long j_forward(j + 1);
                    long j_backward(j - 1);

                    ///Respect periodic boundaries:
                    if(j_forward >= x_max)
                        j_forward = j_forward - x_max;
                    if(j_backward < 0)
                        j_backward = j_backward + x_max;
                    if(i_forward >= y_max)
                        i_forward = i_forward - y_max;
                    if(i_backward < 0)
                        i_backward = i_backward + y_max;

                    ///Perform streaming and collision:
                    result(i_forward, j_backward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                }
            }
        }
    };

    /**
     * \brief Collision and streaming module for LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::PERIODIC, lbm_lattice_types::D2Q9::DIR_5>
    {
        /**
         * \name Collision and Streaming for direction 5.
         *
         * \brief Solves the LB equation.
         *
         * \param result The destination matrix.
         * \param dist The temporary distribution matrix.
         * \param eq_dist The equilibrium distribution matrix..
         * \param s_x Source matrix in x direction.
         * \param s_y Source matrix in y direction..
         * \param e_x Corresponding distribution scalar.
         * \param e_y Corresponding distribution scalar.
         * \param tau The relaxation time.
         */
        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_>& result,
                          DenseMatrix<DT1_>& dist,
                          DenseMatrix<DT1_>& eq_dist,
                          DenseMatrix<DT1_>& s_x,
                          DenseMatrix<DT1_>& s_y,
                          DT2_ e_x,
                          DT2_ e_y,
                          DT2_ tau)
        {

            CONTEXT("When performing collision and streaming in DIR 5:");
            long y_max(result.rows());
            long x_max(result.columns());

            for(long i(0); i < y_max; ++i)
            {
                long i_forward(i + 1);
                long i_backward(i - 1);

                for(long j(0); j < x_max; ++j)
                {
                    long j_forward(j + 1);
                    long j_backward(j - 1);

                    ///Respect periodic boundaries:
                    if(j_forward >= x_max)
                        j_forward = j_forward - x_max;
                    if(j_backward < 0)
                        j_backward = j_backward + x_max;
                    if(i_forward >= y_max)
                        i_forward = i_forward - y_max;
                    if(i_backward < 0)
                        i_backward = i_backward + y_max;

                    ///Perform streaming and collision:
                    result(i, j_backward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                }
            }
        }
    };

    /**
     * \brief Collision and streaming module for LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::PERIODIC, lbm_lattice_types::D2Q9::DIR_6>
    {
        /**
         * \name Collision and Streaming for direction 6.
         *
         * \brief Solves the LB equation.
         *
         * \param result The destination matrix.
         * \param dist The temporary distribution matrix.
         * \param eq_dist The equilibrium distribution matrix..
         * \param s_x Source matrix in x direction.
         * \param s_y Source matrix in y direction..
         * \param e_x Corresponding distribution scalar.
         * \param e_y Corresponding distribution scalar.
         * \param tau The relaxation time.
         */
        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_>& result,
                          DenseMatrix<DT1_>& dist,
                          DenseMatrix<DT1_>& eq_dist,
                          DenseMatrix<DT1_>& s_x,
                          DenseMatrix<DT1_>& s_y,
                          DT2_ e_x,
                          DT2_ e_y,
                          DT2_ tau)
        {

            CONTEXT("When performing collision and streaming in DIR 6:");
            long y_max(result.rows());
            long x_max(result.columns());

            for(long i(0); i < y_max; ++i)
            {
                long i_forward(i + 1);
                long i_backward(i - 1);

                for(long j(0); j < x_max; ++j)
                {
                    long j_forward(j + 1);
                    long j_backward(j - 1);

                    ///Respect periodic boundaries:
                    if(j_forward >= x_max)
                        j_forward = j_forward - x_max;
                    if(j_backward < 0)
                        j_backward = j_backward + x_max;
                    if(i_forward >= y_max)
                        i_forward = i_forward - y_max;
                    if(i_backward < 0)
                        i_backward = i_backward + y_max;

                    ///Perform streaming and collision:
                    result(i_backward, j_backward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                }
            }
        }
    };

    /**
     * \brief Collision and streaming module for LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::PERIODIC, lbm_lattice_types::D2Q9::DIR_7>
    {
        /**
         * \name Collision and Streaming for direction 7.
         *
         * \brief Solves the LB equation.
         *
         * \param result The destination matrix.
         * \param dist The temporary distribution matrix.
         * \param eq_dist The equilibrium distribution matrix..
         * \param s_x Source matrix in x direction.
         * \param s_y Source matrix in y direction..
         * \param e_x Corresponding distribution scalar.
         * \param e_y Corresponding distribution scalar.
         * \param tau The relaxation time.
         */
        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_>& result,
                          DenseMatrix<DT1_>& dist,
                          DenseMatrix<DT1_>& eq_dist,
                          DenseMatrix<DT1_>& s_x,
                          DenseMatrix<DT1_>& s_y,
                          DT2_ e_x,
                          DT2_ e_y,
                          DT2_ tau)
        {

            CONTEXT("When performing collision and streaming in DIR 7:");
            long y_max(result.rows());
            long x_max(result.columns());

            for(long i(0); i < y_max; ++i)
            {
                long i_forward(i + 1);
                long i_backward(i - 1);

                for(long j(0); j < x_max; ++j)
                {
                    long j_forward(j + 1);
                    long j_backward(j - 1);

                    ///Respect periodic boundaries:
                    if(j_forward >= x_max)
                        j_forward = j_forward - x_max;
                    if(j_backward < 0)
                        j_backward = j_backward + x_max;
                    if(i_forward >= y_max)
                        i_forward = i_forward - y_max;
                    if(i_backward < 0)
                        i_backward = i_backward + y_max;

                    ///Perform streaming and collision:
                    result(i_backward, j) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                }
            }
        }
    };

    /**
     * \brief Collision and streaming module for LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::PERIODIC, lbm_lattice_types::D2Q9::DIR_8>
    {
        /**
         * \name Collision and Streaming for direction 8.
         *
         * \brief Solves the LB equation.
         *
         * \param result The destination matrix.
         * \param dist The temporary distribution matrix.
         * \param eq_dist The equilibrium distribution matrix..
         * \param s_x Source matrix in x direction.
         * \param s_y Source matrix in y direction..
         * \param e_x Corresponding distribution scalar.
         * \param e_y Corresponding distribution scalar.
         * \param tau The relaxation time.
         */
        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_>& result,
                          DenseMatrix<DT1_>& dist,
                          DenseMatrix<DT1_>& eq_dist,
                          DenseMatrix<DT1_>& s_x,
                          DenseMatrix<DT1_>& s_y,
                          DT2_ e_x,
                          DT2_ e_y,
                          DT2_ tau)
        {

            CONTEXT("When performing collision and streaming in DIR 8:");
            long y_max(result.rows());
            long x_max(result.columns());

            for(long i(0); i < y_max; ++i)
            {
                long i_forward(i + 1);
                long i_backward(i - 1);

                for(long j(0); j < x_max; ++j)
                {
                    long j_forward(j + 1);
                    long j_backward(j - 1);

                    ///Respect periodic boundaries:
                    if(j_forward >= x_max)
                        j_forward = j_forward - x_max;
                    if(j_backward < 0)
                        j_backward = j_backward + x_max;
                    if(i_forward >= y_max)
                        i_forward = i_forward - y_max;
                    if(i_backward < 0)
                        i_backward = i_backward + y_max;

                    ///Perform streaming and collision:
                    result(i_backward, j_forward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                }
            }
        }
    };

    /**
     * \brief Collision and streaming module for LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::PERIODIC, lbm_lattice_types::D2Q9::DIR_0>
    {
        /**
         * \name Collision and Streaming for direction 0.
         *
         * \brief Solves the LB equation.
         *
         * \param result The destination matrix.
         * \param dist The temporary distribution matrix.
         * \param eq_dist The equilibrium distribution matrix..
         * \param s_x Source matrix in x direction.
         * \param s_y Source matrix in y direction..
         * \param e_x Corresponding distribution scalar.
         * \param e_y Corresponding distribution scalar.
         * \param tau The relaxation time.
         */
        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_>& result,
                          DenseMatrix<DT1_>& dist,
                          DenseMatrix<DT1_>& eq_dist,
                          HONEI_UNUSED DenseMatrix<DT1_>& s_x,
                          HONEI_UNUSED DenseMatrix<DT1_>& s_y,
                          HONEI_UNUSED DT2_ e_x,
                          HONEI_UNUSED DT2_ e_y,
                          DT2_ tau)
        {

            CONTEXT("When performing collision and streaming in DIR 0:");
            long y_max(result.rows());
            long x_max(result.columns());

            for(long i(0); i < y_max; ++i)
            {
                for(long j(0); j < x_max; ++j)
                {
                   ///Perform streaming and collision:
                    result(i, j) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau;
                }
            }
        }
    };

//---------------------------------------------------------------------------------------------------------------
    /**
     * \brief Collision and streaming module for LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP_2, lbm_lattice_types::D2Q9::DIR_1>
    {
        /**
         * \name Collision and Streaming for direction 1..
         *
         * \brief Solves the LB equation.
         *
         * \param result The destination matrix.
         * \param dist The temporary distribution matrix.
         * \param eq_dist The equilibrium distribution matrix..
         * \param s_x Source matrix in x direction.
         * \param s_y Source matrix in y direction..
         * \param e_x Corresponding distribution scalar.
         * \param e_y Corresponding distribution scalar.
         * \param tau The relaxation time.
         */
        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_>& result,
                          DenseMatrix<DT1_>& dist,
                          DenseMatrix<DT1_>& eq_dist,
                          DenseMatrix<DT1_>& s_x,
                          DenseMatrix<DT1_>& s_y,
                          DT2_ e_x,
                          DT2_ e_y,
                          DT2_ tau)
        {
            CONTEXT("When performing collision and streaming in DIR 1:");
            long y_max(result.rows());
            long x_max(result.columns());

            for(long i(0); i < y_max; ++i)
            {
                for(long j(0); j < x_max; ++j)
                {
                    long j_forward(j + 1);
                    long j_backward(j - 1);

                    if(j_forward >= x_max)
                        result(i,j_backward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                    else
                        result(i,j_forward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                }
            }
        }
    };

    /**
     * \brief Collision and streaming module for LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP_2, lbm_lattice_types::D2Q9::DIR_2>
    {
        /**
         * \name Collision and Streaming for direction 2.
         *
         * \brief Solves the LB equation.
         *
         * \param result The destination matrix.
         * \param dist The temporary distribution matrix.
         * \param eq_dist The equilibrium distribution matrix..
         * \param s_x Source matrix in x direction.
         * \param s_y Source matrix in y direction..
         * \param e_x Corresponding distribution scalar.
         * \param e_y Corresponding distribution scalar.
         * \param tau The relaxation time.
         */
        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_>& result,
                          DenseMatrix<DT1_>& dist,
                          DenseMatrix<DT1_>& eq_dist,
                          DenseMatrix<DT1_>& s_x,
                          DenseMatrix<DT1_>& s_y,
                          DT2_ e_x,
                          DT2_ e_y,
                          DT2_ tau)
        {
            CONTEXT("When performing collision and streaming in DIR 2:");
            long y_max(result.rows());
            long x_max(result.columns());

            for(long i(0); i < y_max; ++i)
            {
                long i_forward(i + 1);
                long i_backward(i - 1);

                for(long j(0); j < x_max; ++j)
                {
                    long j_forward(j + 1);
                    long j_backward(j - 1);

                    if(i_forward >= y_max || j_forward >= x_max)
                    {
                        if(i_backward >= 0 && j_backward >= 0)
                            result(i_backward, j_backward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                        else
                        {
                            if(j_backward < 0)
                                result(i_backward, j_forward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                            else
                                result(i_forward, j_backward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                        }
                    }
                    else
                        result(i_forward, j_forward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                }
            }
        }
    };

    /**
     * \brief Collision and streaming module for LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP_2, lbm_lattice_types::D2Q9::DIR_3>
    {
        /**
         * \name Collision and Streaming for direction 3.
         *
         * \brief Solves the LB equation.
         *
         * \param result The destination matrix.
         * \param dist The temporary distribution matrix.
         * \param eq_dist The equilibrium distribution matrix..
         * \param s_x Source matrix in x direction.
         * \param s_y Source matrix in y direction..
         * \param e_x Corresponding distribution scalar.
         * \param e_y Corresponding distribution scalar.
         * \param tau The relaxation time.
         */
        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_>& result,
                          DenseMatrix<DT1_>& dist,
                          DenseMatrix<DT1_>& eq_dist,
                          DenseMatrix<DT1_>& s_x,
                          DenseMatrix<DT1_>& s_y,
                          DT2_ e_x,
                          DT2_ e_y,
                          DT2_ tau)
        {

            CONTEXT("When performing collision and streaming in DIR 3:");
            long y_max(result.rows());
            long x_max(result.columns());

            for(long i(0); i < y_max; ++i)
            {
                long i_forward(i + 1);
                long i_backward(i - 1);

                for(long j(0); j < x_max; ++j)
                {
                    if(i_forward >= y_max)
                        result(i_backward, j) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                    else
                        result(i_forward, j) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                }
            }
        }
    };

    /**
     * \brief Collision and streaming module for LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP_2, lbm_lattice_types::D2Q9::DIR_4>
    {
        /**
         * \name Collision and Streaming for direction 4.
         *
         * \brief Solves the LB equation.
         *
         * \param result The destination matrix.
         * \param dist The temporary distribution matrix.
         * \param eq_dist The equilibrium distribution matrix..
         * \param s_x Source matrix in x direction.
         * \param s_y Source matrix in y direction..
         * \param e_x Corresponding distribution scalar.
         * \param e_y Corresponding distribution scalar.
         * \param tau The relaxation time.
         */
        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_>& result,
                          DenseMatrix<DT1_>& dist,
                          DenseMatrix<DT1_>& eq_dist,
                          DenseMatrix<DT1_>& s_x,
                          DenseMatrix<DT1_>& s_y,
                          DT2_ e_x,
                          DT2_ e_y,
                          DT2_ tau)
        {

            CONTEXT("When performing collision and streaming in DIR 4:");
            long y_max(result.rows());
            long x_max(result.columns());

            for(long i(0); i < y_max; ++i)
            {
                long i_forward(i + 1);
                long i_backward(i - 1);

                for(long j(0); j < x_max; ++j)
                {
                    long j_forward(j + 1);
                    long j_backward(j - 1);

                    if(i_forward >= y_max || j_backward < 0)
                    {
                        if(i_backward >= 0 && j_forward < x_max)
                            result(i_backward, j_forward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                        else
                        {
                            if(i_backward < 0)
                                result(i_forward, j_forward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                            else

                                result(i_backward, j_backward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                        }
                    }
                    else
                        result(i_forward, j_backward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                }
            }
        }
    };

    /**
     * \brief Collision and streaming module for LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP_2, lbm_lattice_types::D2Q9::DIR_5>
    {
        /**
         * \name Collision and Streaming for direction 5.
         *
         * \brief Solves the LB equation.
         *
         * \param result The destination matrix.
         * \param dist The temporary distribution matrix.
         * \param eq_dist The equilibrium distribution matrix..
         * \param s_x Source matrix in x direction.
         * \param s_y Source matrix in y direction..
         * \param e_x Corresponding distribution scalar.
         * \param e_y Corresponding distribution scalar.
         * \param tau The relaxation time.
         */
        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_>& result,
                          DenseMatrix<DT1_>& dist,
                          DenseMatrix<DT1_>& eq_dist,
                          DenseMatrix<DT1_>& s_x,
                          DenseMatrix<DT1_>& s_y,
                          DT2_ e_x,
                          DT2_ e_y,
                          DT2_ tau)
        {

            CONTEXT("When performing collision and streaming in DIR 5:");
            long y_max(result.rows());
            long x_max(result.columns());

            for(long i(0); i < y_max; ++i)
            {
                for(long j(0); j < x_max; ++j)
                {
                    long j_forward(j + 1);
                    long j_backward(j - 1);

                    if(j_backward < 0)
                        result(i, j_forward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                    else
                        result(i, j_backward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                }
            }
        }
    };

    /**
     * \brief Collision and streaming module for LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP_2, lbm_lattice_types::D2Q9::DIR_6>
    {
        /**
         * \name Collision and Streaming for direction 6.
         *
         * \brief Solves the LB equation.
         *
         * \param result The destination matrix.
         * \param dist The temporary distribution matrix.
         * \param eq_dist The equilibrium distribution matrix..
         * \param s_x Source matrix in x direction.
         * \param s_y Source matrix in y direction..
         * \param e_x Corresponding distribution scalar.
         * \param e_y Corresponding distribution scalar.
         * \param tau The relaxation time.
         */
        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_>& result,
                          DenseMatrix<DT1_>& dist,
                          DenseMatrix<DT1_>& eq_dist,
                          DenseMatrix<DT1_>& s_x,
                          DenseMatrix<DT1_>& s_y,
                          DT2_ e_x,
                          DT2_ e_y,
                          DT2_ tau)
        {

            CONTEXT("When performing collision and streaming in DIR 6:");
            long y_max(result.rows());
            long x_max(result.columns());

            for(long i(0); i < y_max; ++i)
            {
                long i_forward(i + 1);
                long i_backward(i - 1);

                for(long j(0); j < x_max; ++j)
                {
                    long j_forward(j + 1);
                    long j_backward(j - 1);

                    if(i_backward < 0 || j_backward < 0)
                    {
                        if(i_forward < y_max && j_forward < y_max)
                            result(i_forward, j_forward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                        else
                        {
                            if(i_forward >= y_max)
                                result(i_backward, j_forward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                            else
                                result(i_forward, j_backward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                        }
                    }
                    else
                        result(i_backward, j_backward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                }
            }
        }
    };

    /**
     * \brief Collision and streaming module for LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP_2, lbm_lattice_types::D2Q9::DIR_7>
    {
        /**
         * \name Collision and Streaming for direction 7.
         *
         * \brief Solves the LB equation.
         *
         * \param result The destination matrix.
         * \param dist The temporary distribution matrix.
         * \param eq_dist The equilibrium distribution matrix..
         * \param s_x Source matrix in x direction.
         * \param s_y Source matrix in y direction..
         * \param e_x Corresponding distribution scalar.
         * \param e_y Corresponding distribution scalar.
         * \param tau The relaxation time.
         */
        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_>& result,
                          DenseMatrix<DT1_>& dist,
                          DenseMatrix<DT1_>& eq_dist,
                          DenseMatrix<DT1_>& s_x,
                          DenseMatrix<DT1_>& s_y,
                          DT2_ e_x,
                          DT2_ e_y,
                          DT2_ tau)
        {

            CONTEXT("When performing collision and streaming in DIR 7:");
            long y_max(result.rows());
            long x_max(result.columns());

            for(long i(0); i < y_max; ++i)
            {
                long i_forward(i + 1);
                long i_backward(i - 1);

                for(long j(0); j < x_max; ++j)
                {
                    if(i_backward < 0)
                        result(i_forward, j) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                    else
                        result(i_backward, j) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                }
            }
        }
    };

    /**
     * \brief Collision and streaming module for LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP_2, lbm_lattice_types::D2Q9::DIR_8>
    {
        /**
         * \name Collision and Streaming for direction 8.
         *
         * \brief Solves the LB equation.
         *
         * \param result The destination matrix.
         * \param dist The temporary distribution matrix.
         * \param eq_dist The equilibrium distribution matrix..
         * \param s_x Source matrix in x direction.
         * \param s_y Source matrix in y direction..
         * \param e_x Corresponding distribution scalar.
         * \param e_y Corresponding distribution scalar.
         * \param tau The relaxation time.
         */
        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_>& result,
                          DenseMatrix<DT1_>& dist,
                          DenseMatrix<DT1_>& eq_dist,
                          DenseMatrix<DT1_>& s_x,
                          DenseMatrix<DT1_>& s_y,
                          DT2_ e_x,
                          DT2_ e_y,
                          DT2_ tau)
        {

            CONTEXT("When performing collision and streaming in DIR 8:");
            long y_max(result.rows());
            long x_max(result.columns());

            for(long i(0); i < y_max; ++i)
            {
                long i_forward(i + 1);
                long i_backward(i - 1);

                for(long j(0); j < x_max; ++j)
                {
                    long j_forward(j + 1);
                    long j_backward(j - 1);

                    if(i_backward < 0 || j_forward >= x_max)
                    {
                        if(i_forward < y_max && j_backward >= 0)
                            result(i_forward, j_backward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                        else
                        {
                            if(i_forward >= y_max)
                                result(i_backward, j_backward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                            else
                                result(i_forward, j_forward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                        }
                    }
                    else
                        result(i_backward, j_forward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                }
            }
        }
    };

    /**
     * \brief Collision and streaming module for LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP_2, lbm_lattice_types::D2Q9::DIR_0>
    {
        /**
         * \name Collision and Streaming for direction 0.
         *
         * \brief Solves the LB equation.
         *
         * \param result The destination matrix.
         * \param dist The temporary distribution matrix.
         * \param eq_dist The equilibrium distribution matrix..
         * \param s_x Source matrix in x direction.
         * \param s_y Source matrix in y direction..
         * \param e_x Corresponding distribution scalar.
         * \param e_y Corresponding distribution scalar.
         * \param tau The relaxation time.
         */
        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_>& result,
                          DenseMatrix<DT1_>& dist,
                          DenseMatrix<DT1_>& eq_dist,
                          DenseMatrix<DT1_>& /*s_x*/,
                          DenseMatrix<DT1_>& /*s_y*/,
                          DT2_ /*e_x*/,
                          DT2_ /*e_y*/,
                          DT2_ tau)
        {

            CONTEXT("When performing collision and streaming in DIR 0:");
            long y_max(result.rows());
            long x_max(result.columns());

            for(long i(0); i < y_max; ++i)
            {
                for(long j(0); j < x_max; ++j)
                {
                   ///Perform streaming and collision:
                    result(i, j) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau;
                }
            }
        }
    };
//---------------------------------------------------------------------------------------------------------------
    /**
     * \brief Collision and streaming module for LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9::DIR_1>
    {
        /**
         * \name Collision and Streaming for direction 1..
         *
         * \brief Solves the LB equation.
         *
         * \param result The destination matrix.
         * \param dist The temporary distribution matrix.
         * \param eq_dist The equilibrium distribution matrix..
         * \param s_x Source matrix in x direction.
         * \param s_y Source matrix in y direction..
         * \param e_x Corresponding distribution scalar.
         * \param e_y Corresponding distribution scalar.
         * \param tau The relaxation time.
         */
        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_>& result,
                          DenseMatrix<DT1_>& dist,
                          DenseMatrix<DT1_>& eq_dist,
                          DenseMatrix<DT1_>& s_x,
                          DenseMatrix<DT1_>& s_y,
                          DT2_ e_x,
                          DT2_ e_y,
                          DT2_ tau)
        {
            CONTEXT("When performing collision and streaming in DIR 1:");
            long y_max(result.rows());
            long x_max(result.columns());

            for(long i(0); i < y_max; ++i)
            {
                //long i_forward(i + 1);
                //long i_backward(i - 1);

                for(long j(0); j < x_max; ++j)
                {
                    long j_forward(j + 1);
                    //long j_backward(j - 1);

                    if(j_forward >= x_max)
                        ;//result(i,j) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                    else
                        result(i,j_forward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                }
            }
        }
    };

    /**
     * \brief Collision and streaming module for LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9::DIR_2>
    {
        /**
         * \name Collision and Streaming for direction 2.
         *
         * \brief Solves the LB equation.
         *
         * \param result The destination matrix.
         * \param dist The temporary distribution matrix.
         * \param eq_dist The equilibrium distribution matrix..
         * \param s_x Source matrix in x direction.
         * \param s_y Source matrix in y direction..
         * \param e_x Corresponding distribution scalar.
         * \param e_y Corresponding distribution scalar.
         * \param tau The relaxation time.
         */
        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_>& result,
                          DenseMatrix<DT1_>& dist,
                          DenseMatrix<DT1_>& eq_dist,
                          DenseMatrix<DT1_>& s_x,
                          DenseMatrix<DT1_>& s_y,
                          DT2_ e_x,
                          DT2_ e_y,
                          DT2_ tau)
        {
            CONTEXT("When performing collision and streaming in DIR 2:");
            long y_max(result.rows());
            long x_max(result.columns());

            for(long i(0); i < y_max; ++i)
            {
                long i_forward(i + 1);
                //long i_backward(i - 1);

                for(long j(0); j < x_max; ++j)
                {
                    long j_forward(j + 1);
                    //long j_backward(j - 1);

                    if(i_forward >= y_max || j_forward >= x_max)
                        ;//result(i, j) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                    else
                        result(i_forward, j_forward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                }
            }
        }
    };

    /**
     * \brief Collision and streaming module for LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9::DIR_3>
    {
        /**
         * \name Collision and Streaming for direction 3.
         *
         * \brief Solves the LB equation.
         *
         * \param result The destination matrix.
         * \param dist The temporary distribution matrix.
         * \param eq_dist The equilibrium distribution matrix..
         * \param s_x Source matrix in x direction.
         * \param s_y Source matrix in y direction..
         * \param e_x Corresponding distribution scalar.
         * \param e_y Corresponding distribution scalar.
         * \param tau The relaxation time.
         */
        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_>& result,
                          DenseMatrix<DT1_>& dist,
                          DenseMatrix<DT1_>& eq_dist,
                          DenseMatrix<DT1_>& s_x,
                          DenseMatrix<DT1_>& s_y,
                          DT2_ e_x,
                          DT2_ e_y,
                          DT2_ tau)
        {

            CONTEXT("When performing collision and streaming in DIR 3:");
            long y_max(result.rows());
            long x_max(result.columns());

            for(long i(0); i < y_max; ++i)
            {
                long i_forward(i + 1);
                //long i_backward(i - 1);

                for(long j(0); j < x_max; ++j)
                {
                    //long j_forward(j + 1);
                    //long j_backward(j - 1);

                    if(i_forward >= y_max)
                        ;//result(i, j) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                    else
                        result(i_forward, j) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                }
            }
        }
    };

    /**
     * \brief Collision and streaming module for LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9::DIR_4>
    {
        /**
         * \name Collision and Streaming for direction 4.
         *
         * \brief Solves the LB equation.
         *
         * \param result The destination matrix.
         * \param dist The temporary distribution matrix.
         * \param eq_dist The equilibrium distribution matrix..
         * \param s_x Source matrix in x direction.
         * \param s_y Source matrix in y direction..
         * \param e_x Corresponding distribution scalar.
         * \param e_y Corresponding distribution scalar.
         * \param tau The relaxation time.
         */
        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_>& result,
                          DenseMatrix<DT1_>& dist,
                          DenseMatrix<DT1_>& eq_dist,
                          DenseMatrix<DT1_>& s_x,
                          DenseMatrix<DT1_>& s_y,
                          DT2_ e_x,
                          DT2_ e_y,
                          DT2_ tau)
        {

            CONTEXT("When performing collision and streaming in DIR 4:");
            long y_max(result.rows());
            long x_max(result.columns());

            for(long i(0); i < y_max; ++i)
            {
                long i_forward(i + 1);
                //long i_backward(i - 1);

                for(long j(0); j < x_max; ++j)
                {
                    //long j_forward(j + 1);
                    long j_backward(j - 1);

                    if(i_forward >= y_max || j_backward < 0)
                        ;//result(i, j) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                    else
                        result(i_forward, j_backward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                }
            }
        }
    };

    /**
     * \brief Collision and streaming module for LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9::DIR_5>
    {
        /**
         * \name Collision and Streaming for direction 5.
         *
         * \brief Solves the LB equation.
         *
         * \param result The destination matrix.
         * \param dist The temporary distribution matrix.
         * \param eq_dist The equilibrium distribution matrix..
         * \param s_x Source matrix in x direction.
         * \param s_y Source matrix in y direction..
         * \param e_x Corresponding distribution scalar.
         * \param e_y Corresponding distribution scalar.
         * \param tau The relaxation time.
         */
        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_>& result,
                          DenseMatrix<DT1_>& dist,
                          DenseMatrix<DT1_>& eq_dist,
                          DenseMatrix<DT1_>& s_x,
                          DenseMatrix<DT1_>& s_y,
                          DT2_ e_x,
                          DT2_ e_y,
                          DT2_ tau)
        {

            CONTEXT("When performing collision and streaming in DIR 5:");
            long y_max(result.rows());
            long x_max(result.columns());

            for(long i(0); i < y_max; ++i)
            {
                //long i_forward(i + 1);
                //long i_backward(i - 1);

                for(long j(0); j < x_max; ++j)
                {
                    //long j_forward(j + 1);
                    long j_backward(j - 1);

                    if(j_backward < 0)
                        ;//result(i, j) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                    else
                        result(i, j_backward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                }
            }
        }
    };

    /**
     * \brief Collision and streaming module for LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9::DIR_6>
    {
        /**
         * \name Collision and Streaming for direction 6.
         *
         * \brief Solves the LB equation.
         *
         * \param result The destination matrix.
         * \param dist The temporary distribution matrix.
         * \param eq_dist The equilibrium distribution matrix..
         * \param s_x Source matrix in x direction.
         * \param s_y Source matrix in y direction..
         * \param e_x Corresponding distribution scalar.
         * \param e_y Corresponding distribution scalar.
         * \param tau The relaxation time.
         */
        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_>& result,
                          DenseMatrix<DT1_>& dist,
                          DenseMatrix<DT1_>& eq_dist,
                          DenseMatrix<DT1_>& s_x,
                          DenseMatrix<DT1_>& s_y,
                          DT2_ e_x,
                          DT2_ e_y,
                          DT2_ tau)
        {

            CONTEXT("When performing collision and streaming in DIR 6:");
            long y_max(result.rows());
            long x_max(result.columns());

            for(long i(0); i < y_max; ++i)
            {
                //long i_forward(i + 1);
                long i_backward(i - 1);

                for(long j(0); j < x_max; ++j)
                {
                    //long j_forward(j + 1);
                    long j_backward(j - 1);

                    if(i_backward < 0 || j_backward < 0)
                        ;//result(i, j) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                    else
                        result(i_backward, j_backward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                }
            }
        }
    };

    /**
     * \brief Collision and streaming module for LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9::DIR_7>
    {
        /**
         * \name Collision and Streaming for direction 7.
         *
         * \brief Solves the LB equation.
         *
         * \param result The destination matrix.
         * \param dist The temporary distribution matrix.
         * \param eq_dist The equilibrium distribution matrix..
         * \param s_x Source matrix in x direction.
         * \param s_y Source matrix in y direction..
         * \param e_x Corresponding distribution scalar.
         * \param e_y Corresponding distribution scalar.
         * \param tau The relaxation time.
         */
        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_>& result,
                          DenseMatrix<DT1_>& dist,
                          DenseMatrix<DT1_>& eq_dist,
                          DenseMatrix<DT1_>& s_x,
                          DenseMatrix<DT1_>& s_y,
                          DT2_ e_x,
                          DT2_ e_y,
                          DT2_ tau)
        {

            CONTEXT("When performing collision and streaming in DIR 7:");
            long y_max(result.rows());
            long x_max(result.columns());

            for(long i(0); i < y_max; ++i)
            {
                //long i_forward(i + 1);
                long i_backward(i - 1);

                for(long j(0); j < x_max; ++j)
                {
                    //long j_forward(j + 1);
                    //long j_backward(j - 1);

                    if(i_backward < 0)
                        ;//result(i, j) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                    else
                        result(i_backward, j) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                }
            }
        }
    };

    /**
     * \brief Collision and streaming module for LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9::DIR_8>
    {
        /**
         * \name Collision and Streaming for direction 8.
         *
         * \brief Solves the LB equation.
         *
         * \param result The destination matrix.
         * \param dist The temporary distribution matrix.
         * \param eq_dist The equilibrium distribution matrix..
         * \param s_x Source matrix in x direction.
         * \param s_y Source matrix in y direction..
         * \param e_x Corresponding distribution scalar.
         * \param e_y Corresponding distribution scalar.
         * \param tau The relaxation time.
         */
        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_>& result,
                          DenseMatrix<DT1_>& dist,
                          DenseMatrix<DT1_>& eq_dist,
                          DenseMatrix<DT1_>& s_x,
                          DenseMatrix<DT1_>& s_y,
                          DT2_ e_x,
                          DT2_ e_y,
                          DT2_ tau)
        {

            CONTEXT("When performing collision and streaming in DIR 8:");
            long y_max(result.rows());
            long x_max(result.columns());

            for(long i(0); i < y_max; ++i)
            {
                //long i_forward(i + 1);
                long i_backward(i - 1);

                for(long j(0); j < x_max; ++j)
                {
                    long j_forward(j + 1);
                    //long j_backward(j - 1);

                    if(i_backward < 0 || j_forward >= x_max)
                        ;//result(i, j) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                    else
                        result(i_backward, j_forward) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau + DT1_(1./6.) * (e_x * s_x(i,j) + e_y * s_y(i,j));
                }
            }
        }
    };

    /**
     * \brief Collision and streaming module for LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct CollideStream<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9::DIR_0>
    {
        /**
         * \name Collision and Streaming for direction 0.
         *
         * \brief Solves the LB equation.
         *
         * \param result The destination matrix.
         * \param dist The temporary distribution matrix.
         * \param eq_dist The equilibrium distribution matrix..
         * \param s_x Source matrix in x direction.
         * \param s_y Source matrix in y direction..
         * \param e_x Corresponding distribution scalar.
         * \param e_y Corresponding distribution scalar.
         * \param tau The relaxation time.
         */
        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_>& result,
                          DenseMatrix<DT1_>& dist,
                          DenseMatrix<DT1_>& eq_dist,
                          HONEI_UNUSED DenseMatrix<DT1_>& s_x,
                          HONEI_UNUSED DenseMatrix<DT1_>& s_y,
                          HONEI_UNUSED DT2_ e_x,
                          HONEI_UNUSED DT2_ e_y,
                          DT2_ tau)
        {

            CONTEXT("When performing collision and streaming in DIR 0:");
            long y_max(result.rows());
            long x_max(result.columns());

            for(long i(0); i < y_max; ++i)
            {
                for(long j(0); j < x_max; ++j)
                {
                   ///Perform streaming and collision:
                    result(i, j) = dist(i,j) - (dist(i,j) - eq_dist(i,j))/tau;
                }
            }
        }
    };

}
#endif
