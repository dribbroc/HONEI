/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the SWE C++ library. LibSWE is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibSWE is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef LIBSWE_GUARD_SCENARIO_HH
#define LIBSWE_GUARD_SCENARIO_HH 1

//#include <libswe/implicit_solver.hh>
#include <libswe/boundary_types.hh>
#include <libmath/methods.hh>


namespace honei
{
    namespace swe_solvers
    {
        class RELAX
        {
        };

        class IMPLICIT
        {
        };
    }

    using namespace swe_solvers;
    template<typename ResPrec_, typename SWESolver_, typename BoundaryType_>
        class Scenario
        {
        };

    template<typename ResPrec_>
        class Scenario<ResPrec_, swe_solvers::IMPLICIT, boundaries::REFLECT>
        {

            private:
                ///Flags for validation by ScenarioManager:
                bool scalarfields_set;
                bool boundaries_set;
                bool delta_t_set;

                bool allocation;

            public:
                ///Our Problem size.
                unsigned long n;

                ///Stepsize in x direction.
                ResPrec_ delta_x;
                ///Stepsize in y direction.
                ResPrec_ delta_y;
                ///Size of timestep.
                ResPrec_ delta_t;

                ///Dimensions of OMEGA:
                ResPrec_ d_width;
                ResPrec_ d_height;
                unsigned long grid_width;
                unsigned long grid_height;

                ///The input- and to-be-updated - data.
                DenseMatrix<ResPrec_> * bottom;
                DenseMatrix<ResPrec_> * height;
                DenseMatrix<ResPrec_> * x_veloc;
                DenseMatrix<ResPrec_> * y_veloc;

                ///The data to work on.
                BandedMatrix<ResPrec_> * system_matrix;
                DenseVector<ResPrec_> * right_hand_side;
                DenseVector<ResPrec_> * u_temp;
                DenseVector<ResPrec_> * v_temp;

                ///The boundary maps of the scalarfields:
                DenseMatrix<ResPrec_>* height_bound;
                DenseMatrix<ResPrec_>* bottom_bound;
                DenseMatrix<ResPrec_>* x_veloc_bound;
                DenseMatrix<ResPrec_>* y_veloc_bound;


                /**
                 * Constructor for square shaped OMEGA.
                 *
                 * \param size Our problem size.
                 *
                 **/
                Scenario(unsigned long size)
                {
                    n = size;
                    grid_width = n;
                    grid_height = n;
                }

                /**
                 * Constructor for rectangular OMEGA.
                 *
                 * \param d_width The width of OMEGA (OMEGA x parameter interval range).
                 * \param d_height The height of OMEGA (OMEGA y parameter interval range).
                 * \param delta_x The stepsize in x direction.
                 * \param delta_y The stepsize in y direction.
                 *
                 **/
                Scenario(ResPrec_ dwidth, ResPrec_ dheight, ResPrec_ deltax, ResPrec_ deltay)
                {
                    d_width = dwidth;
                    d_height = dheight;
                    delta_x = deltax;
                    delta_y = deltay;

                    //\TODO: catch error: d_width%delta_x !=0
                    grid_width = d_width/delta_x;
                    grid_height = d_height/delta_y;
                }

                /**
                 * Constructor for rectangular grid.
                 *
             * \param grid_width The width of the grid.
             * \param grid_height The height of the grid.
             *
             **/
            Scenario(unsigned long gridwidth, unsigned long gridheight)
            {
                grid_width = gridwidth;
                grid_height = gridheight;
            }
    };

    template<typename ResPrec_>
    class Scenario<ResPrec_, swe_solvers::RELAX, boundaries::REFLECT>
    {
        public:
            ///Stepsize in x direction.
            ResPrec_ delta_x;
            ///Stepsize in y direction.
            ResPrec_ delta_y;
            ///Size of timestep.
            ResPrec_ delta_t;
            ///Current timestep.
            unsigned int solve_time;

            ///Relaxation parameter epsilon.
            ResPrec_ eps;

            ///Dimensions of OMEGA:
            unsigned long d_width;
            unsigned long d_height;

            ///The input- and to-be-updated - data.
            DenseMatrix<ResPrec_> * bottom;
            DenseMatrix<ResPrec_> * height;
            DenseMatrix<ResPrec_> * x_veloc;
            DenseMatrix<ResPrec_> * y_veloc;

            ///The data to work on
            DenseVector<ResPrec_> * c;
            DenseVector<ResPrec_> * d;

            DenseVector<ResPrec_> * u;
            DenseVector<ResPrec_> * v;
            DenseVector<ResPrec_> * w;

            ///The boundary maps of the scalarfields:
            DenseMatrix<ResPrec_>* height_bound;
            DenseMatrix<ResPrec_>* bottom_bound;
            DenseMatrix<ResPrec_>* x_veloc_bound;
            DenseMatrix<ResPrec_>* y_veloc_bound;

            ///Vectors for the bottom slopes.
            DenseVector<ResPrec_> * bottom_slopes_x;
            DenseVector<ResPrec_> * bottom_slopes_y;

            ///Manning constant:
            ResPrec_ manning_n;

            /**
             * Constructor for rectangular grid.
             *
             * \param grid_width The width of the grid.
             * \param grid_height The height of the grid.
             *
             **/
            Scenario(unsigned long gridwidth, unsigned long gridheight)
            {
                d_width = gridwidth;
                d_height = gridheight;
            }

    };

}

#endif
