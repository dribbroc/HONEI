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

#include <libswe/implicit_solver.hh>
#include <libswe/boundary_types.hh>
#include <libmath/methods.hh>

namespace swe_solvers
{
    class RELAX
    {
    };

    class IMPLICIT
    {
    };
}
namespace honei
{

    template<typename ResPrec_, typename SWESolver_, typename BoundaryType_>
    class Scenario
    {
    };

    template<typename ResPrec_>
    class Scenario<ResPrec_, swe_solvers::IMPLICIT, boundaries::REFLECT>
    {

        public:
            ///Our Problem size.
            unsigned long _n;

            ///Stepsize in x direction.
            ResPrec_ _delta_x;
            ///Stepsize in y direction.
            ResPrec_ _delta_y;
            ///Size of timestep.
            ResPrec_ _delta_t;
            ///Current timestep.
            unsigned int _solve_time;

            ///Dimensions of OMEGA:
            ResPrec_ _d_width;
            ResPrec_ _d_height;
            unsigned long _grid_width;
            unsigned long _grid_height;

            ///The input- and to-be-updated - data.
            DenseMatrix<ResPrec_> * _bottom;
            DenseMatrix<ResPrec_> * _height;
            DenseMatrix<ResPrec_> * _x_veloc;
            DenseMatrix<ResPrec_> * _y_veloc;

            ///The data to work on.
            DenseMatrix<ResPrec_> * _system_matrix;
            DenseVector<ResPrec_> * _right_hand_side;

            /**
             * Constructor for square shaped OMEGA.
             *
             * \param size Our problem size.
             *
             **/
            Scenario(unsigned long size)
            {
                _n = size;
                _grid_width = _n;
                _grid_height = _n;
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
            Scenario(ResPrec_ d_width, ResPrec_ d_height, ResPrec_ delta_x, ResPrec_ delta_y)
            {
                _d_width = d_width;
                _d_height = d_height;
                _delta_x = delta_x;
                _delta_y = delta_y;

                //\TODO: catch error: d_width%delta_x !=0
                _grid_width = _d_width/_delta_x;
                _grid_height = _d_height/_delta_y;
            }

            /**
             * Constructor for rectangular grid.
             *
             * \param grid_width The width of the grid.
             * \param grid_height The height of the grid.
             *
             **/
            Scenario(ResPrec_ grid_width, ResPrec_ grid_height)
            {
                _grid_width = grid_width;
                _grid_height = grid_height;
            }
    };

    template<typename ResPrec_>
    class Scenario<ResPrec_, swe_solvers::RELAX, boundaries::REFLECT>
    {
        public:
            ///Stepsize in x direction.
            ResPrec_ _delta_x;
            ///Stepsize in y direction.
            ResPrec_ _delta_y;
            ///Size of timestep.
            ResPrec_ _delta_t;
            ///Current timestep.
            unsigned int _solve_time;

            ///Dimensions of OMEGA:
            unsigned long _d_width;
            unsigned long _d_height;

            ///The input- and to-be-updated - data.
            DenseMatrix<ResPrec_> * _bottom;
            DenseMatrix<ResPrec_> * _height;
            DenseMatrix<ResPrec_> * _x_veloc;
            DenseMatrix<ResPrec_> * _y_veloc;

            ///The data to work on 
            // \TODO: c, d, etc...
    };

}

#endif
