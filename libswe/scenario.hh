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
        ///Stepsize in x direction.
        ResPrec_ _delta_x;
        ///Stepsize in y direction.
        ResPrec_ _delta_y;
        ///Size of timestep.
        ResPrec_ _delta_t;
        ///Current timestep.
        unsigned int _solve_time;

        ///Dimensions of OMEGA:
        double _d_width;
        double _d_height;

        ///The input- and to-be-updated - data.
        DenseMatrix<ResPrec_> * _bottom;
        DenseMatrix<ResPrec_> * _height;
        DenseMatrix<ResPrec_> * _x_veloc;
        DenseMatrix<ResPrec_> * _y_veloc;

        ///The data to work on.
        DenseMatrix<ResPrec_> * _system_matrix;
        DenseVector<ResPrec_> * _right_hand_side;

    };

    template<typename ResPrec_>
    class Scenario<ResPrec_, swe_solvers::RELAX, boundaries::REFLECT>
    {
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
