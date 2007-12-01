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

/**
 * \file
 *
 * Implementation of a semi-implicit solver for the 2D Shallow Water Equations.
 * Solver is arbitrary in precision. In contrast to the RelaxSolver, which can
 * be classified as a pipelined architecture, ImplicitSolver mainly consists of
 * the assembly and solution of a linear system. Due to this, mixed precision
 * can be obtained in an intuitive way, by using an iterative refinement solver
 * for the system. Therefore, you can specify the precision of the "inner"
 * workhorse by WorkPrec_ (system assembly).
 *
 * Since we have a finite differences - discetization here, we have a slightly
 * different point of view than with RelaxSolver. Either the stepsizes in the
 * two dimensions are determined by given input fields and the dimension sizes,
 * or the other way around. For the user it will be the most popular usecase
 * to simply have some inputfields and give the dimension sizes of OMEGA,
 * assuming (0,0) to be the lower left anchor of the rectangular region.
 *
 * \ingroup grplibswe
 **/

#ifndef LIBSWE_GUARD_IMPLICIT_SOLVER_HH
#define LIBSWE_GUARD_IMPLICIT_SOLVER_HH 1

#include <libla/dense_matrix.hh>
#include <libla/dense_vector.hh>
#include <libutil/tags.hh>
#include <libmath/methods.hh>
#include <libswe/boundary_types.hh>
#include <libswe/scenario.hh>
using namespace std;
using namespace methods;
using namespace swe_solvers;
using namespace boundaries;

namespace honei {
    template<typename Tag_, typename ResPrec_, typename SolverType_, typename BoundaryType_>
    class ImplicitSolver
    {
    };

    template<typename Tag_, typename ResPrec_>
    class ImplicitSolver<Tag_, ResPrec_, CG, REFLECT>
    {
        ///Private members:
        private:
            Scenario<ResPrec_, IMPLICIT, REFLECT>* scenario;
            /**
             * System assembly: A.
             **/
            template<typename WorkPrec_>
            void assemble_matrix();

            /**
             * System assembly: b.
             **/
            template<typename WorkPrec_>
            void assemble_right_hand_side();

        public:
            /**
             * Constructor.
             * */
            ImplicitSolver(Scenario<ResPrec_, IMPLICIT, REFLECT> & s)
            {
                scenario = &s;
            }

            /**
             * Solution capsule for one timestep.
             **/
            template<typename WorkPrec_>
            void solve(unsigned long iter_numbers);

            /**
             * Solution capsule for one timestep.
             **/
            template<typename WorkPrec_>
            void solve(double conv_rad);
    };
}
#endif
