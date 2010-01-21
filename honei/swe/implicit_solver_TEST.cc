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

//#define SOLVER_DEBUG 1
//#define SOLVER_VERBOSE 1
#ifdef DEBUG
#define SOLVER_VERBOSE
#endif
#include "implicit_solver.hh"
#include <unittest/unittest.hh>
#include <honei/util/stringify.hh>
#include <string>
#include <honei/swe/scenario.hh>

#include <sys/time.h>
#include <iostream>
using namespace honei;
using namespace tests;
using namespace std;
using namespace swe_solvers;
using namespace boundaries;

template <typename Tag_, typename DataType_>
class ImplicitSolverCreationTest :
    public BaseTest
{
    public:
        ImplicitSolverCreationTest(const std::string & type) :
            BaseTest("implicit_solver_creation_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            Scenario<DataType_, IMPLICIT, REFLECT> scenario(2000);
            ImplicitSolver<Tag_, DataType_, CG, REFLECT> solver(scenario);
            TEST_CHECK(true);
        }
};
template <typename Tag_, typename DataType_>
class ImplicitSolverPreprocessingTest :
    public BaseTest
{
    public:
        ImplicitSolverPreprocessingTest(const std::string & type) :
            BaseTest("implicit_solver_preprocessing_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            Scenario<DataType_, IMPLICIT, REFLECT> scenario(4);
            DenseMatrix<DataType_> h_b(6, 6, DataType_(0));
            DenseMatrix<DataType_> h(4 , 4, DataType_(1));
            DenseMatrix<DataType_> xv_b(6, 6, DataType_(0));
            DenseMatrix<DataType_> xv(4 , 4, DataType_(1));
            DenseMatrix<DataType_> yv_b(6, 6, DataType_(0));
            DenseMatrix<DataType_> yv(4 , 4, DataType_(1));
            DenseMatrix<DataType_> b_b(6, 6, DataType_(0));
            DenseMatrix<DataType_> b(4 , 4, DataType_(1));

            scenario.height = &h;
            scenario.x_veloc = &xv;
            scenario.y_veloc = &yv;
            scenario.bottom = &b;
            scenario.height_bound = &h_b;
            scenario.x_veloc_bound = &xv_b;
            scenario.y_veloc_bound = &yv_b;
            scenario.bottom_bound = &b_b;

            ImplicitSolver<Tag_, DataType_, CG, REFLECT> solver(scenario);
            solver.do_preprocessing();
            std::cout << *(scenario.height_bound) << std::endl;
            std::cout << *(scenario.x_veloc_bound) << std::endl;
            std::cout << *(scenario.y_veloc_bound) << std::endl;
            std::cout << *(scenario.bottom_bound)<< std::endl;
            TEST_CHECK(true);

        }
};
template <typename Tag_, typename DataType_>
class ImplicitSolverMatrixAssTest :
    public BaseTest
{
    public:
        ImplicitSolverMatrixAssTest(const std::string & type) :
            BaseTest("implicit_solver_matrix_ass_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            Scenario<DataType_, IMPLICIT, REFLECT> scenario(4);
            DenseMatrix<DataType_> h_b(6, 6, DataType_(0));
            DenseMatrix<DataType_> h(4 , 4, DataType_(1));
            DenseMatrix<DataType_> xv_b(6, 6, DataType_(0));
            DenseMatrix<DataType_> xv(4 , 4, DataType_(1));
            DenseMatrix<DataType_> yv_b(6, 6, DataType_(0));
            DenseMatrix<DataType_> yv(4 , 4, DataType_(1));
            DenseMatrix<DataType_> b_b(6, 6, DataType_(0));
            DenseMatrix<DataType_> b(4 , 4, DataType_(1));
            DenseVector<DataType_> u_t(36, DataType_(0));
            DenseVector<DataType_> v_t(36, DataType_(0));

            DenseVector<DataType_> rhs(36, DataType_(0));
            BandedMatrix<DataType_> A(36);
            scenario.height = &h;
            scenario.x_veloc = &xv;
            scenario.y_veloc = &yv;
            scenario.bottom = &b;
            scenario.height_bound = &h_b;
            scenario.x_veloc_bound = &xv_b;
            scenario.y_veloc_bound = &yv_b;
            scenario.bottom_bound = &b_b;
            scenario.system_matrix = &A;
            scenario.right_hand_side = &rhs;


            DataType_ delta_t(1);
            DataType_ delta_x(1);
            DataType_ delta_y(1);

            scenario.delta_t = delta_t;
            scenario.delta_x = delta_x;
            scenario.delta_y = delta_y;
            scenario.u_temp = &u_t;
            scenario.v_temp = &v_t;

            ImplicitSolver<Tag_, DataType_, CG, REFLECT> solver(scenario);

            solver.do_preprocessing();
            solver.solve(1);
            std::cout << "After solve:" << std::endl;
            std::cout << "A:" << std::endl;
            std::cout << A << std::endl;

            TEST_CHECK(true);

        }

};
template <typename Tag_, typename DataType_>
class ImplicitSolverRHSAssTest :
    public BaseTest
{
    public:
        ImplicitSolverRHSAssTest(const std::string & type) :
            BaseTest("implicit_solver_rhs_ass_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            Scenario<DataType_, IMPLICIT, REFLECT> scenario(4);
            DenseMatrix<DataType_> h_b(6, 6, DataType_(0));
            DenseMatrix<DataType_> h(4 , 4, DataType_(0));
            DenseMatrix<DataType_> xv_b(6, 6, DataType_(0));
            DenseMatrix<DataType_> xv(4 , 4, DataType_(0));
            DenseMatrix<DataType_> yv_b(6, 6, DataType_(0));
            DenseMatrix<DataType_> yv(4 , 4, DataType_(0));
            DenseMatrix<DataType_> b_b(6, 6, DataType_(0));
            DenseMatrix<DataType_> b(4 , 4, DataType_(0));

            DenseVector<DataType_> rhs(36, DataType_(0));
            BandedMatrix<DataType_> A(36);
            scenario.height = &h;
            scenario.x_veloc = &xv;
            scenario.y_veloc = &yv;
            scenario.bottom = &b;
            scenario.height_bound = &h_b;
            scenario.x_veloc_bound = &xv_b;
            scenario.y_veloc_bound = &yv_b;
            scenario.bottom_bound = &b_b;
            scenario.system_matrix = &A;
            scenario.right_hand_side = &rhs;
            DenseVector<DataType_> u_t(36, DataType_(0));
            DenseVector<DataType_> v_t(36, DataType_(0));

            DataType_ delta_t(1);
            DataType_ delta_x(1);
            DataType_ delta_y(1);

            scenario.delta_t = delta_t;
            scenario.delta_x = delta_x;
            scenario.delta_y = delta_y;

            scenario.u_temp = &u_t;
            scenario.v_temp = &v_t;

            ImplicitSolver<Tag_, DataType_, CG, REFLECT> solver(scenario);
            solver.do_preprocessing();
            solver.solve(1);
            std::cout << "b:" << std::endl;
            std::cout << rhs << std::endl;

            TEST_CHECK(true);

        }
};
template <typename Tag_, typename DataType_>
class ImplicitSolverRuntimeTest :
    public BaseTest
{
    public:
        ImplicitSolverRuntimeTest(const std::string & type) :
            BaseTest("implicit_solver_runtime_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            Scenario<DataType_, IMPLICIT, REFLECT> scenario(4);
            DenseMatrix<DataType_> h_b(6, 6, DataType_(1));
            DenseMatrix<DataType_> h(4 , 4, DataType_(5));
            h[1][1] = 10;
            h[1][2] = 10;
            h[2][1] = 10;
            h[2][2] = 10;
            DenseMatrix<DataType_> xv_b(6, 6, DataType_(0));
            DenseMatrix<DataType_> xv(4 , 4, DataType_(0));
            DenseMatrix<DataType_> yv_b(6, 6, DataType_(0));
            DenseMatrix<DataType_> yv(4 , 4, DataType_(0));
            DenseMatrix<DataType_> b_b(6, 6, DataType_(0));
            DenseMatrix<DataType_> b(4 , 4, DataType_(0));

            DenseVector<DataType_> rhs(36, DataType_(0));
            BandedMatrix<DataType_> A(36);
            scenario.height = &h;
            scenario.x_veloc = &xv;
            scenario.y_veloc = &yv;
            scenario.bottom = &b;
            scenario.height_bound = &h_b;
            scenario.x_veloc_bound = &xv_b;
            scenario.y_veloc_bound = &yv_b;
            scenario.bottom_bound = &b_b;
            scenario.system_matrix = &A;
            scenario.right_hand_side = &rhs;
            DenseVector<DataType_> u_t(36, DataType_(0));
            DenseVector<DataType_> v_t(36, DataType_(0));

            DataType_ delta_t(0.05);
            DataType_ delta_x(1);
            DataType_ delta_y(1);

            scenario.delta_t = delta_t;
            scenario.delta_x = delta_x;
            scenario.delta_y = delta_y;

            scenario.u_temp = &u_t;
            scenario.v_temp = &v_t;

            ImplicitSolver<Tag_, DataType_, CG, REFLECT> solver(scenario);

            solver.do_preprocessing();

            int timesteps = 2;
            for(int timestep = 0; timestep < timesteps; ++timestep)
            {
                std::cout << "Time: " << timestep << std::endl;
                solver.solve(20);
            }
            std::cout << "After solve:" << std::endl;

            std::cout << h_b;
            std::cout << "u: " << xv_b << std::endl;
            std::cout << "v: " << yv_b << std::endl;
            TEST_CHECK(true);

        }
};
//ImplicitSolverCreationTest<tags::CPU, float> implicit_solver_creation_test_float("float");
//ImplicitSolverPreprocessingTest<tags::CPU, float> implicit_solver_preprocessing_test_float("float");
//ImplicitSolverMatrixAssTest<tags::CPU, float> implicit_solver_matrix_test_float("float");
//ImplicitSolverRHSAssTest<tags::CPU, float> implicit_solver_rhs_test_float("float");
//ImplicitSolverCreationTest<tags::CPU, double> implicit_solver_creation_test_double("double");
//ImplicitSolverPreprocessingTest<tags::CPU, double> implicit_solver_preprocessing_test_double("double");
//ImplicitSolverMatrixAssTest<tags::CPU, double> implicit_solver_matrix_test_double("double");
//ImplicitSolverRHSAssTest<tags::CPU, double> implicit_solver_rhs_test_double("double");
ImplicitSolverRuntimeTest<tags::CPU, float> implicit_solver_runtime_test_float("float");
//ImplicitSolverRuntimeTest<tags::CPU, double> implicit_solver_runtime_test_double("double");
/*#ifdef HONEI_SSE
ImplicitSolverRuntimeTest<tags::CPU::SSE, float> sse_implicit_solver_runtime_test_float("SSE float");
ImplicitSolverRuntimeTest<tags::CPU::SSE, double> sse_implicit_solver_runtime_test_double("SSE double");
#endif*/
#ifdef HONEI_CELL
ImplicitSolverRuntimeTest<tags::Cell, float> cell_implicit_solver_runtime_test_float("cell float");
#endif
