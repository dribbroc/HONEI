/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the SWE C++ library. LiSWE is free software;
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
#include <honei/lbm/solver_labswe.hh>
#include <honei/swe/post_processing.hh>
#include <unittest/unittest.hh>
#include <iostream>

using namespace honei;
using namespace tests;
using namespace std;
using namespace output_types;

template <typename Tag_, typename DataType_>
class SolverLABSWETest :
    public TaggedTest<Tag_>
{
    public:
        SolverLABSWETest(const std::string & type) :
            TaggedTest<Tag_>("solver_labswe_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(20);
            unsigned long g_w(20);
            unsigned long timesteps(5);

            DenseMatrix<DataType_> h(g_h, g_w, DataType_(0.05));
            DenseMatrix<DataType_> b(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> u(g_h, g_w, DataType_(0.25));
            DenseMatrix<DataType_> v(g_h, g_w, DataType_(0.));

            //All needed distribution functions:

            DenseMatrix<DataType_> d_0(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_1(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_2(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_3(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_4(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_5(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_6(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_7(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_8(g_h, g_w, DataType_(0.));

            DenseMatrix<DataType_> e_d_0(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_1(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_2(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_3(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_4(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_5(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_6(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_7(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> e_d_8(g_h, g_w, DataType_(0.));

            DenseMatrix<DataType_> t_d_0(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_1(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_2(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_3(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_4(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_5(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_6(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_7(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> t_d_8(g_h, g_w, DataType_(0.));

            //All needed vectors:
            DenseVector<DataType_> v_x(9, DataType_(0));
            DenseVector<DataType_> v_y(9, DataType_(0));

            //Other matrices needed by solver:

            DenseMatrix<DataType_> s_x(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> s_y(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_x(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> d_y(g_h, g_w, DataType_(0.));

            SolverLABSWE<Tag_, DataType_,lbm_source_types::SIMPLE, lbm_source_schemes::BASIC, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP_PERIODIC> solver(1.,1.,1., g_h, g_w, &h, &b, &u, &v);

            solver.set_distribution(&d_0, &d_1, &d_2, &d_3, &d_4, &d_5, &d_6, &d_7, &d_8);
            solver.set_eq_distribution(&e_d_0, &e_d_1, &e_d_2, &e_d_3, &e_d_4, &e_d_5, &e_d_6, &e_d_7, &e_d_8);
            solver.set_temp_distribution(&t_d_0, &t_d_1, &t_d_2, &t_d_3, &t_d_4, &t_d_5, &t_d_6, &t_d_7, &t_d_8);
            solver.set_vectors(&v_x, &v_y);
            solver.set_source(&s_x, &s_y);
            solver.set_slopes(&d_x, &d_y);
            solver.do_preprocessing();

            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                //std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                solver.solve();
                PostProcessing<GNUPLOT>::value(h, 100, g_w, g_h, i);
            }
            std::cout << h << std::endl;
            TEST_CHECK(true);
        }

};
/*SolverLABSWETest<tags::CPU, float> solver_test_float("float");
SolverLABSWETest<tags::CPU, double> solver_test_double("double");
SolverLABSWETest<tags::CPU::MultiCore, float> solver_test_float_mc("float");
SolverLABSWETest<tags::CPU::MultiCore, double> solver_test_double_mc("double");*/
#ifdef HONEI_SSE
//SolverLABSWETest<tags::CPU::SSE, float> solver_test_float_sse("float");
SolverLABSWETest<tags::CPU::SSE, double> solver_test_double_sse("double");
//SolverLABSWETest<tags::CPU::MultiCore::SSE, float> solver_test_float_mc_sse("float");
//SolverLABSWETest<tags::CPU::MultiCore::SSE, double> solver_test_double_mc_sse("double");
#endif
#ifdef HONEI_CELL
SolverLABSWETest<tags::Cell, float> solver_test_float_cell("float");
SolverLABSWETest<tags::Cell, double> solver_test_double_cell("double");
#endif
