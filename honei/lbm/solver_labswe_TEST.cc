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
#include <unittest/unittest.hh>
#include <iostream>

using namespace honei;
using namespace tests;
using namespace std;

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
            unsigned long timesteps(100);

            DenseMatrix<DataType_>h(g_h, g_w, DataType_(0.05));
            DenseMatrix<DataType_>b(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_>u(g_h, g_w, DataType_(0.25));
            DenseMatrix<DataType_>v(g_h, g_w, DataType_(0.));

            SolverLABSWE<Tag_, DataType_,lbm_source_types::SIMPLE, lbm_source_schemes::BASIC, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP_PERIODIC> solver(1.,1.,1., g_h, g_w, &h, &b, &u, &v);

            solver.do_preprocessing();

            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                //std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                solver.solve();
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
