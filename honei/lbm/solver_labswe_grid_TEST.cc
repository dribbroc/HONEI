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
#include <honei/lbm/solver_labswe_grid.hh>
#include <honei/lbm/partial_derivative.hh>
#include <honei/swe/post_processing.hh>
#include <unittest/unittest.hh>
#include <iostream>
#include <honei/swe/volume.hh>
#include <honei/math/quadrature.hh>
#include <honei/lbm/grid.hh>
#include <honei/lbm/grid_packer.hh>

using namespace honei;
using namespace tests;
using namespace std;
using namespace output_types;
using namespace lbm::lbm_lattice_types;

//#define SOLVER_VERBOSE
//#define SOLVER_POSTPROCESSING

template <typename Tag_, typename DataType_>
class SolverLABSWEGridTest :
    public TaggedTest<Tag_>
{
    public:
        SolverLABSWEGridTest(const std::string & type) :
            TaggedTest<Tag_>("solver_labswe_grid_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(50);
            unsigned long g_w(50);
            unsigned long timesteps(100);

            DenseMatrix<DataType_> h(g_h, g_w, DataType_(0.05));
            Cylinder<DataType_> c1(h, DataType_(0.02), 25, 25);
            c1.value();

            DenseMatrix<DataType_> u(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> v(g_h, g_w, DataType_(0.));

            DenseMatrix<DataType_> b(g_h, g_w, DataType_(0.));

            Cylinder<DataType_> b1(b, DataType_(0.04), 15, 15);
            b1.value();

            DenseMatrix<DataType_> b_x(PartialDerivative<Tag_, X, CENTRALDIFF>::value(b , DataType_(1)));
            DenseMatrix<DataType_> b_y(PartialDerivative<Tag_, Y, CENTRALDIFF>::value(b , DataType_(1)));

            Grid<D2Q9, DataType_> grid;
            DenseMatrix<bool> obstacles(g_h, g_w, false);
            grid.obstacles = &obstacles;
            grid.h = &h;
            grid.u = &u;
            grid.v = &v;
            grid.b_x = &b_x;
            grid.b_y = &b_y;

            PackedGridData<D2Q9, DataType_>  data;
            PackedGridInfo<D2Q9> info;

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);

            SolverLABSWEGrid<Tag_, DataType_,lbm_source_types::CENTRED, lbm_source_schemes::CENTRALDIFF, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP> solver(&data, &info, 1., 1., 1.);

            solver.do_preprocessing();

            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                solver.solve();
#ifdef SOLVER_POSTPROCESSING
                GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
                PostProcessing<GNUPLOT>::value(*grid.h, 1, g_w, g_h, i);
#endif
            }
#ifdef SOLVER_VERBOSE
            GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
            std::cout << *grid.h << std::endl;
#endif
            TEST_CHECK(true);
        }

};
SolverLABSWEGridTest<tags::CPU, float> solver_test_float("float");
SolverLABSWEGridTest<tags::CPU, double> solver_test_double("double");
SolverLABSWEGridTest<tags::CPU::MultiCore, float> mc_solver_test_float("float");
SolverLABSWEGridTest<tags::CPU::MultiCore, double> mc_solver_test_double("double");
#ifdef HONEI_SSE
SolverLABSWEGridTest<tags::CPU::SSE, float> sse_solver_test_float("float");
SolverLABSWEGridTest<tags::CPU::SSE, double> sse_solver_test_double("double");
SolverLABSWEGridTest<tags::CPU::MultiCore::SSE, float> mcsse_solver_test_float("float");
SolverLABSWEGridTest<tags::CPU::MultiCore::SSE, double> mcsse_solver_test_double("double");
#endif
#ifdef HONEI_CUDA
SolverLABSWEGridTest<tags::GPU::CUDA, float> cuda_solver_test_float("float");
#endif
