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
#include <honei/lbm/solver_lbm_grid.hh>
#include <honei/lbm/partial_derivative.hh>
#include <honei/swe/post_processing.hh>
#include <honei/util/unittest.hh>
#include <iostream>
#include <honei/swe/volume.hh>
#include <honei/math/quadrature.hh>
#include <honei/lbm/grid.hh>
#include <honei/lbm/grid_packer.hh>
#include <honei/lbm/scenario_collection.hh>

#ifdef DEBUG
#define SOLVER_VERBOSE
#endif

using namespace honei;
using namespace tests;
using namespace std;
using namespace output_types;
using namespace lbm::lbm_lattice_types;

//#define SOLVER_VERBOSE
//#define SOLVER_POSTPROCESSING

template <typename Tag_, typename DataType_>
class SolverLBMGridTest :
    public TaggedTest<Tag_>
{
    public:
        SolverLBMGridTest(const std::string & type) :
            TaggedTest<Tag_>("solver_lbm_grid_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(50);
            unsigned long g_w(50);
            unsigned long timesteps(100);


            for (unsigned long i(0) ; i < ScenarioCollection::get_stable_scenario_count() ; ++i)
            {
                Grid<D2Q9, DataType_> grid;

                ScenarioCollection::get_scenario(i, g_h, g_w, grid);

                PackedGridData<D2Q9, DataType_>  data;
                PackedGridInfo<D2Q9> info;

                GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);

                SolverLBMGrid<Tag_, lbm_applications::LABSWE, DataType_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> solver(&info, &data, grid.d_x, grid.d_y, grid.d_t, grid.tau);

                solver.do_preprocessing();
                std::cout << "Solving: " << grid.description << std::endl;
                for(unsigned long i(0); i < timesteps; ++i)
                {
#ifdef SOLVER_VERBOSE
                    std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                    solver.solve();
#ifdef SOLVER_POSTPROCESSING
                    solver.do_postprocessing();
                    GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
                    PostProcessing<GNUPLOT>::value(*grid.h, 1, g_w, g_h, i);
#endif
                }
                solver.do_postprocessing();
                GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
#ifdef SOLVER_VERBOSE
                std::cout << *grid.h << std::endl;
#endif
                for (unsigned long i(0) ; i < (*grid.h).rows() ; ++i)
                    for(unsigned long j(0) ; j < (*grid.h).columns() ; ++j)
                        TEST_CHECK_EQUAL_WITHIN_EPS((*grid.h)( i , j), DataType_(0.02), DataType_(0.1));
            }
        }

};

SolverLBMGridTest<tags::CPU, float> solver_test_float("float");
SolverLBMGridTest<tags::CPU, double> solver_test_double("double");
SolverLBMGridTest<tags::CPU::MultiCore, float> mc_solver_test_float("float");
SolverLBMGridTest<tags::CPU::MultiCore, double> mc_solver_test_double("double");
#ifdef HONEI_SSE
SolverLBMGridTest<tags::CPU::SSE, float> sse_solver_test_float("float");
SolverLBMGridTest<tags::CPU::SSE, double> sse_solver_test_double("double");
SolverLBMGridTest<tags::CPU::MultiCore::SSE, float> mcsse_solver_test_float("float");
SolverLBMGridTest<tags::CPU::MultiCore::SSE, double> mcsse_solver_test_double("double");
#endif
#ifdef HONEI_CUDA
SolverLBMGridTest<tags::GPU::CUDA, float> cuda_solver_test_float("float");
SolverLBMGridTest<tags::GPU::MultiCore::CUDA, float> mc_cuda_solver_test_float("float");
#endif
#ifdef HONEI_CELL
SolverLBMGridTest<tags::Cell, float> cell_solver_test_float("float");
#endif


template <typename Tag_, typename DataType_>
class SolverLBMGridMassConservationTest :
    public TaggedTest<Tag_>
{
    public:
        SolverLBMGridMassConservationTest(const std::string & type) :
            TaggedTest<Tag_>("solver_lbm_grid_mass_cons_test<" + type + ">")
    {
    }

        virtual void run() const
        {
            unsigned long g_h(50);
            unsigned long g_w(50);
            unsigned long timesteps(100);

            DenseMatrix<DataType_> h(g_h, g_w, DataType_(0.05));
            Cuboid<DataType_> a(h, 5, 5, DataType_(0.02), 25, 25);
            a.value();

            DenseMatrix<DataType_> u(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> v(g_h, g_w, DataType_(0.));

            DenseMatrix<DataType_> b(g_h, g_w, DataType_(0.));

            Grid<D2Q9, DataType_> grid;
            DenseMatrix<bool> obstacles(g_h, g_w, false);
            grid.obstacles = &obstacles;
            grid.h = &h;
            grid.u = &u;
            grid.v = &v;
            grid.b = &b;

            PackedGridData<D2Q9, DataType_>  data;
            PackedGridInfo<D2Q9> info;

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);

            SolverLBMGrid<Tag_, lbm_applications::LABSWE, DataType_, lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> solver(&info, &data, 1., 1., 1., 1.5);

            solver.do_preprocessing();

            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                solver.solve();
#ifdef SOLVER_POSTPROCESSING
                solver.do_postprocessing();
                GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
                PostProcessing<GNUPLOT>::value(*grid.h, 1, g_w, g_h, i);
#endif
            }
            solver.do_postprocessing();
            GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
            DataType_ ana_vol(5. * 5. * 0.02 + g_w * g_h * 0.05);
            std::cout << "Analytical Vol.: " << ana_vol << " ";
            DataType_ vol = GaussianQuadrature2D<tags::CPU, tags::Trapezoid>::value(h, DataType_(0), DataType_(g_w), DataType_(1.), DataType_(1.));
            std::cout << "Vol.: " << vol << std::endl;
            TEST_CHECK_EQUAL_WITHIN_EPS(vol, ana_vol, 0.1);
        }

};
SolverLBMGridMassConservationTest<tags::CPU, float> solver_grid_mc_test_float("float");
SolverLBMGridMassConservationTest<tags::CPU, double> solver_grid_mc_test_double("double");
SolverLBMGridMassConservationTest<tags::CPU::MultiCore, float> mc_solver_grid_mc_test_float("float");
SolverLBMGridMassConservationTest<tags::CPU::MultiCore, double> mc_solver_grid_mc_test_double("double");
#ifdef HONEI_SSE
SolverLBMGridMassConservationTest<tags::CPU::SSE, float> sse_solver_grid_mc_test_float("float");
SolverLBMGridMassConservationTest<tags::CPU::SSE, double> sse_solver_grid_mc_test_double("double");
SolverLBMGridMassConservationTest<tags::CPU::MultiCore::SSE, float> mcsse_solver_grid_mc_test_float("float");
SolverLBMGridMassConservationTest<tags::CPU::MultiCore::SSE, double> mcsse_solver_grid_mc_test_double("double");
#endif
#ifdef HONEI_CUDA
SolverLBMGridMassConservationTest<tags::GPU::CUDA, float> cuda_solver_grid_mc_test_float("float");
SolverLBMGridMassConservationTest<tags::GPU::MultiCore::CUDA, float> mc_cuda_solver_grid_mc_test_float("float");
#endif
#ifdef HONEI_Cell
SolverLBMGridMassConservationTest<tags::Cell, float> cell_solver_grid_mc_test_float("float");
#endif


template <typename Tag_, typename DataType_>
class SimpleSolverLBMGridMassConservationTest :
    public TaggedTest<Tag_>
{
    public:
        SimpleSolverLBMGridMassConservationTest(const std::string & type) :
            TaggedTest<Tag_>("simple_solver_lbm_grid_mass_cons_test<" + type + ">")
    {
    }

        virtual void run() const
        {
            unsigned long g_h(50);
            unsigned long g_w(50);
            unsigned long timesteps(100);

            DenseMatrix<DataType_> h(g_h, g_w, DataType_(0.05));
            Cuboid<DataType_> a(h, 5, 5, DataType_(0.02), 25, 25);
            a.value();

            DenseMatrix<DataType_> u(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> v(g_h, g_w, DataType_(0.));

            DenseMatrix<DataType_> b(g_h, g_w, DataType_(0.));

            Grid<D2Q9, DataType_> grid;
            DenseMatrix<bool> obstacles(g_h, g_w, false);
            grid.obstacles = &obstacles;
            grid.h = &h;
            grid.u = &u;
            grid.v = &v;
            grid.b = &b;

            PackedGridData<D2Q9, DataType_>  data;
            PackedGridInfo<D2Q9> info;

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);

            SolverLBMGrid<Tag_, lbm_applications::LABSWE, DataType_, lbm_force::NONE, lbm_source_schemes::NONE, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::WET> solver(&info, &data, 1., 1., 1., 1.5);

            solver.do_preprocessing();

            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                solver.solve();
#ifdef SOLVER_POSTPROCESSING
                solver.do_postprocessing();
                GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
                PostProcessing<GNUPLOT>::value(*grid.h, 1, g_w, g_h, i);
#endif
            }
            solver.do_postprocessing();
            GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
            DataType_ ana_vol(5. * 5. * 0.02 + g_w * g_h * 0.05);
            std::cout << "Analytical Vol.: " << ana_vol << " ";
            DataType_ vol = GaussianQuadrature2D<tags::CPU, tags::Trapezoid>::value(h, DataType_(0), DataType_(g_w), DataType_(1.), DataType_(1.));
            std::cout << "Vol.: " << vol << std::endl;
            TEST_CHECK_EQUAL_WITHIN_EPS(vol, ana_vol, 0.1);
        }

};
SimpleSolverLBMGridMassConservationTest<tags::CPU, float> simple_solver_grid_mc_test_float("float");
SimpleSolverLBMGridMassConservationTest<tags::CPU, double> simple_solver_grid_mc_test_double("double");
SimpleSolverLBMGridMassConservationTest<tags::CPU::MultiCore, float> mc_simple_solver_grid_mc_test_float("float");
SimpleSolverLBMGridMassConservationTest<tags::CPU::MultiCore, double> mc_simple_solver_grid_mc_test_double("double");
#ifdef HONEI_SSE
SimpleSolverLBMGridMassConservationTest<tags::CPU::SSE, float> sse_simple_solver_grid_mc_test_float("float");
SimpleSolverLBMGridMassConservationTest<tags::CPU::SSE, double> sse_simple_solver_grid_mc_test_double("double");
SimpleSolverLBMGridMassConservationTest<tags::CPU::MultiCore::SSE, float> mcsse_simple_solver_grid_mc_test_float("float");
SimpleSolverLBMGridMassConservationTest<tags::CPU::MultiCore::SSE, double> mcsse_simple_solver_grid_mc_test_double("double");
#endif
#ifdef HONEI_CUDA
SimpleSolverLBMGridMassConservationTest<tags::GPU::CUDA, float> cuda_simple_solver_grid_mc_test_float("float");
SimpleSolverLBMGridMassConservationTest<tags::GPU::MultiCore::CUDA, float> mc_cuda_simple_solver_grid_mc_test_float("float");
#endif
#ifdef HONEI_Cell
SimpleSolverLBMGridMassConservationTest<tags::Cell, float> cell_simple_solver_grid_mc_test_float("float");
#endif
