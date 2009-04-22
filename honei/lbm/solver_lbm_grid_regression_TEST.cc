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
#include <honei/swe/post_processing.hh>
#include <honei/swe/volume.hh>
#include <unittest/unittest.hh>
#include <iostream>
#include <honei/swe/volume.hh>
#include <honei/math/quadrature.hh>
#include <honei/lbm/grid.hh>
#include <honei/lbm/grid_packer.hh>
#include <honei/la/norm.hh>
#include <honei/la/difference.hh>
#include <honei/lbm/scenario_collection.hh>

using namespace honei;
using namespace tests;
using namespace std;
using namespace output_types;
using namespace lbm::lbm_lattice_types;

//#define SOLVER_VERBOSE 1

template <typename Tag_, typename DataType_>
class SolverLBMGridRegressionTest :
    public TaggedTest<Tag_>
{
    public:
        SolverLBMGridRegressionTest(const std::string & type) :
            TaggedTest<Tag_>("solver_lbm_grid_regression_test<" + type + ">")
    {
    }

        virtual void run() const
        {
            for (unsigned long scen(0) ; scen < ScenarioCollection::get_stable_scenario_count() ; ++scen)
            {
                unsigned long g_h(50);
                unsigned long g_w(50);
                unsigned long timesteps(200);

                Grid<D2Q9, DataType_> grid;
                ScenarioCollection::get_scenario(scen, g_h, g_w, grid);

                PackedGridData<D2Q9, DataType_>  data;
                PackedGridInfo<D2Q9> info;

                GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);

                SolverLBMGrid<Tag_, lbm_applications::LABSWE, DataType_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> solver(&info, &data, grid.d_x, grid.d_y, grid.d_t, grid.tau);
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

#ifdef SOLVER_VERBOSE
                std::cout << *grid.h << std::endl;
#endif
                //Standard solver using tags::CPU:
                unsigned long g_h_standard(50);
                unsigned long g_w_standard(50);
                unsigned long timesteps_standard(200);

                Grid<D2Q9, DataType_> grid_standard;
                ScenarioCollection::get_scenario(scen, g_h_standard, g_w_standard, grid_standard);

                PackedGridData<D2Q9, DataType_>  data_standard;
                PackedGridInfo<D2Q9> info_standard;

                GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid_standard, info_standard, data_standard);

                SolverLBMGrid<tags::CPU, lbm_applications::LABSWE, DataType_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> solver_standard(&info_standard, &data_standard, grid_standard.d_x, grid_standard.d_y, grid_standard.d_t, grid_standard.tau);

                solver_standard.do_preprocessing();

                for(unsigned long i(0); i < timesteps_standard; ++i)
                {
#ifdef SOLVER_VERBOSE
                    std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                    solver_standard.solve();
#ifdef SOLVER_POSTPROCESSING
                    solver.do_postprocessing();
                    GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid_standard, info_standard, data_standard);
                    PostProcessing<GNUPLOT>::value(*grid_standard.h, 1, g_w_standard, g_h_standard, i);
#endif
                }
                solver.do_postprocessing();
                GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid_standard, info_standard, data_standard);

#ifdef SOLVER_VERBOSE
                std::cout << *grid_standard.h << std::endl;
#endif

                std::cout << grid.description <<": ";

                TEST_CHECK_EQUAL(g_h, g_h_standard);
                TEST_CHECK_EQUAL(g_w, g_w_standard);

                //Compare CPU results of both solvers:
                for(unsigned long i(0) ; i < g_h ; ++i)
                {
                    for(unsigned long j(0) ; j < g_w ; ++j)
                    {
#ifdef SOLVER_VERBOSE
                        std::cout << "(" << i << " , " << j << ")" << std::endl;
                        std::cout << (*grid.h)(i , j) << " " << (*grid_standard.h)(i , j) << std::endl;
#endif
                        TEST_CHECK_EQUAL_WITHIN_EPS((*grid.h)(i , j) , (*grid_standard.h)(i , j), std::numeric_limits<DataType_>::epsilon() * 2e2);
                    }
                }

                //Save matrices to vectors, compute norm:
                DenseVector<double> result_grid(g_h*g_w);
                DenseVector<double> result_standard(g_h*g_w);

                unsigned long inner(0);
                for(unsigned long i(0) ; i < g_h ; ++i)
                {
                    for(unsigned long j(0) ; j < g_w ; ++j)
                    {
                        result_grid[inner] = (*grid.h)(i , j);
                        result_standard[inner] = (*grid_standard.h)(i , j);
                        ++inner;
                    }
                }

                Difference<tags::CPU>::value(result_grid, result_standard);
                double l2 = Norm<vnt_l_two, false, tags::CPU>::value(result_grid);
                TEST_CHECK_EQUAL_WITHIN_EPS(l2, DataType_(0.), std::numeric_limits<DataType_>::epsilon() * 2);

                std::cout << "L2 norm " << l2 << std::endl;
            }
        }
};
SolverLBMGridRegressionTest<tags::CPU::MultiCore, float> mc_solver_test_float("float");
SolverLBMGridRegressionTest<tags::CPU::MultiCore, double> mc_solver_test_double("double");
#ifdef HONEI_SSE
SolverLBMGridRegressionTest<tags::CPU::SSE, float> sse_solver_test_float("float");
SolverLBMGridRegressionTest<tags::CPU::SSE, double> sse_solver_test_double("double");
SolverLBMGridRegressionTest<tags::CPU::MultiCore::SSE, float> mcsse_solver_test_float("float");
SolverLBMGridRegressionTest<tags::CPU::MultiCore::SSE, double> mcsse_solver_test_double("double");
#endif
#ifdef HONEI_CUDA
SolverLBMGridRegressionTest<tags::GPU::CUDA, float> cuda_solver_test_float("float");
#endif
#ifdef HONEI_CELL
SolverLBMGridRegressionTest<tags::Cell, float> cell_solver_test_float("float");
#endif

