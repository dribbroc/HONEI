/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2009 Markus Geveler <apryde@gmx.de>
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
#include <honei/lbm/solver_lbm_fsi.hh>
#include <honei/lbm/scan_conversion_fsi.hh>
#include <honei/lbm/partial_derivative.hh>
#include <honei/swe/post_processing.hh>
#include <honei/util/unittest.hh>
#include <iostream>
#include <honei/swe/volume.hh>
#include <honei/math/quadrature.hh>
#include <honei/lbm/grid.hh>
#include <honei/lbm/grid_packer.hh>
#include <honei/lbm/scenario_collection.hh>
#include <honei/lbm/solid_emulation_fsi.hh>

using namespace honei;
using namespace tests;
using namespace std;
using namespace output_types;
using namespace lbm::lbm_lattice_types;

//#define SOLVER_VERBOSE
//#define SOLVER_POSTPROCESSING

template <typename Tag_, typename DataType_>
class MinAngleTest :
    public TaggedTest<Tag_>
{
    public:
        MinAngleTest(const std::string & type) :
            TaggedTest<Tag_>("MinAngle test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(10);
            unsigned long g_w(10);

            Grid<D2Q9, DataType_> grid;

            ScenarioCollection::get_scenario(6, g_h, g_w, grid);

            PackedGridData<D2Q9, DataType_>  data;
            PackedSolidData<D2Q9, DataType_>  solids;
            PackedGridInfo<D2Q9> info;

            DenseMatrix<bool> line(g_h, g_w, false);
            DenseMatrix<bool> bound(g_h, g_w, false);
            DenseMatrix<bool> stf(g_h, g_w, false);
            DenseMatrix<bool> sol(g_h, g_w, false);

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);
            GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DataType_>::cuda_pack(info, data);
            GridPackerFSI<D2Q9, NOSLIP, DataType_>::allocate(data, solids);
            GridPackerFSI<D2Q9, NOSLIP, DataType_>::pack(grid, data, solids, line, bound, stf, sol, *grid.obstacles);

            SolverLBMFSI<Tag_, lbm_applications::LABSWE, DataType_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> solver(&info, &data, &solids, grid.d_x, grid.d_y, grid.d_t, grid.tau);

            solver.do_preprocessing();

            solids.current_u = 1.;
            solids.current_v = 0.;
            unsigned long result(MinAngle2D<tags::CPU>::value(solids.current_u, solids.current_v, solids, data));
            std::cout << result << std::endl;;

            solids.current_u = 1.;
            solids.current_v = 1.;
            result = (MinAngle2D<tags::CPU>::value(solids.current_u, solids.current_v, solids, data));
            std::cout << result<< std::endl;;

            solids.current_u = 0.;
            solids.current_v = 1.;
            result = (MinAngle2D<tags::CPU>::value(solids.current_u, solids.current_v, solids, data));
            std::cout << result << std::endl;;

            solids.current_u = -1.;
            solids.current_v = 1.;
            result = (MinAngle2D<tags::CPU>::value(solids.current_u, solids.current_v, solids, data));
            std::cout << result << std::endl;;

            solids.current_u = -1.;
            solids.current_v = 0.;
            result = (MinAngle2D<tags::CPU>::value(solids.current_u, solids.current_v, solids, data));
            std::cout << result << std::endl;;

            solids.current_u = -1.;
            solids.current_v = -1.;
            result = (MinAngle2D<tags::CPU>::value(solids.current_u, solids.current_v, solids, data));
            std::cout << result << std::endl;;

            solids.current_u =  0.;
            solids.current_v = -1.;
            result = (MinAngle2D<tags::CPU>::value(solids.current_u, solids.current_v, solids, data));
            std::cout << result << std::endl;;

            solids.current_u =  1.;
            solids.current_v = -1.;
            result = (MinAngle2D<tags::CPU>::value(solids.current_u, solids.current_v, solids, data));
            std::cout << result << std::endl;;

            solids.current_u =  0.1;
            solids.current_v =  0.1;
            result = (MinAngle2D<tags::CPU>::value(solids.current_u, solids.current_v, solids, data));
            TEST_CHECK_EQUAL(result, 6ul);

            solids.current_u =  0.2;
            solids.current_v =  0.1;
            result = (MinAngle2D<tags::CPU>::value(solids.current_u, solids.current_v, solids, data));
            TEST_CHECK_EQUAL(result, 5ul);

            solids.current_u =  0.1;
            solids.current_v =  3;
            result = (MinAngle2D<tags::CPU>::value(solids.current_u, solids.current_v, solids, data));
            TEST_CHECK_EQUAL(result, 7ul);

            solids.current_u =  -0.1;
            solids.current_v =  -3;
            result = (MinAngle2D<tags::CPU>::value(solids.current_u, solids.current_v, solids, data));
            TEST_CHECK_EQUAL(result, 3ul);

            solids.current_u =  0.1;
            solids.current_v =  -3;
            result = (MinAngle2D<tags::CPU>::value(solids.current_u, solids.current_v, solids, data));
            TEST_CHECK_EQUAL(result, 3ul);

            solids.current_u =  -0.1;
            solids.current_v =  3;
            result = (MinAngle2D<tags::CPU>::value(solids.current_u, solids.current_v, solids, data));
            TEST_CHECK_EQUAL(result, 7ul);

            solids.current_u =  5;
            solids.current_v =  -1;
            result = (MinAngle2D<tags::CPU>::value(solids.current_u, solids.current_v, solids, data));
            TEST_CHECK_EQUAL(result, 5ul);

            solids.current_u = -5;
            solids.current_v =  -1;
            result = (MinAngle2D<tags::CPU>::value(solids.current_u, solids.current_v, solids, data));
            TEST_CHECK_EQUAL(result, 1ul);
        }
};
MinAngleTest<tags::CPU, float> min_angle("float");

template <typename Tag_, typename DataType_>
class SolverLBMFSINonStationaryREACTIONTest :
    public TaggedTest<Tag_>
{
    public:
        SolverLBMFSINonStationaryREACTIONTest(const std::string & type) :
            TaggedTest<Tag_>("solver_lbm_fsi_stationary_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(100);
            unsigned long g_w(100);
            unsigned long timesteps(100);


            Grid<D2Q9, DataType_> grid;

            ScenarioCollection::get_scenario(6, g_h, g_w, grid);

            PackedGridData<D2Q9, DataType_>  data;
            PackedSolidData<D2Q9, DataType_>  solids;
            PackedGridInfo<D2Q9> info;

            DenseMatrix<bool> line(g_h, g_w, false);
            DenseMatrix<bool> bound(g_h, g_w, false);
            DenseMatrix<bool> stf(g_h, g_w, false);
            DenseMatrix<bool> sol(g_h, g_w, false);

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);
            GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DataType_>::cuda_pack(info, data);
            GridPackerFSI<D2Q9, NOSLIP, DataType_>::allocate(data, solids);
            GridPackerFSI<D2Q9, NOSLIP, DataType_>::pack(grid, data, solids, line, bound, stf, sol, *grid.obstacles);

            SolverLBMFSI<Tag_, lbm_applications::LABSWE, DataType_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> solver(&info, &data, &solids, grid.d_x, grid.d_y, grid.d_t, grid.tau);

            //Directly dealing with omega-coordinates
            Line<DataType_, lbm_solid_dims::D2> line_1_0(DataType_(45) * grid.d_x, DataType_(45) * grid.d_y, DataType_(45)* grid.d_x, DataType_(50) * grid.d_y);
            Line<DataType_, lbm_solid_dims::D2> line_2_0(DataType_(45)* grid.d_x, DataType_(50) * grid.d_y, DataType_(50)* grid.d_x, DataType_(50) * grid.d_y);
            Line<DataType_, lbm_solid_dims::D2> line_3_0(DataType_(50)* grid.d_x, DataType_(50) * grid.d_y, DataType_(50)* grid.d_x, DataType_(45) * grid.d_y);
            Line<DataType_, lbm_solid_dims::D2> line_4_0(DataType_(50)* grid.d_x, DataType_(45) * grid.d_y, DataType_(45)* grid.d_x, DataType_(45) * grid.d_y);

            Polygon<DataType_, lbm_solid_dims::D2> tri_0(4);
            tri_0.add_line(line_1_0);
            tri_0.add_line(line_2_0);
            tri_0.add_line(line_3_0);
            tri_0.add_line(line_4_0);
            tri_0.value();

            ScanConversionFSI<tags::CPU>::value(grid, info, data, solids, tri_0, true);

            solids.current_u = 0;
            solids.current_v = 0;
            solver.do_preprocessing();
            std::cout << "Solving: " << grid.description << std::endl;
            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                //Directly dealing with omega-coordinates
                ScanConversionFSI<tags::CPU>::value(grid, info, data, solids, tri_0, true);
                SolidEmulation2D<tags::CPU>::value(tri_0, solids, data, grid);
                std::cout << "(" << solids.current_u << "," << solids.current_v << "), " << MinAngle2D<tags::CPU>::value(solids.current_u, solids.current_v, solids, data) << std::endl;
                solids.current_u = DataType_(0);
                solver.solve(MinAngle2D<tags::CPU>::value(solids.current_u, solids.current_v, solids, data));
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
        }

};
//SolverLBMFSINonStationaryREACTIONTest<tags::CPU, float> solver_test_nstat_float("float");
#ifdef HONEI_CUDA
//SolverLBMFSINonStationaryREACTIONTest<tags::GPU::CUDA, float> solver_test_nstat_cuda_float(" CUDA float");
#endif
