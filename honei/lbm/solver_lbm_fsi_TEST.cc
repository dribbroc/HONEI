/* vim: set number sw=4 sts=4 et nofoldenable : */

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
#include <unittest/unittest.hh>
#include <iostream>
#include <honei/swe/volume.hh>
#include <honei/math/quadrature.hh>
#include <honei/lbm/grid.hh>
#include <honei/lbm/grid_packer.hh>
#include <honei/lbm/scenario_collection.hh>

using namespace honei;
using namespace tests;
using namespace std;
using namespace output_types;
using namespace lbm::lbm_lattice_types;

//#define SOLVER_VERBOSE
//#define SOLVER_POSTPROCESSING

template <typename Tag_, typename DataType_>
class SolverLBMFSITest :
    public TaggedTest<Tag_>
{
    public:
        SolverLBMFSITest(const std::string & type) :
            TaggedTest<Tag_>("solver_lbm_fsi_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(50);
            unsigned long g_w(50);
            unsigned long timesteps(100);


            Grid<D2Q9, DataType_> grid;

            ScenarioCollection::get_scenario(0, g_h, g_w, grid);

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
            std::cout << "Solving: " << grid.description << std::endl;
            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                solver.solve(1ul);
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

};
SolverLBMFSITest<tags::CPU, float> solver_test_float("float");


template <typename Tag_, typename DataType_>
class SolverLBMFSIStationaryTest :
    public TaggedTest<Tag_>
{
    public:
        SolverLBMFSIStationaryTest(const std::string & type) :
            TaggedTest<Tag_>("solver_lbm_fsi_stationary_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(50);
            unsigned long g_w(50);
            unsigned long timesteps(100);


            Grid<D2Q9, DataType_> grid;

            ScenarioCollection::get_scenario(0, g_h, g_w, grid);

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
            Line<DataType_, lbm_solid_dims::D2> line_1_0(DataType_(5) * grid.d_x, DataType_(20) * grid.d_y, DataType_(10)* grid.d_x, DataType_(20) * grid.d_y);
            Line<DataType_, lbm_solid_dims::D2> line_2_0(DataType_(10)* grid.d_x, DataType_(20) * grid.d_y, DataType_(10)* grid.d_x, DataType_(25) * grid.d_y);
            Line<DataType_, lbm_solid_dims::D2> line_3_0(DataType_(10)* grid.d_x, DataType_(25) * grid.d_y, DataType_(5)* grid.d_x, DataType_(25) * grid.d_y);
            Line<DataType_, lbm_solid_dims::D2> line_4_0(DataType_(5)* grid.d_x, DataType_(25) * grid.d_y, DataType_(5)* grid.d_x, DataType_(20) * grid.d_y);

            Polygon<DataType_, lbm_solid_dims::D2> tri_0(4);
            tri_0.add_line(line_1_0);
            tri_0.add_line(line_2_0);
            tri_0.add_line(line_3_0);
            tri_0.add_line(line_4_0);
            tri_0.value();

            ScanConversionFSI<Tag_>::value(grid, info, data, solids, tri_0, true);
            solver.do_preprocessing();
            std::cout << "Solving: " << grid.description << std::endl;
            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                solver.solve(1ul);
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
SolverLBMFSIStationaryTest<tags::CPU, float> solver_test_stat_float("float");

template <typename Tag_, typename DataType_>
class SolverLBMFSINonStationaryTest :
    public TaggedTest<Tag_>
{
    public:
        SolverLBMFSINonStationaryTest(const std::string & type) :
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
            Line<DataType_, lbm_solid_dims::D2> line_1_0(DataType_(5) * grid.d_x, DataType_(20) * grid.d_y, DataType_(10)* grid.d_x, DataType_(20) * grid.d_y);
            Line<DataType_, lbm_solid_dims::D2> line_2_0(DataType_(10)* grid.d_x, DataType_(20) * grid.d_y, DataType_(10)* grid.d_x, DataType_(25) * grid.d_y);
            Line<DataType_, lbm_solid_dims::D2> line_3_0(DataType_(10)* grid.d_x, DataType_(25) * grid.d_y, DataType_(5)* grid.d_x, DataType_(25) * grid.d_y);
            Line<DataType_, lbm_solid_dims::D2> line_4_0(DataType_(5)* grid.d_x, DataType_(25) * grid.d_y, DataType_(5)* grid.d_x, DataType_(20) * grid.d_y);

            Polygon<DataType_, lbm_solid_dims::D2> tri_0(4);
            tri_0.add_line(line_1_0);
            tri_0.add_line(line_2_0);
            tri_0.add_line(line_3_0);
            tri_0.add_line(line_4_0);
            tri_0.value();

            ScanConversionFSI<tags::CPU>::value(grid, info, data, solids, tri_0, true);

            solids.current_u = DataType_(1./2.* grid.d_x);
            solids.current_v = -DataType_(1./2.* grid.d_x);//DataType_(0.);
            solver.do_preprocessing();
            std::cout << "Solving: " << grid.description << std::endl;
            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                //Directly dealing with omega-coordinates
                if (i < 80)
                {
                    Line<DataType_, lbm_solid_dims::D2> line_1_i(DataType_(5.+ i/2.) * grid.d_x, DataType_(20+ i/2.) * grid.d_y, DataType_(10.+ i/2.)* grid.d_x, DataType_(20+ i/2.) * grid.d_y);
                    Line<DataType_, lbm_solid_dims::D2> line_2_i(DataType_(10.+ i/2.)* grid.d_x, DataType_(20+ i/2.) * grid.d_y, DataType_(10.+ i/2.)* grid.d_x, DataType_(25+ i/2.) * grid.d_y);
                    Line<DataType_, lbm_solid_dims::D2> line_3_i(DataType_(10.+ i/2.)* grid.d_x, DataType_(25+ i/2.) * grid.d_y, DataType_(5.+ i/2.)* grid.d_x, DataType_(25+ i/2.) * grid.d_y);
                    Line<DataType_, lbm_solid_dims::D2> line_4_i(DataType_(5.+ i/2.)* grid.d_x, DataType_(25+ i/2.) * grid.d_y, DataType_(5.+ i/2.)* grid.d_x, DataType_(20+ i/2.) * grid.d_y);

                    Polygon<DataType_, lbm_solid_dims::D2> tri_i(4);
                    tri_i.add_line(line_1_i);
                    tri_i.add_line(line_2_i);
                    tri_i.add_line(line_3_i);
                    tri_i.add_line(line_4_i);
                    tri_i.value();
                    ScanConversionFSI<tags::CPU>::value(grid, info, data, solids, tri_i, true);
                }
                else
                {
                    solids.current_u = DataType_(0);
                }
                solver.solve(2ul);
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
SolverLBMFSINonStationaryTest<tags::CPU, float> solver_test_nstat_float("float");
#ifdef HONEI_CUDA
SolverLBMFSINonStationaryTest<tags::GPU::CUDA, float> solver_test_nstat_cuda_float(" CUDA float");
#endif
