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
#include <honei/lbm/solver_lbm_grid.hh>
#include <honei/swe/post_processing.hh>
#include <unittest/unittest.hh>
#include <iostream>
#include <honei/lbm/grid.hh>
#include <honei/lbm/grid_packer.hh>
#include <honei/lbm/scenario_collection.hh>
#include <honei/lbm/solid.hh>

using namespace honei;
using namespace lbm;
using namespace tests;
using namespace std;
using namespace output_types;
using namespace lbm_lattice_types;

//#define SOLVER_VERBOSE
//#define SOLVER_POSTPROCESSING

template <typename Tag_, typename DataType_>
class FSISolverLBMGridTest :
    public TaggedTest<Tag_>
{
    public:
        FSISolverLBMGridTest(const std::string & type) :
            TaggedTest<Tag_>("fsi_lbm_grid_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(50);
            unsigned long g_w(50);
            unsigned long timesteps(100);

            Grid<D2Q9, DataType_> grid;

            ///Use STEADY STATE values:
            ScenarioCollection::get_scenario(6, g_h, g_w, grid);

            PackedGridData<D2Q9, DataType_>  data;
            PackedGridInfo<D2Q9> info;
            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);

            SolverLBMGrid<Tag_, lbm_applications::LABSWE, DataType_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> solver(&info, &data, grid.d_x, grid.d_y, grid.d_t, grid.tau);


            DenseMatrix<bool> obstacles_new(g_h, g_w, false);
            DenseMatrix<bool> fts(g_h, g_w, false);
            DenseMatrix<bool> stf(g_h, g_w, false);
            std::cout << "Solving: " << grid.description << std::endl;
            for(unsigned long i(0); i < timesteps; ++i)
            {
                for(unsigned long j(0); j < obstacles_new.rows() ; ++j)
                {
                    for(unsigned long k(0) ; k < obstacles_new.columns() ; ++k)
                    {
                        obstacles_new[j][k] = false;
                        fts[j][k] = false;
                        stf[j][k] = false;
                    }
                }
                GridPacker<D2Q9, NOSLIP, DataType_>::deflate(grid, data, data.h, grid.h);
                GridPacker<D2Q9, NOSLIP, DataType_>::deflate(grid, data, data.u, grid.u);
                GridPacker<D2Q9, NOSLIP, DataType_>::deflate(grid, data, data.v, grid.v);

                //Directly dealing with omega-coordinates
                Line<DataType_, lbm_solid_dims::D2> line_1(DataType_(5 + i) * grid.d_x, DataType_(20) * grid.d_y, DataType_(10 + i)* grid.d_x, DataType_(20) * grid.d_y);
                Line<DataType_, lbm_solid_dims::D2> line_2(DataType_(10 + i)* grid.d_x, DataType_(20) * grid.d_y, DataType_(10 + i)* grid.d_x, DataType_(25) * grid.d_y);
                Line<DataType_, lbm_solid_dims::D2> line_3(DataType_(10 + i)* grid.d_x, DataType_(25) * grid.d_y, DataType_(5 + i)* grid.d_x, DataType_(25) * grid.d_y);
                Line<DataType_, lbm_solid_dims::D2> line_4(DataType_(5 + i)* grid.d_x, DataType_(25) * grid.d_y, DataType_(5 + i)* grid.d_x, DataType_(20) * grid.d_y);

                Polygon<DataType_, lbm_solid_dims::D2> tri(4);
                tri.add_line(line_1);
                tri.add_line(line_2);
                tri.add_line(line_3);
                tri.add_line(line_4);
                tri.value();

                ScanConversion<Tag_>::value(tri, obstacles_new, grid.d_x, grid.d_y, true);
#ifdef SOLVER_VERBOSE
                if(i == 0)
                    std::cout << obstacles_new << std::endl;
#endif
                FluidToSolidCells<Tag_>::value(*grid.obstacles, obstacles_new, fts);
                FTSExtrapolation<Tag_, lbm_solid_extrapolation_methods::SIMPLE>::value(obstacles_new, fts, *grid.h, *grid.u, *grid.v, DataType_(1), DataType_(0), DataType_(0.001), DataType_(grid.d_x));

                SolidToFluidCells<Tag_>::value(*grid.obstacles, obstacles_new, stf);
                STFExtrapolation<Tag_, lbm_solid_extrapolation_methods::SIMPLE>::value(obstacles_new, stf, *grid.h, *grid.u, *grid.v, DataType_(1), DataType_(0), DataType_(0.001), DataType_(grid.d_x));

                info.destroy();
                data.destroy();
                GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);
                solver.do_preprocessing();
                solver.solve();
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
#ifdef SOLVER_POSTPROCESSING
                solver.do_postprocessing();
                PostProcessing<GNUPLOT>::value(*grid.h, 1, g_w, g_h, i);
#endif
            }
#ifdef SOLVER_VERBOSE
            std::cout << *grid.h << std::endl;
#endif
        }
};

FSISolverLBMGridTest<tags::CPU, float> fs_solver_test_float("float");
