/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2011 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
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
#include <honei/util/time_stamp.hh>

#include <honei/woolb3/grid3.hh>
#include <honei/woolb3/packed_grid3.hh>
#include <honei/woolb3/solver_lbm3.hh>

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
class SolverLBM3Test :
    public TaggedTest<Tag_>
{
    public:
        SolverLBM3Test(const std::string & type) :
            TaggedTest<Tag_>("solver_lbm3_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(128);
            unsigned long g_w(128);
            unsigned long timesteps(250);


            Grid<D2Q9, DataType_> grid;

            ScenarioCollection::get_scenario(1, g_h, g_w, grid);

            PackedGridData<D2Q9, DataType_>  data;
            PackedGridInfo<D2Q9> info;

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);

            SolverLBMGrid<Tag_, lbm_applications::LABSWE, DataType_,lbm_force::NONE, lbm_source_schemes::NONE, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::WET> solver(&info, &data, grid.d_x, grid.d_y, grid.d_t, grid.tau);

            solver.do_preprocessing();


            grid.destroy();
            ScenarioCollection::get_scenario(1, g_h, g_w, grid);
            Grid3<DataType_, 9> grid3(*grid.obstacles, *grid.h, *grid.b, *grid.u, *grid.v);
            PackedGrid3<DataType_, 9> pgrid3(grid3);
            SolverLBM3<Tag_, DataType_, 9> solver3(grid3, pgrid3, grid.d_x, grid.d_y, grid.d_t, grid.tau);
            solver3.do_preprocessing();

            /*TEST_CHECK_EQUAL( *pgrid3.f_eq[0],*(data.f_eq_0));
            TEST_CHECK_EQUAL( *pgrid3.f_eq[1],*(data.f_eq_1));
            TEST_CHECK_EQUAL( *pgrid3.f_eq[2],*(data.f_eq_2));
            TEST_CHECK_EQUAL( *pgrid3.f_eq[3],*(data.f_eq_3));
            TEST_CHECK_EQUAL( *pgrid3.f_eq[4],*(data.f_eq_4));
            TEST_CHECK_EQUAL( *pgrid3.f_eq[5],*(data.f_eq_5));
            TEST_CHECK_EQUAL( *pgrid3.f_eq[6],*(data.f_eq_6));
            TEST_CHECK_EQUAL( *pgrid3.f_eq[7],*(data.f_eq_7));
            TEST_CHECK_EQUAL( *pgrid3.f_eq[8],*(data.f_eq_8));

            TEST_CHECK_EQUAL( *pgrid3.f_temp[0],*(data.f_temp_0));
            TEST_CHECK_EQUAL( *pgrid3.f_temp[1],*(data.f_temp_1));
            TEST_CHECK_EQUAL( *pgrid3.f_temp[2],*(data.f_temp_2));
            TEST_CHECK_EQUAL( *pgrid3.f_temp[3],*(data.f_temp_3));
            TEST_CHECK_EQUAL( *pgrid3.f_temp[4],*(data.f_temp_4));
            TEST_CHECK_EQUAL( *pgrid3.f_temp[5],*(data.f_temp_5));
            TEST_CHECK_EQUAL( *pgrid3.f_temp[6],*(data.f_temp_6));
            TEST_CHECK_EQUAL( *pgrid3.f_temp[7],*(data.f_temp_7));
            TEST_CHECK_EQUAL( *pgrid3.f_temp[8],*(data.f_temp_8));*/


            TimeStamp at, bt;
            at.take();
            for (unsigned long i(0) ; i < timesteps ; ++i)
            {
                solver.solve();
            }
            bt.take();
            std::cout<<"LBMGrid TOE: "<<bt.total()-at.total()<<std::endl;

            at.take();
            for (unsigned long i(0) ; i < timesteps ; ++i)
            {
                solver3.solve();
            }
            bt.take();
            std::cout<<"LBM3 TOE: "<<bt.total()-at.total()<<std::endl;

            /*TEST_CHECK_EQUAL( *pgrid3.f_temp[0],*(data.f_temp_0));
            TEST_CHECK_EQUAL( *pgrid3.f_temp[1],*(data.f_temp_1));
            TEST_CHECK_EQUAL( *pgrid3.f_temp[2],*(data.f_temp_2));
            TEST_CHECK_EQUAL( *pgrid3.f_temp[3],*(data.f_temp_3));
            TEST_CHECK_EQUAL( *pgrid3.f_temp[4],*(data.f_temp_4));
            TEST_CHECK_EQUAL( *pgrid3.f_temp[5],*(data.f_temp_5));
            TEST_CHECK_EQUAL( *pgrid3.f_temp[6],*(data.f_temp_6));
            TEST_CHECK_EQUAL( *pgrid3.f_temp[7],*(data.f_temp_7));
            TEST_CHECK_EQUAL( *pgrid3.f_temp[8],*(data.f_temp_8));

            TEST_CHECK_EQUAL( *pgrid3.f_eq[0],*(data.f_eq_0));
            TEST_CHECK_EQUAL( *pgrid3.f_eq[1],*(data.f_eq_1));
            TEST_CHECK_EQUAL( *pgrid3.f_eq[2],*(data.f_eq_2));
            TEST_CHECK_EQUAL( *pgrid3.f_eq[3],*(data.f_eq_3));
            TEST_CHECK_EQUAL( *pgrid3.f_eq[4],*(data.f_eq_4));
            TEST_CHECK_EQUAL( *pgrid3.f_eq[5],*(data.f_eq_5));
            TEST_CHECK_EQUAL( *pgrid3.f_eq[6],*(data.f_eq_6));
            TEST_CHECK_EQUAL( *pgrid3.f_eq[7],*(data.f_eq_7));
            TEST_CHECK_EQUAL( *pgrid3.f_eq[8],*(data.f_eq_8));

            TEST_CHECK_EQUAL( *pgrid3.h,*(data.h));*/
            //TEST_CHECK_EQUAL( *pgrid3.u,*(data.u));
            //TEST_CHECK_EQUAL( *pgrid3.v,*(data.v));


            solver.do_postprocessing();
            GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);

            DenseMatrix<DataType_> h_3(grid.h->rows(), grid.h->columns(), DataType_(0));
            grid3.fill_h(h_3, *pgrid3.h);

            for (unsigned long row(0) ; row < h_3.rows() ; ++row)
                for (unsigned long col(0) ; col < h_3.columns() ; ++col)
                    TEST_CHECK_EQUAL_WITHIN_EPS(h_3(row, col), (*grid.h)(row, col), 1e-6);


            info.destroy();
            data.destroy();
            grid.destroy();
        }

};
SolverLBM3Test<tags::CPU, float> solver_test_float("float");

template <typename Tag_, typename DataType_>
class MultiSolverLBM3Test :
    public TaggedTest<Tag_>
{
    public:
        MultiSolverLBM3Test(const std::string & type) :
            TaggedTest<Tag_>("multi_solver_lbm3_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(128);
            unsigned long g_w(128);
            unsigned long timesteps(250);


            Grid<D2Q9, DataType_> grid;
            ScenarioCollection::get_scenario(1, g_h, g_w, grid);
            DenseMatrix<DataType_> h_p(grid.h->rows(), grid.h->columns(), 0);
            DenseMatrix<DataType_> h_s(grid.h->rows(), grid.h->columns(), 0);

            // serial solver
            Grid3<DataType_, 9> grid_s(*grid.obstacles, *grid.h, *grid.b, *grid.u, *grid.v);
            PackedGrid3<DataType_, 9> pgrid_s(grid_s);
            SolverLBM3<Tag_, DataType_, 9> solver_s(grid_s, pgrid_s, grid.d_x, grid.d_y, grid.d_t, grid.tau);
            solver_s.do_preprocessing();

            for (unsigned long i(0) ; i < timesteps ; ++i)
            {
                solver_s.solve();
            }


            // parallel solver
            //  TODO use stl vector for multiple process count
            Grid3<DataType_, 9> grid_p1(*grid.obstacles, *grid.h, *grid.b, *grid.u, *grid.v, 0, 2);
            PackedGrid3<DataType_, 9> pgrid_p1(grid_p1);
            SolverLBM3<Tag_, DataType_, 9> solver_p1(grid_p1, pgrid_p1, grid.d_x, grid.d_y, grid.d_t, grid.tau);
            solver_p1.do_preprocessing();

            Grid3<DataType_, 9> grid_p2(*grid.obstacles, *grid.h, *grid.b, *grid.u, *grid.v, 1, 2);
            PackedGrid3<DataType_, 9> pgrid_p2(grid_p2);
            SolverLBM3<Tag_, DataType_, 9> solver_p2(grid_p2, pgrid_p2, grid.d_x, grid.d_y, grid.d_t, grid.tau);
            solver_p2.do_preprocessing();


            std::list<SyncData<DataType_> > p1_data = pgrid_p1.export_synch_data();
            std::list<SyncData<DataType_> > p2_data = pgrid_p2.export_synch_data();
            pgrid_p1.import_synch_data(p2_data);
            pgrid_p2.import_synch_data(p1_data);


            for (unsigned long i(0) ; i < timesteps ; ++i)
            {
                solver_p1.solve_outer();
                solver_p2.solve_outer();

                std::list<SyncData<DataType_> > p1_data = pgrid_p1.export_synch_data();
                std::list<SyncData<DataType_> > p2_data = pgrid_p2.export_synch_data();

                solver_p1.solve_inner();
                solver_p2.solve_inner();

                pgrid_p1.import_synch_data(p2_data);
                pgrid_p2.import_synch_data(p1_data);
            }

            grid_s.fill_h(h_s, *pgrid_s.h);
            grid_p1.fill_h(h_p, *pgrid_p1.h);
            grid_p2.fill_h(h_p, *pgrid_p2.h);

            for (unsigned long row(0) ; row < h_s.rows() ; ++row)
                for (unsigned long col(0) ; col < h_s.columns() ; ++col)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(h_p(row, col), h_s(row, col), 1e-6);
                }


            grid.destroy();
        }

};
MultiSolverLBM3Test<tags::CPU, float> multi_solver_test_float("float");
