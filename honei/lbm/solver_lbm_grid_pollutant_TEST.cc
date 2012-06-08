/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2012 Markus Geveler <apryde@gmx.de>
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
#include <honei/lbm/solver_lbm_grid_pollutant.hh>
#include <honei/lbm/partial_derivative.hh>
#include <honei/swe/post_processing.hh>
#include <honei/lbm/post_processing.hh>
#include <honei/util/unittest.hh>
#include <iostream>
#include <honei/swe/volume.hh>
#include <honei/math/quadrature.hh>
#include <honei/lbm/grid.hh>
#include <honei/lbm/grid_packer.hh>
#include <honei/lbm/scenario_collection.hh>
#include <honei/lbm/extraction_grid_pollutant.hh>
#include <honei/lbm/equilibrium_distribution_grid_pollutant.hh>
#include <honei/lbm/force_grid_pollutant.hh>
#include <honei/lbm/collide_stream_grid_pollutant.hh>
#include <honei/lbm/function.hh>

#ifdef DEBUG
#define SOLVER_VERBOSE
#endif

using namespace honei;
using namespace tests;
using namespace std;
using namespace output_types;
using namespace lbm::lbm_lattice_types;

//#define SOLVER_VERBOSE
#define SOLVER_POSTPROCESSING

template <typename Tag_, typename DataType_>
class SolverLBMGridPollutantTest :
    public TaggedTest<Tag_>
{
    public:
        SolverLBMGridPollutantTest(const std::string & type) :
            TaggedTest<Tag_>("solver_lbm_grid_pollutant_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(50);
            unsigned long g_w(50);
            unsigned long timesteps(100);
            DataType_ tau_c(1.5);
            DataType_ k(0.0005);
            DataType_ s_0(0.);

            Grid<D2Q9, DataType_> grid;
            Grid<D2Q9, DataType_> grid_poll;

            ScenarioCollection::get_scenario(6, g_h, g_w, grid);
            ScenarioCollection::get_scenario(6, g_h, g_w, grid_poll);
            delete grid_poll.h;
            grid_poll.h = new DenseMatrix<DataType_>(g_h, g_w, DataType_(0.1));

            PackedGridData<D2Q9, DataType_>  data_flow;
            PackedGridData<D2Q9, DataType_>  data_poll;
            PackedGridInfo<D2Q9> info;

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data_flow);
            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid_poll, info, data_poll);

            std::cout << "D = " << (grid_poll.d_t / 6.) * (2 * tau_c - 1) * (grid_poll.d_x / grid_poll.d_t) * (grid_poll.d_x / grid_poll.d_t) << std::endl;

            SolverLBMGrid<Tag_, lbm_applications::LABSWE, DataType_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> solver(&info, &data_flow, grid.d_x, grid.d_y, grid.d_t, grid.tau);

            SolverLBMGridPollutant<Tag_, lbm_applications::LABSWE, DataType_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> poll_solver(&info, &data_flow, &data_poll, grid_poll.d_x, grid_poll.d_y, grid_poll.d_t, DataType_(tau_c), DataType_(k), DataType_(s_0), 0.);

            solver.do_preprocessing();
            poll_solver.do_preprocessing();

            std::cout << "Solving: " << grid.description << std::endl;
            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                solver.solve();

                poll_solver.solve();
#ifdef SOLVER_POSTPROCESSING
                solver.do_postprocessing();
                poll_solver.do_postprocessing();
                //GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data_flow);
                //PostProcessing<GNUPLOT>::value(*grid.h, 1, g_w, g_h, i);
                GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid_poll, info, data_poll);
                PostProcessing<GNUPLOT>::value(*grid_poll.h, 1, g_w, g_h, i);

                //if(i==0)
                //    std::cout << *grid_poll.h;
#endif
            }
            solver.do_postprocessing();
            poll_solver.do_postprocessing();
            GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data_flow);
#ifdef SOLVER_VERBOSE
            std::cout << *grid.h << std::endl;
            std::cout << *grid_poll.h << std::endl;
#endif
            /*for (unsigned long i(0) ; i < (*grid.h).rows() ; ++i)
                for(unsigned long j(0) ; j < (*grid.h).columns() ; ++j)
                    TEST_CHECK_EQUAL_WITHIN_EPS((*grid.h)( i , j), DataType_(0.02), DataType_(0.1));
            */
            info.destroy();
            data_flow.destroy();
            data_poll.destroy();
            grid.destroy();
            grid_poll.destroy();
        }

};
//SolverLBMGridPollutantTest<tags::CPU::Generic, float> solver_test_float("float");

template <typename Tag_, typename DataType_>
class SolverLBMGridPollutantFlowTest :
    public TaggedTest<Tag_>
{
    public:
        SolverLBMGridPollutantFlowTest(const std::string & type) :
            TaggedTest<Tag_>("solver_lbm_grid_pollutant_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(50);
            unsigned long g_w(50);
            unsigned long timesteps(1000);
            DataType_ tau_c(1);
            DataType_ k(0.5);
            DataType_ s_0(0.);

            Grid<D2Q9, DataType_> grid;
            Grid<D2Q9, DataType_> grid_poll;

            ScenarioCollection::get_scenario(10, g_h, g_w, grid);
            ScenarioCollection::get_scenario(10, g_h, g_w, grid_poll);
            delete grid_poll.h;
            grid_poll.h = new DenseMatrix<DataType_>(g_h, g_w, DataType_(0.1));

            PackedGridData<D2Q9, DataType_>  data_flow;
            PackedGridData<D2Q9, DataType_>  data_poll;
            PackedGridInfo<D2Q9> info;

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data_flow);
            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid_poll, info, data_poll);

            std::cout << "D = " << (grid_poll.d_t / 6.) * (2 * tau_c) * (grid_poll.d_x / grid_poll.d_t) * (grid_poll.d_x / grid_poll.d_t) << std::endl;

            SolverLBMGrid<Tag_, lbm_applications::LABSWE, DataType_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> solver(&info, &data_flow, grid.d_x, grid.d_y, grid.d_t, grid.tau);

            SolverLBMGridPollutant<Tag_, lbm_applications::LABSWE, DataType_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> poll_solver(&info, &data_flow, &data_poll, grid_poll.d_x, grid_poll.d_y, grid_poll.d_t, DataType_(tau_c), DataType_(k), DataType_(s_0), 0.1);

            solver.do_preprocessing();
            poll_solver.do_preprocessing();

            std::cout << "Solving: " << grid.description << std::endl;
            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                solver.solve();

                poll_solver.solve();
#ifdef SOLVER_POSTPROCESSING
                solver.do_postprocessing();
                poll_solver.do_postprocessing();
                GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data_flow);
                PostProcessing<GNUPLOT>::value(*grid.h, 10, g_w, g_h, i, "flow" + stringify(i) + ".dat");
                GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid_poll, info, data_poll);
                PostProcessing<GNUPLOT>::value(*grid_poll.h, 10, g_w, g_h, i, "poll" + stringify(i) + ".dat");

                //if(i==0)
                //    std::cout << *grid_poll.h;
#endif
            }
            solver.do_postprocessing();
            poll_solver.do_postprocessing();
            GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data_flow);
#ifdef SOLVER_VERBOSE
            std::cout << *grid.h << std::endl;
            std::cout << *grid_poll.h << std::endl;
#endif
            /*for (unsigned long i(0) ; i < (*grid.h).rows() ; ++i)
                for(unsigned long j(0) ; j < (*grid.h).columns() ; ++j)
                    TEST_CHECK_EQUAL_WITHIN_EPS((*grid.h)( i , j), DataType_(0.02), DataType_(0.1));
            */
            info.destroy();
            data_flow.destroy();
            data_poll.destroy();
            grid.destroy();
            grid_poll.destroy();
        }

};
//SolverLBMGridPollutantFlowTest<tags::CPU::Generic, float> solver_test_float("float");

template <typename Tag_, typename DataType_>
class Showcase :
    public TaggedTest<Tag_>
{
    public:
        Showcase(const std::string & type) :
            TaggedTest<Tag_>("solver_lbm_grid_pollutant_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned print_every(10);
            unsigned long g_h(50);
            unsigned long g_w(50);
            unsigned long timesteps(5000);
            DataType_ tau_c(1.);
            DataType_ k(1e-10);
            DataType_ s_0(0.);

            Grid<D2Q9, DataType_> grid;
            Grid<D2Q9, DataType_> grid_poll;

            ScenarioCollection::get_scenario(10, g_h, g_w, grid);
            //override flow parameters
            grid.d_x = 5.;
            grid.d_y = 5.;
            grid.d_t = 1.;
            grid.tau = 0.9;

            delete grid.h;
            grid.h = new DenseMatrix<DataType_>(g_h, g_w, DataType_(1.));
            /*for(unsigned long i(0); i < g_h ; ++i)
                for(unsigned long j(0); j < g_w ; ++j)
                {
                    (*grid.h)[i][j] += Gaussian2D<90>::value(grid.d_x * i,  //x
                                                             grid.d_y * j,  //y
                                                             1.f,          //amplitude
                                                             grid.d_x * 25, //center x
                                                             grid.d_y * 15, //center y
                                                             35.f,          //x spread
                                                             35.f);         //y spread
                }*/

            ScenarioCollection::get_scenario(10, g_h, g_w, grid_poll);
            grid_poll.d_x = grid.d_x;
            grid_poll.d_y = grid.d_y;
            grid_poll.d_t = grid.d_t;
            delete grid_poll.h;
            grid_poll.h = new DenseMatrix<DataType_>(g_h, g_w, DataType_(0.));
            for(unsigned long i(0); i < g_h ; ++i)
                for(unsigned long j(0); j < g_w ; ++j)
                   (*grid_poll.h)[i][j] = Gaussian2D<0>::value(grid_poll.d_x * i,
                                                               grid_poll.d_y * j,
                                                               0.5f,
                                                               grid_poll.d_x * 25,
                                                               grid_poll.d_y * 25,
                                                               3.f,
                                                               3.f);

            //LBMPostProcessing::value(*grid_poll.h, 1, g_w, g_h, grid.d_x, grid.d_y, 1, "poll_init.dat");

            PackedGridData<D2Q9, DataType_>  data_flow;
            PackedGridData<D2Q9, DataType_>  data_poll;
            PackedGridInfo<D2Q9> info;

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data_flow);
            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid_poll, info, data_poll);

            std::cout << "D = " << (grid_poll.d_t / 6.) * (2 * tau_c - 1) * (grid_poll.d_x / grid_poll.d_t) * (grid_poll.d_x / grid_poll.d_t) << std::endl;

            SolverLBMGrid<Tag_, lbm_applications::LABSWE, DataType_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> solver(&info, &data_flow, grid.d_x, grid.d_y, grid.d_t, grid.tau);

            SolverLBMGridPollutant<Tag_, lbm_applications::LABSWE, DataType_,lbm_force::CENTRED, lbm_source_schemes::BED_FULL, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> poll_solver(&info, &data_flow, &data_poll, grid_poll.d_x, grid_poll.d_y, grid_poll.d_t, DataType_(tau_c), DataType_(k), DataType_(s_0), 0.005);

            DenseVector<DataType_> const_u(data_flow.u->size(), DataType_(0.25));
            DenseVector<DataType_> const_v(data_flow.v->size(), DataType_(0));

            *data_flow.u = const_u.copy();
            *data_flow.v = const_v.copy();

            solver.do_preprocessing();
            poll_solver.do_preprocessing();

            std::cout << "Solving: " << grid.description << std::endl;
            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                *data_flow.u = const_u.copy();
                *data_flow.v = const_v.copy();
                solver.solve();
                *data_flow.u = const_u.copy();
                *data_flow.v = const_v.copy();
                poll_solver.solve();
#ifdef SOLVER_POSTPROCESSING
                solver.do_postprocessing();
                poll_solver.do_postprocessing();
                GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data_flow);

                unsigned I(i), digits(0);
                while(I > 0)
                {
                    digits++;
                    I/=10;
                }
                std::string i_;
                for(unsigned j(0) ; j < (i != 0 ? 4 : 3) - digits ; j++)
                {
                    i_ += "0";
                }
                i_ += stringify(i);
                std::cout << i_ << std::endl;

                LBMPostProcessing::value(*grid.h, print_every, g_w, g_h, grid.d_x, grid.d_y, i, "flow" + i_ + ".dat");
                GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid_poll, info, data_poll);
                LBMPostProcessing::value(*grid_poll.h, print_every, g_w, g_h, grid.d_x, grid.d_y, i, "poll" + i_ + ".dat");

                //if(i==0)
                //    std::cout << *grid_poll.h;
#endif
            }
            solver.do_postprocessing();
            poll_solver.do_postprocessing();
            GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data_flow);
#ifdef SOLVER_VERBOSE
            std::cout << *grid.h << std::endl;
            std::cout << *grid_poll.h << std::endl;
#endif
            /*for (unsigned long i(0) ; i < (*grid.h).rows() ; ++i)
                for(unsigned long j(0) ; j < (*grid.h).columns() ; ++j)
                    TEST_CHECK_EQUAL_WITHIN_EPS((*grid.h)( i , j), DataType_(0.02), DataType_(0.1));
            */
            info.destroy();
            data_flow.destroy();
            data_poll.destroy();
            grid.destroy();
            grid_poll.destroy();
        }
};
Showcase<tags::GPU::CUDA, float> solver_test_float("float");
