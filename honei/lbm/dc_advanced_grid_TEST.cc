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
#include <honei/lbm/partial_derivative.hh>
#include <honei/swe/post_processing.hh>
#include <unittest/unittest.hh>
#include <iostream>
#include <honei/swe/volume.hh>
#include <honei/math/quadrature.hh>
#include <honei/lbm/grid.hh>
#include <honei/lbm/grid_packer.hh>
#include <honei/lbm/scenario_collection.hh>
#include <honei/lbm/dc_util.hh>

using namespace honei;
using namespace tests;
using namespace std;
using namespace output_types;
using namespace lbm::lbm_lattice_types;
using namespace lbm;

//#define SOLVER_VERBOSE
//#define SOLVER_POSTPROCESSING

template <typename Tag_, typename DataType_>
class SolverLABNAVSTOGridDCTest :
    public TaggedTest<Tag_>
{
    public:
        SolverLABNAVSTOGridDCTest(const std::string & type) :
            TaggedTest<Tag_>("solver_lbm_grid_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(129);
            unsigned long g_w(129);
            unsigned long timesteps(10000);
            DataType_ dx(1);
            DataType_ dy(1);
            DataType_ dt(1.2);
            DataType_ tau(1);
            DataType_ lid_U(Reynolds::adjust_veloc(double(100.), double(0.), g_w, dx, dt, tau));

            std::cout << "U: " << lid_U << std::endl;
            std::cout << "Reynolds: " << Reynolds::value(lid_U, g_w, dx, dt, tau) << std::endl;

            Grid<D2Q9, DataType_> grid;

            //just take some scenario and adjust it;
            ScenarioCollection::get_scenario(7, g_h, g_w, grid);
            grid.d_x = dx;
            grid.d_y = dy;
            grid.d_t = dt;
            grid.tau = tau;

            //init lid velocity
            for (unsigned long i(0) ; i < g_w ; ++i)
            {
                (*grid.u)[0][i] = lid_U;
            }

            PackedGridData<D2Q9, DataType_>  data;
            PackedGridInfo<D2Q9> info;

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);


            SolverLBMGrid<Tag_, lbm_applications::LABNAVSTO, DataType_,lbm_force::NONE, lbm_source_schemes::NONE, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::WET> solver(&info, &data, grid.d_x, grid.d_y, grid.d_t, grid.tau);

            solver.do_preprocessing();
            std::cout << "Solving: " << grid.description << std::endl;
            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                //after solving, reset lid velocity by using h_index
                for (unsigned long j(0) ; j < g_w ; ++j)
                {
                    (*data.u)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, 0, j)] = lid_U;
                    (*data.v)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, 0, j)] = DataType_(0.);

                    (*data.u)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, g_h - 1, j)] = DataType_(0.);
                    (*data.v)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, g_h - 1, j)] = DataType_(0.);

                    if(j > 0)
                    {
                        (*data.u)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, j, 0)] = DataType_(0.);
                        (*data.v)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, j, 0)] = DataType_(0.);

                        (*data.u)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, j, g_w - 1)] = DataType_(0.);
                        (*data.v)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, j, g_w - 1)] = DataType_(0.);
                    }
                }

                //solver.solve();
                EquilibriumDistributionGrid<Tag_, lbm_applications::LABNAVSTO>::
                    value(DataType_(9.81), (grid.d_x / grid.d_t) * (grid.d_x / grid.d_t), info, data);

                CollideStreamGrid<Tag_, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>::
                    value(info,
                          data,
                          tau);
                ///Boundary correction:
                UpdateVelocityDirectionsGrid<Tag_, lbm_boundary_types::NOSLIP>::
                    value(info, data);

                //extract velocities out of h from previous timestep:
                ExtractionGrid<Tag_, lbm_modes::WET>::value(info, data, DataType_(10e-5));

                //after solving, reset lid velocity by using h_index
                for (unsigned long j(0) ; j < g_w ; ++j)
                {
                    (*data.u)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, 0, j)] = lid_U;
                    (*data.v)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, 0, j)] = DataType_(0.);

                    (*data.u)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, g_h - 1, j)] = DataType_(0.);
                    (*data.v)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, g_h - 1, j)] = DataType_(0.);

                    if(j > 0)
                    {
                        (*data.u)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, j, 0)] = DataType_(0.);
                        (*data.v)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, j, 0)] = DataType_(0.);

                        (*data.u)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, j, g_w - 1)] = DataType_(0.);
                        (*data.v)[GridPacker<D2Q9, NOSLIP, DataType_>::h_index(grid, j, g_w - 1)] = DataType_(0.);
                    }
                }
#ifdef SOLVER_POSTPROCESSING
                GridPacker<D2Q9, NOSLIP, DataType_>::unpack_u(grid, info, data);
                PostProcessing<GNUPLOT>::value(*grid.u, 99, g_w, g_h, i);
#endif
            }
            GridPacker<D2Q9, NOSLIP, DataType_>::unpack_u(grid, info, data);
            DenseVector<DataType_> test_line(g_h);
            std::cout<<"Index: " << (g_w - 1) / 2 << std::endl;
            for(unsigned long i(0); i < g_h; ++i)
            {
                test_line[i] = DataType_( (*grid.u)(i,( g_w - 1 ) / 2)/ lid_U);
            }
            std::cout<<"Result:"<<test_line<<std::endl;

            //Reference data by Ghia et al. 1982 for reynolds number of 100:
            DenseVector<double> ref_result_100(17);
            unsigned long indices_100[17];

            ref_result_100[0] = double(1);
            indices_100[0] = 0;

            ref_result_100[1] = double(0.84123);
            indices_100[1] = 3;

            ref_result_100[2] = double(0.78871);
            indices_100[2] = 4;

            ref_result_100[3] = double(0.73722);
            indices_100[3] = 5;

            ref_result_100[4] = double(0.68717);
            indices_100[4] = 6;

            ref_result_100[5] = double(0.23151);
            indices_100[5] = 19;

            ref_result_100[6] = double(0.00332);
            indices_100[6] = 34;

            ref_result_100[7] = double(-0.13641);
            indices_100[7] = 49;

            ref_result_100[8] = double(-0.20581);
            indices_100[8] = 64;

            ref_result_100[9] = double(-0.21090);
            indices_100[9] = 70;

            ref_result_100[10] = double(-0.15662);
            indices_100[10] = 92;

            ref_result_100[11] = double(-0.10150);
            indices_100[11] = 106;

            ref_result_100[12] = double(-0.06434);
            indices_100[12] = 115;

            ref_result_100[13] = double(-0.04775);
            indices_100[13] = 119;

            ref_result_100[14] = double(-0.04192);
            indices_100[14] = 120;

            ref_result_100[15] = double(-0.03717);
            indices_100[15] = 121;

            ref_result_100[16] = double(0.);
            indices_100[16] = 128;

            DenseVector<DataType_> diff(test_line.copy());

            for(unsigned long i(0); i < 17; ++i)
            {
                diff[indices_100[i]] = ref_result_100[i];
                //TEST_CHECK_EQUAL_WITHIN_EPS(test_line[indices_100[i]], ref_result_100[i], 0.5);
            }

            Difference<>::value(diff, test_line);

            std::cout <<"Difference vector: " << diff << std::endl;

            double norm = Norm<vnt_l_two, false, Tag_>::value(diff);
            std::cout << "L2 norm: " << norm << std::endl;
        }
};

/*SolverLABNAVSTOGridDCTest<tags::CPU, double> solver_test_double("double");*/
SolverLABNAVSTOGridDCTest<tags::CPU, float> solver_test_float("float");
