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
#include <honei/lbm/solver_labswe_grid.hh>
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
class SolverLABNAVSTOGridDCTest :
    public TaggedTest<Tag_>
{
    public:
        SolverLABNAVSTOGridDCTest(const std::string & type) :
            TaggedTest<Tag_>("solver_labswe_grid_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(129);
            unsigned long g_w(129);
            unsigned long timesteps(1000);
            DataType_ lid_U(0.1);
            DataType_ dx(1);
            DataType_ dy(1);
            DataType_ dt(1);
            DataType_ tau(1.5);

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


            SolverLABSWEGrid<Tag_, lbm_applications::LABNAVSTO, DataType_,lbm_force::NONE, lbm_source_schemes::NONE, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::WET> solver(&info, &data, grid.d_x, grid.d_y, grid.d_t, grid.tau);

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
                    value(9.81, (grid.d_x / grid.d_t) * (grid.d_x / grid.d_t), info, data);

                CollideStreamGrid<Tag_, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>::
                    value(info,
                          data,
                          tau);
                ///Boundary correction:
                UpdateVelocityDirectionsGrid<Tag_, lbm_boundary_types::NOSLIP>::
                    value(data, info);

                //extract velocities out of h from previous timestep:
                ExtractionGrid<Tag_, lbm_modes::WET>::value(info, data, DataType_(10e-5));

#ifdef SOLVER_POSTPROCESSING
                GridPacker<D2Q9, NOSLIP, DataType_>::unpack_u(grid, info, data);
                PostProcessing<GNUPLOT>::value(*grid.u, 1, g_w, g_h, i);
#endif
            }
            GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
#ifdef SOLVER_VERBOSE
            std::cout << *grid.h << std::endl;
#endif
            for (unsigned long i(0) ; i < (*grid.h).rows() ; ++i)
                for(unsigned long j(0) ; j < (*grid.h).columns() ; ++j)
                    TEST_CHECK_EQUAL_WITHIN_EPS((*grid.h)( i , j), DataType_(0.02), DataType_(0.1));

        }
};

SolverLABNAVSTOGridDCTest<tags::CPU, double> solver_test_double("double");
