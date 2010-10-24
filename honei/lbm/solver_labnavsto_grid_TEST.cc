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
#include <honei/util/unittest.hh>
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
class SolverLABNAVSTOGridTest :
    public TaggedTest<Tag_>
{
    public:
        SolverLABNAVSTOGridTest(const std::string & type) :
            TaggedTest<Tag_>("solver_lbm_grid_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(50);
            unsigned long g_w(50);
            unsigned long timesteps(100);


            for (unsigned long i(0) ; i < ScenarioCollection::get_scenario_count() - 2 ; ++i)
            {
                Grid<D2Q9, DataType_> grid;

                ScenarioCollection::get_scenario(i, g_h, g_w, grid);

                PackedGridData<D2Q9, DataType_>  data;
                PackedGridInfo<D2Q9> info;

                GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);

                SolverLBMGrid<Tag_, lbm_applications::LABNAVSTO, DataType_,lbm_force::NONE, lbm_source_schemes::NONE, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> solver(&info, &data, grid.d_x, grid.d_y, grid.d_t, grid.tau);

                solver.do_preprocessing();
                std::cout << "Solving: " << grid.description << std::endl;
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

SolverLABNAVSTOGridTest<tags::CPU, float> solver_test_float("float");
SolverLABNAVSTOGridTest<tags::CPU, double> solver_test_double("double");
SolverLABNAVSTOGridTest<tags::CPU::MultiCore, float> mc_solver_test_float("float");
SolverLABNAVSTOGridTest<tags::CPU::MultiCore, double> mc_solver_test_double("double");


template <typename Tag_, typename DataType_>
class SolverLABNAVSTOGridMassConservationTest :
    public TaggedTest<Tag_>
{
    public:
        SolverLABNAVSTOGridMassConservationTest(const std::string & type) :
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

            SolverLBMGrid<Tag_, lbm_applications::LABNAVSTO, DataType_, lbm_force::NONE, lbm_source_schemes::NONE, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::DRY> solver(&info, &data, 1., 1., 1., 1.5);

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
            GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
            DataType_ ana_vol(5. * 5. * 0.02 + g_w * g_h * 0.05);
            std::cout << "Analytical Vol.: " << ana_vol << " ";
            DataType_ vol = GaussianQuadrature2D<tags::CPU, tags::Trapezoid>::value(h, DataType_(0), DataType_(g_w), DataType_(1.), DataType_(1.));
            std::cout << "Vol.: " << vol << std::endl;
            TEST_CHECK_EQUAL_WITHIN_EPS(vol, ana_vol, 0.1);
        }

};
SolverLABNAVSTOGridMassConservationTest<tags::CPU, float> solver_grid_mc_test_float("float");
SolverLABNAVSTOGridMassConservationTest<tags::CPU, double> solver_grid_mc_test_double("double");
SolverLABNAVSTOGridMassConservationTest<tags::CPU::MultiCore, float> mc_solver_grid_mc_test_float("float");
SolverLABNAVSTOGridMassConservationTest<tags::CPU::MultiCore, double> mc_solver_grid_mc_test_double("double");
