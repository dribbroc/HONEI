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
#include <honei/lbm/solver_labswe_grid.hh>
#include <honei/lbm/partial_derivative.hh>
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

using namespace honei;
using namespace tests;
using namespace std;
using namespace output_types;
using namespace lbm::lbm_lattice_types;

//#define SOLVER_VERBOSE 1

template <typename Tag_, typename DataType_>
class SolverLABSWEGridMultiRegressionTest :
    public TaggedTest<Tag_>
{
    public:
        SolverLABSWEGridMultiRegressionTest(const std::string & type) :
            TaggedTest<Tag_>("solver_labswe_grid_multi_regression_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(50);
            unsigned long g_w(50);
            unsigned long timesteps(200);

            DenseMatrix<DataType_> h(g_h, g_w, DataType_(0.05));
            Cylinder<DataType_> c1(h, DataType_(0.02), 25, 25);
            c1.value();

            DenseMatrix<DataType_> u(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> v(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> b(g_h, g_w, DataType_(0.));

            Cylinder<DataType_> b1(b, DataType_(0.001), 15, 15);
            b1.value();

            Grid<D2Q9, DataType_> grid;
            DenseMatrix<bool> obstacles(g_h, g_w, false);
            Cuboid<bool> q2(obstacles, 15, 5, 1, 10, 0);
            q2.value();
            Cuboid<bool> q3(obstacles, 40, 5, 1, 10, 30);
            q3.value();
            grid.obstacles = &obstacles;
            grid.h = &h;
            grid.u = &u;
            grid.v = &v;
            grid.b = &b;
            PackedGridData<D2Q9, DataType_>  data;
            PackedGridInfo<D2Q9> info;

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);

            SolverLABSWEGrid<Tag_, DataType_,lbm_force::CENTRED, lbm_source_schemes::BED_SLOPE, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP> solver(&data, &info, 1., 1., 1., 1.5);

            solver.do_preprocessing();

            for(unsigned long i(0); i < timesteps; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                solver.solve();
#ifdef SOLVER_POSTPROCESSING
                GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);
                PostProcessing<GNUPLOT>::value(h, 1, g_w, g_h, i);
#endif
            }
            GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid, info, data);

#ifdef SOLVER_VERBOSE
            std::cout << *grid.h << std::endl;
#endif
            //Standard solver not using partitioner:
            unsigned long g_h_standard(50);
            unsigned long g_w_standard(50);
            unsigned long timesteps_standard(200);

            DenseMatrix<DataType_> h_standard(g_h_standard, g_w_standard, DataType_(0.05));
            Cylinder<DataType_> c1_standard(h_standard, DataType_(0.02), 25, 25);
            c1_standard.value();

            DenseMatrix<DataType_> u_standard(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> v_standard(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> b_standard(g_h, g_w, DataType_(0.));

            Cylinder<DataType_> b1_standard(b_standard, DataType_(0.001), 15, 15);
            b1_standard.value();

            Grid<D2Q9, DataType_> grid_standard;
            DenseMatrix<bool> obstacles_standard(g_h_standard, g_w_standard, false);
            Cuboid<bool> q2_standard(obstacles_standard, 15, 5, 1, 10, 0);
            q2_standard.value();
            Cuboid<bool> q3_standard(obstacles_standard, 40, 5, 1, 10, 30);
            q3_standard.value();
            grid_standard.obstacles = &obstacles_standard;
            grid_standard.h = &h_standard;
            grid_standard.u = &u_standard;
            grid_standard.v = &v_standard;
            grid_standard.b = &b_standard;
            PackedGridData<D2Q9, DataType_>  data_standard;
            PackedGridInfo<D2Q9> info_standard;

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid_standard, info_standard, data_standard);

            SolverLABSWEGrid<tags::CPU, DataType_,lbm_force::CENTRED, lbm_source_schemes::BED_SLOPE, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP> solver_standard(&data_standard, &info_standard, 1., 1., 1., 1.5);

            solver_standard.do_preprocessing();

            for(unsigned long i(0); i < timesteps_standard; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps << std::endl;
#endif
                solver_standard.solve();
#ifdef SOLVER_POSTPROCESSING
                GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid_standard, info_standard, data_standard);
                PostProcessing<GNUPLOT>::value(h_standard, 1, g_w_standard, g_h_standard, i);
#endif
            }
            GridPacker<D2Q9, NOSLIP, DataType_>::unpack(grid_standard, info_standard, data_standard);

#ifdef SOLVER_VERBOSE
            std::cout << *grid_standard.h << std::endl;
#endif


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
                    TEST_CHECK_EQUAL_WITHIN_EPS((*grid.h)(i , j) , (*grid_standard.h)(i , j), std::numeric_limits<DataType_>::epsilon());
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
                    result_standard[inner] = h_standard(i , j);
                    ++inner;
                }
            }

            Difference<tags::CPU>::value(result_grid, result_standard);
            double l2 = Norm<vnt_l_two, false, tags::CPU>::value(result_grid);
            TEST_CHECK_EQUAL_WITHIN_EPS(l2, DataType_(0.), std::numeric_limits<DataType_>::epsilon());

            std::cout << "L2 norm: " << l2 << std::endl;
        }


};
SolverLABSWEGridMultiRegressionTest<tags::CPU::MultiCore, float> mc_solver_multi_test_float("float");
SolverLABSWEGridMultiRegressionTest<tags::CPU::MultiCore, double> mc_solver_multi_test_double("double");
#ifdef HONEI_SSE
SolverLABSWEGridMultiRegressionTest<tags::CPU::SSE, float> sse_solver_multi_test_float("float");
SolverLABSWEGridMultiRegressionTest<tags::CPU::SSE, double> sse_solver_multi_test_double("double");
SolverLABSWEGridMultiRegressionTest<tags::CPU::MultiCore::SSE, float> mcsse_solver_multi_test_float("float");
SolverLABSWEGridMultiRegressionTest<tags::CPU::MultiCore::SSE, double> mcsse_solver_multi_test_double("double");
#endif
#ifdef HONEI_CUDA
SolverLABSWEGridMultiRegressionTest<tags::GPU::CUDA, float> cuda_solver_multi_test_float("float");
#endif

