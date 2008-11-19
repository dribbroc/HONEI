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
#include <honei/lbm/solver_labswe.hh>
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
class SolverLABSWEGridRegressionTest :
    public TaggedTest<Tag_>
{
    public:
        SolverLABSWEGridRegressionTest(const std::string & type) :
            TaggedTest<Tag_>("solver_labswe_grid_regression_test<" + type + ">")
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
            DenseMatrix<DataType_> b_x(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> b_y(g_h, g_w, DataType_(0.));

            Grid<D2Q9, DataType_> grid;
            DenseMatrix<bool> obstacles(g_h, g_w, false);
            grid.obstacles = &obstacles;
            grid.h = &h;
            grid.u = &u;
            grid.v = &v;
            grid.b_x = &b_x;
            grid.b_y = &b_y;
            PackedGridData<D2Q9, DataType_>  data;
            PackedGridInfo<D2Q9> info;

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);

            SolverLABSWEGrid<Tag_, DataType_,lbm_source_types::CENTRED, lbm_source_schemes::CENTRALDIFF, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP> solver(&data, &info, 1., 1., 1.);

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
            //Standard solver not using grid:

            unsigned long g_h_standard(50);
            unsigned long g_w_standard(50);
            unsigned long timesteps_standard(200);

            DenseMatrix<DataType_> h_standard(g_h_standard, g_w_standard, DataType_(0.05));
            Cylinder<DataType_> c1_standard(h_standard, DataType_(0.02), 25, 25);
            c1_standard.value();

            DenseMatrix<DataType_> b_standard(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> u_standard(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> v_standard(g_h_standard, g_w_standard, DataType_(0.));

            //All needed distribution functions:

            DenseMatrix<DataType_> d_0(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> d_1(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> d_2(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> d_3(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> d_4(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> d_5(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> d_6(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> d_7(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> d_8(g_h_standard, g_w_standard, DataType_(0.));

            DenseMatrix<DataType_> e_d_0(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> e_d_1(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> e_d_2(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> e_d_3(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> e_d_4(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> e_d_5(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> e_d_6(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> e_d_7(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> e_d_8(g_h_standard, g_w_standard, DataType_(0.));

            DenseMatrix<DataType_> t_d_0(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> t_d_1(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> t_d_2(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> t_d_3(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> t_d_4(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> t_d_5(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> t_d_6(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> t_d_7(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> t_d_8(g_h_standard, g_w_standard, DataType_(0.));

            //All needed vectors:
            DenseVector<DataType_> v_x(9, DataType_(0));
            DenseVector<DataType_> v_y(9, DataType_(0));

            //Other matrices needed by solver:

            DenseMatrix<DataType_> s_x_standard(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> s_y_standard(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> d_x(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> d_y(g_h_standard, g_w_standard, DataType_(0.));

            SolverLABSWE<tags::CPU, DataType_,lbm_source_types::SIMPLE, lbm_source_schemes::BASIC, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP_PERIODIC> solver_standard(1.,1.,1., g_w_standard, g_h_standard, &h_standard, &b_standard, &u_standard, &v_standard);

            solver_standard.set_distribution(&d_0, &d_1, &d_2, &d_3, &d_4, &d_5, &d_6, &d_7, &d_8);
            solver_standard.set_eq_distribution(&e_d_0, &e_d_1, &e_d_2, &e_d_3, &e_d_4, &e_d_5, &e_d_6, &e_d_7, &e_d_8);
            solver_standard.set_temp_distribution(&t_d_0, &t_d_1, &t_d_2, &t_d_3, &t_d_4, &t_d_5, &t_d_6, &t_d_7, &t_d_8);
            solver_standard.set_vectors(&v_x, &v_y);
            solver_standard.set_source(&s_x_standard, &s_y_standard);
            solver_standard.set_slopes(&d_x, &d_y);
            solver_standard.do_preprocessing();

            for(unsigned long i(0); i < timesteps_standard; ++i)
            {
#ifdef SOLVER_VERBOSE
                std::cout<<"Timestep: " << i << "/" << timesteps_standard << std::endl;
#endif
                solver_standard.solve();
#ifdef SOLVER_POSTPROCESSING
                PostProcessing<GNUPLOT>::value(h_standard, 1, g_w_standard, g_h_standard, i);
#endif
            }
#ifdef SOLVER_VERBOSE
            std::cout << h_standard << std::endl;
#endif
            TEST_CHECK(true);

            TEST_CHECK_EQUAL(g_h, g_h_standard);
            TEST_CHECK_EQUAL(g_w, g_w_standard);

            //Compare CPU results of both solvers:
            for(unsigned long i(0) ; i < g_h ; ++i)
            {
                for(unsigned long j(0) ; j < g_w ; ++j)
                {
#ifdef SOLVER_VERBOSE
                    std::cout << "(" << i << " , " << j << ")" << std::endl;
                    std::cout << (*grid.h)(i , j) << " " << h_standard(i , j) << std::endl;
#endif
                    TEST_CHECK_EQUAL_WITHIN_EPS((*grid.h)(i , j) , h_standard(i , j), std::numeric_limits<float>::epsilon());
                }
            }

            //Save matrices to vectors, compute norm:
            DenseVector<DataType_> result_grid(g_h*g_w);
            DenseVector<DataType_> result_standard(g_h*g_w);

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
            DataType_ l2 = Norm<vnt_l_two, false, tags::CPU>::value(result_grid);
            TEST_CHECK_EQUAL_WITHIN_EPS(l2, DataType_(0.), std::numeric_limits<float>::epsilon());

            std::cout << "L2 norm: " << l2 << std::endl;
        }


};
SolverLABSWEGridRegressionTest<tags::CPU, float> solver_test_float("float");
SolverLABSWEGridRegressionTest<tags::CPU, double> solver_test_double("double");
SolverLABSWEGridRegressionTest<tags::CPU::MultiCore, float> mc_solver_test_float("float");
SolverLABSWEGridRegressionTest<tags::CPU::MultiCore, double> mc_solver_test_double("double");
#ifdef HONEI_CUDA
SolverLABSWEGridRegressionTest<tags::GPU::CUDA, float> cuda_solver_test_float("float");
#endif

