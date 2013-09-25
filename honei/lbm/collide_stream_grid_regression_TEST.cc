/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2009 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
#include <honei/swe/volume.hh>
#include <honei/util/unittest.hh>
#include <iostream>
#include <honei/swe/volume.hh>
#include <honei/math/quadrature.hh>
#include <honei/lbm/grid.hh>
#include <honei/lbm/grid_packer.hh>

using namespace honei;
using namespace tests;
using namespace std;
using namespace lbm::lbm_lattice_types;

//#define SOLVER_VERBOSE 1

template <typename Tag_, typename DataType_>
class CollideStreamGridRegressionTest :
    public TaggedTest<Tag_>
{
    public:
        CollideStreamGridRegressionTest(const std::string & type) :
            TaggedTest<Tag_>("collide_stream_grid_regression_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(50);
            unsigned long g_w(50);

            DenseMatrix<DataType_> h(g_h, g_w, DataType_(0.05));
            Cylinder<DataType_> c1(h, DataType_(0.02), 25, 25);
            c1.value();

            DenseMatrix<DataType_> u(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> v(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> b(g_h, g_w, DataType_(0.));

            Cylinder<DataType_> b1(b, DataType_(0.0001), 15, 15);
            b1.value();

            Grid<D2Q9, DataType_> grid;
            DenseMatrix<bool> obstacles(g_h, g_w, false);
            Cuboid<bool> q2(obstacles, 15, 5, 1, 10, 0);
            q2.value();
            Cuboid<bool> q3(obstacles, 40, 5, 1, 10, 30);
            q3.value();
            grid.obstacles = new DenseMatrix<bool>(obstacles);
            grid.h = new DenseMatrix<DataType_>(h);
            grid.u = new DenseMatrix<DataType_>(u);
            grid.v = new DenseMatrix<DataType_>(v);
            grid.b = new DenseMatrix<DataType_>(b);
            PackedGridData<D2Q9, DataType_>  data;
            PackedGridInfo<D2Q9> info;

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);

            SolverLBMGrid<Tag_, lbm_applications::LABSWE, DataType_,lbm_force::NONE, lbm_source_schemes::NONE, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::WET> solver(&info, &data, 0.01, 0.01, 0.005, 1.1);

            solver.do_preprocessing();
            solver.solve();

            CollideStreamGrid<Tag_, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>::
                value(info, data, DataType_(1.1));

            //Standard solver using tags::CPU:
            unsigned long g_h_standard(50);
            unsigned long g_w_standard(50);

            DenseMatrix<DataType_> h_standard(g_h_standard, g_w_standard, DataType_(0.05));
            Cylinder<DataType_> c1_standard(h_standard, DataType_(0.02), 25, 25);
            c1_standard.value();

            DenseMatrix<DataType_> u_standard(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> v_standard(g_h_standard, g_w_standard, DataType_(0.));
            DenseMatrix<DataType_> b_standard(g_h, g_w, DataType_(0.));

            Cylinder<DataType_> b1_standard(b_standard, DataType_(0.0001), 15, 15);
            b1_standard.value();

            Grid<D2Q9, DataType_> grid_standard;
            DenseMatrix<bool> obstacles_standard(g_h_standard, g_w_standard, false);
            Cuboid<bool> q2_standard(obstacles_standard, 15, 5, 1, 10, 0);
            q2_standard.value();
            Cuboid<bool> q3_standard(obstacles_standard, 40, 5, 1, 10, 30);
            q3_standard.value();
            grid_standard.obstacles = new DenseMatrix<bool>(obstacles_standard);
            grid_standard.h = new DenseMatrix<DataType_>(h_standard);
            grid_standard.u = new DenseMatrix<DataType_>(u_standard);
            grid_standard.v = new DenseMatrix<DataType_>(v_standard);
            grid_standard.b = new DenseMatrix<DataType_>(b_standard);
            PackedGridData<D2Q9, DataType_>  data_standard;
            PackedGridInfo<D2Q9> info_standard;

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid_standard, info_standard, data_standard);

            SolverLBMGrid<tags::CPU, lbm_applications::LABSWE, DataType_,lbm_force::NONE, lbm_source_schemes::NONE, lbm_grid_types::RECTANGULAR, lbm_lattice_types::D2Q9, lbm_boundary_types::NOSLIP, lbm_modes::WET> solver_standard(&info_standard, &data_standard, 0.01, 0.01, 0.005, 1.1);

            solver_standard.do_preprocessing();
            solver_standard.solve();

            CollideStreamGrid<tags::CPU, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>::
                value(info_standard, data_standard, 1.1);


            //Compare CPU results of both solvers:
            data.f_temp_0->lock(lm_read_only);
            data.f_temp_1->lock(lm_read_only);
            data.f_temp_2->lock(lm_read_only);
            data.f_temp_3->lock(lm_read_only);
            data.f_temp_4->lock(lm_read_only);
            data.f_temp_5->lock(lm_read_only);
            data.f_temp_6->lock(lm_read_only);
            data.f_temp_7->lock(lm_read_only);
            data.f_temp_8->lock(lm_read_only);
            for(unsigned long i(0) ; i < data.f_temp_0->size() ; ++i)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS((*data.f_temp_0)[i], (*data_standard.f_temp_0)[i], std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS((*data.f_temp_1)[i], (*data_standard.f_temp_1)[i], std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS((*data.f_temp_2)[i], (*data_standard.f_temp_2)[i], std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS((*data.f_temp_3)[i], (*data_standard.f_temp_3)[i], std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS((*data.f_temp_4)[i], (*data_standard.f_temp_4)[i], std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS((*data.f_temp_5)[i], (*data_standard.f_temp_5)[i], std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS((*data.f_temp_6)[i], (*data_standard.f_temp_6)[i], std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS((*data.f_temp_7)[i], (*data_standard.f_temp_7)[i], std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS((*data.f_temp_8)[i], (*data_standard.f_temp_8)[i], std::numeric_limits<DataType_>::epsilon());
            }
            data.f_temp_0->unlock(lm_read_only);
            data.f_temp_1->unlock(lm_read_only);
            data.f_temp_2->unlock(lm_read_only);
            data.f_temp_3->unlock(lm_read_only);
            data.f_temp_4->unlock(lm_read_only);
            data.f_temp_5->unlock(lm_read_only);
            data.f_temp_6->unlock(lm_read_only);
            data.f_temp_7->unlock(lm_read_only);
            data.f_temp_8->unlock(lm_read_only);

            grid.destroy();
            info.destroy();
            data.destroy();
            grid_standard.destroy();
            info_standard.destroy();
            data_standard.destroy();
        }
};
//CollideStreamGridRegressionTest<tags::CPU::Generic, float> generic_collide_stream_regression_test_float("float");
//CollideStreamGridRegressionTest<tags::CPU::Generic, double> generic_collide_stream_regression_test_double("double");
#ifdef HONEI_SSE
CollideStreamGridRegressionTest<tags::CPU::SSE, float> sse_collide_stream_regression_test_float("float");
CollideStreamGridRegressionTest<tags::CPU::SSE, double> sse_collide_stream_regression_test_double("double");
#endif
#ifdef HONEI_CUDA
CollideStreamGridRegressionTest<tags::GPU::CUDA, float> cuda_collide_stream_regression_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
CollideStreamGridRegressionTest<tags::GPU::CUDA, double> cuda_collide_stream_regression_test_double("double");
#endif
#endif
#ifdef HONEI_CELL
CollideStreamGridRegressionTest<tags::Cell, float> cell_collide_stream_regression_test_float("float");
#endif
#ifdef HONEI_OPENCL
CollideStreamGridRegressionTest<tags::OpenCL::CPU, float> ocl_cpu_collide_stream_regression_test_float("float");
CollideStreamGridRegressionTest<tags::OpenCL::CPU, double> ocl_cpu_collide_stream_regression_test_double("double");
#endif
#ifdef HONEI_OPENCL_GPU
CollideStreamGridRegressionTest<tags::OpenCL::GPU, float> ocl_gpu_collide_stream_regression_test_float("float");
CollideStreamGridRegressionTest<tags::OpenCL::GPU, double> ocl_gpu_collide_stream_regression_test_double("double");
#endif
#ifdef HONEI_OPENCL_ACC
CollideStreamGridRegressionTest<tags::OpenCL::Accelerator, float> ocl_acc_collide_stream_regression_test_float("float");
CollideStreamGridRegressionTest<tags::OpenCL::Accelerator, double> ocl_acc_collide_stream_regression_test_double("double");
#endif
