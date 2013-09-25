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
#include <honei/lbm/tags.hh>
#include <honei/util/unittest.hh>
#include <honei/lbm/grid.hh>
#include <honei/lbm/collide_stream_grid.hh>
#include <honei/lbm/grid_packer.hh>
using namespace honei;
using namespace tests;
using namespace std;
using namespace lbm;
using namespace lbm_boundary_types;

template <typename Tag_, typename DataType_>
class CollideStreamGridLABSWETest :
    public TaggedTest<Tag_>
{
    public:
        CollideStreamGridLABSWETest(const std::string & type) :
            TaggedTest<Tag_>("collideandstream_grid_labswe_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(25);
            unsigned long g_w(25);

            DenseMatrix<DataType_> h(g_h, g_w, DataType_(0.05));

            DenseMatrix<DataType_> u(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> v(g_h, g_w, DataType_(0.));

            DenseMatrix<DataType_> b(g_h, g_w, DataType_(0.));

            Grid<D2Q9, DataType_> grid;
            DenseMatrix<bool> obstacles(g_h, g_w, false);
            grid.obstacles = new DenseMatrix<bool>(obstacles);
            grid.h = new DenseMatrix<DataType_>(h);
            grid.u = new DenseMatrix<DataType_>(u);
            grid.v = new DenseMatrix<DataType_>(v);
            grid.b = new DenseMatrix<DataType_>(b);

            PackedGridData<D2Q9, DataType_>  data;
            PackedGridInfo<D2Q9> info;

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);
            GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DataType_>::cuda_pack(info, data);

            DataType_ tau (1);

            for(unsigned long i(0); i < data.h->size(); i++)
            {
                (*data.f_eq_0)[i] = DataType_(1.234);
                (*data.f_eq_1)[i] = DataType_(1.234);
                (*data.f_eq_2)[i] = DataType_(1.234);
                (*data.f_eq_3)[i] = DataType_(1.234);
                (*data.f_eq_4)[i] = DataType_(1.234);
                (*data.f_eq_5)[i] = DataType_(1.234);
                (*data.f_eq_6)[i] = DataType_(1.234);
                (*data.f_eq_7)[i] = DataType_(1.234);
                (*data.f_eq_8)[i] = DataType_(1.234);
                (*data.f_0)[i] = DataType_(2.234);
                (*data.f_1)[i] = DataType_(2.234);
                (*data.f_2)[i] = DataType_(11.234);
                (*data.f_3)[i] = DataType_(2.234);
                (*data.f_4)[i] = DataType_(2.234);
                (*data.f_5)[i] = DataType_(2.234);
                (*data.f_6)[i] = DataType_(2.234);
                (*data.f_7)[i] = DataType_(2.234);
                (*data.f_8)[i] = DataType_(2.234);
                (*data.f_temp_1)[i] = DataType_(4711);
                (*data.f_temp_2)[i] = DataType_(4711);
                (*data.f_temp_3)[i] = DataType_(4711);
                (*data.f_temp_4)[i] = DataType_(4711);
                (*data.f_temp_5)[i] = DataType_(4711);
                (*data.f_temp_6)[i] = DataType_(4711);
                (*data.f_temp_7)[i] = DataType_(4711);
                (*data.f_temp_8)[i] = DataType_(4711);
            }
            CollideStreamGrid<Tag_, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>::
                value(info, data, tau);
            data.f_temp_1->lock(lm_read_only);
            //std::cout<<*data.f_temp_1;
            data.f_temp_1->unlock(lm_read_only);

            info.destroy();
            data.destroy();
            grid.destroy();
        }
};
CollideStreamGridLABSWETest<tags::CPU, float> collidestream_grid_test_float("float");
CollideStreamGridLABSWETest<tags::CPU, double> collidestream_grid_test_double("double");
CollideStreamGridLABSWETest<tags::CPU::Generic, float> generic_collidestream_grid_test_float("float");
CollideStreamGridLABSWETest<tags::CPU::Generic, double> generic_collidestream_grid_test_double("double");
#ifdef HONEI_SSE
CollideStreamGridLABSWETest<tags::CPU::SSE, float> sse_collidestream_grid_test_float("float");
CollideStreamGridLABSWETest<tags::CPU::SSE, double> sse_collidestream_grid_test_double("double");
#endif
#ifdef HONEI_CUDA
CollideStreamGridLABSWETest<tags::GPU::CUDA, float> cuda_collidestream_grid_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
CollideStreamGridLABSWETest<tags::GPU::CUDA, double> cuda_collidestream_grid_test_double("double");
#endif
#endif
#ifdef HONEI_CELL
CollideStreamGridLABSWETest<tags::Cell, float> cell_collidestream_grid_test_float("float");
#endif
#ifdef HONEI_OPENCL
CollideStreamGridLABSWETest<tags::OpenCL::CPU, float> ocl_cpu_collidestream_grid_test_float("float");
CollideStreamGridLABSWETest<tags::OpenCL::CPU, double> ocl_cpu_collidestream_grid_test_double("double");
#endif
#ifdef HONEI_OPENCL_GPU
CollideStreamGridLABSWETest<tags::OpenCL::GPU, float> ocl_gpu_collidestream_grid_test_float("float");
CollideStreamGridLABSWETest<tags::OpenCL::GPU, double> ocl_gpu_collidestream_grid_test_double("double");
#endif
#ifdef HONEI_OPENCL_ACC
CollideStreamGridLABSWETest<tags::OpenCL::Accelerator, float> ocl_acc_collidestream_grid_test_float("float");
CollideStreamGridLABSWETest<tags::OpenCL::Accelerator, double> ocl_acc_collidestream_grid_test_double("double");
#endif
