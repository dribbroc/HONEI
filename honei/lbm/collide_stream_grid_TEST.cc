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
#include <unittest/unittest.hh>
#include <honei/lbm/grid.hh>
#include <honei/lbm/collide_stream_grid.hh>

using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DataType_>
class CollideStreamLABSWETest :
    public TaggedTest<Tag_>
{
    public:
        CollideStreamLABSWETest(const std::string & type) :
            TaggedTest<Tag_>("collideandstream_labswe_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            PackedGridData<lbm_lattice_types::D2Q9, DataType_> data;
            PackedGridInfo<lbm_lattice_types::D2Q9> info;

            DenseVector<DataType_> e_x(9ul, DataType_(2.));
            DenseVector<DataType_> e_y(9ul, DataType_(2.));
            DenseVector<DataType_> s_x(1000ul, DataType_(2.));
            DenseVector<DataType_> s_y(1000ul, DataType_(2.));
            data.f_eq_0 = new DenseVector<DataType_>(1000);
            data.f_eq_1 = new DenseVector<DataType_>(1000);
            data.f_eq_2 = new DenseVector<DataType_>(1000);
            data.f_eq_3 = new DenseVector<DataType_>(1000);
            data.f_eq_4 = new DenseVector<DataType_>(1000);
            data.f_eq_5 = new DenseVector<DataType_>(1000);
            data.f_eq_6 = new DenseVector<DataType_>(1000);
            data.f_eq_7 = new DenseVector<DataType_>(1000);
            data.f_eq_8 = new DenseVector<DataType_>(1000);
            data.f_0 = new DenseVector<DataType_>(1000);
            data.f_1 = new DenseVector<DataType_>(1000);
            data.f_2 = new DenseVector<DataType_>(1000);
            data.f_3 = new DenseVector<DataType_>(1000);
            data.f_4 = new DenseVector<DataType_>(1000);
            data.f_5 = new DenseVector<DataType_>(1000);
            data.f_6 = new DenseVector<DataType_>(1000);
            data.f_7 = new DenseVector<DataType_>(1000);
            data.f_8 = new DenseVector<DataType_>(1000);
            data.f_temp_0 = new DenseVector<DataType_>(1000);
            data.f_temp_1 = new DenseVector<DataType_>(1000);
            data.f_temp_2 = new DenseVector<DataType_>(1000);
            data.f_temp_3 = new DenseVector<DataType_>(1000);
            data.f_temp_4 = new DenseVector<DataType_>(1000);
            data.f_temp_5 = new DenseVector<DataType_>(1000);
            data.f_temp_6 = new DenseVector<DataType_>(1000);
            data.f_temp_7 = new DenseVector<DataType_>(1000);
            data.f_temp_8 = new DenseVector<DataType_>(1000);
            DataType_ tau (1);
            info.limits = new DenseVector<unsigned long>(3);
            (*info.limits)[0] = 0;
            (*info.limits)[1] = 10;
            (*info.limits)[2] = 1000;

            info.dir_0 = new DenseVector<unsigned long>(3);
            (*info.dir_0)[0] = 0;
            (*info.dir_0)[1] = 10;
            (*info.dir_0)[2] = 100;
            info.dir_1 = new DenseVector<unsigned long>(3);
            (*info.dir_1)[0] = 0;
            (*info.dir_1)[1] = 0;
            (*info.dir_1)[2] = 0;
            info.dir_2 = new DenseVector<unsigned long>(3);
            (*info.dir_2)[0] = 0;
            (*info.dir_2)[1] = 0;
            (*info.dir_2)[2] = 0;
            info.dir_3 = new DenseVector<unsigned long>(3);
            (*info.dir_3)[0] = 0;
            (*info.dir_3)[1] = 0;
            (*info.dir_3)[2] = 0;
            info.dir_4 = new DenseVector<unsigned long>(3);
            (*info.dir_4)[0] = 0;
            (*info.dir_4)[1] = 0;
            (*info.dir_4)[2] = 0;
            info.dir_5 = new DenseVector<unsigned long>(3);
            (*info.dir_5)[0] = 0;
            (*info.dir_5)[1] = 0;
            (*info.dir_5)[2] = 0;
            info.dir_6 = new DenseVector<unsigned long>(3);
            (*info.dir_6)[0] = 0;
            (*info.dir_6)[1] = 0;
            (*info.dir_6)[2] = 0;
            info.dir_7 = new DenseVector<unsigned long>(3);
            (*info.dir_7)[0] = 0;
            (*info.dir_7)[1] = 0;
            (*info.dir_7)[2] = 0;
            info.dir_8 = new DenseVector<unsigned long>(3);
            (*info.dir_8)[0] = 0;
            (*info.dir_8)[1] = 0;
            (*info.dir_8)[2] = 0;

            CollideStreamGrid<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>::
                value(info, data, s_x, s_y, e_x, e_y, tau);
            TEST_CHECK(true);
        }
};
CollideStreamLABSWETest<tags::CPU, float> source_test_float("float");
CollideStreamLABSWETest<tags::CPU, double> source_test_double("double");
CollideStreamLABSWETest<tags::CPU::MultiCore, float> source_test_float_mc("float");
CollideStreamLABSWETest<tags::CPU::MultiCore, double> source_test_double_mc("double");
#ifdef HONEI_SSE
CollideStreamLABSWETest<tags::CPU::SSE, float> source_test_float_sse("float");
CollideStreamLABSWETest<tags::CPU::SSE, double> source_test_double_sse("double");
CollideStreamLABSWETest<tags::CPU::MultiCore::SSE, float> source_test_float_mc_sse("float");
CollideStreamLABSWETest<tags::CPU::MultiCore::SSE, double> source_test_double_mc_sse("double");
#endif
#ifdef HONEI_CELL
CollideStreamLABSWETest<tags::Cell, float> source_test_float_cell("float");
CollideStreamLABSWETest<tags::Cell, double> source_test_double_cell("double");
#endif
