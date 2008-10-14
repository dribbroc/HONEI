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
            unsigned long g_h(50);
            unsigned long g_w(50);

            DenseMatrix<DataType_> h(g_h, g_w, DataType_(0.05));

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

            //Other matrices needed by solver:
            /// \todo
            DenseVector<DataType_> b(data.h->size(), DataType_(0.));

            data.distribution_x = new DenseVector<DataType_>(9ul, DataType_(2.));
            data.distribution_y = new DenseVector<DataType_>(9ul, DataType_(2.));
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

            }
            CollideStreamGrid<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>::
                value(info, data, tau);
            TEST_CHECK(true);

            //std::cout<<GridPacker<D2Q9, NOSLIP, DataType_>::extract_ftemp2(grid, info, data)<<std::endl;
            //std::cout<<*info.limits<<std::endl;
            //std::cout<<*info.dir_2<<std::endl;

        }
};
//CollideStreamGridLABSWETest<tags::CPU, float> collidestream_grid_test_float("float");
CollideStreamGridLABSWETest<tags::CPU, double> collidestream_grid_test_double("double");
/*CollideStreamGridLABSWETest<tags::CPU::MultiCore, float> collidestream_grid_test_float_mc("float");
CollideStreamGridLABSWETest<tags::CPU::MultiCore, double> collidestream_grid_test_double_mc("double");
#ifdef HONEI_SSE
CollideStreamGridLABSWETest<tags::CPU::SSE, float> collidestream_grid_test_float_sse("float");
CollideStreamGridLABSWETest<tags::CPU::SSE, double> collidestream_grid_test_double_sse("double");
CollideStreamGridLABSWETest<tags::CPU::MultiCore::SSE, float> collidestream_grid_test_float_mc_sse("float");
CollideStreamGridLABSWETest<tags::CPU::MultiCore::SSE, double> collidestream_grid_test_double_mc_sse("double");
#endif
#ifdef HONEI_CELL
CollideStreamGridLABSWETest<tags::Cell, float> collidestream_grid_test_float_cell("float");
CollideStreamGridLABSWETest<tags::Cell, double> collidestream_grid_test_double_cell("double");
#endif
*/
