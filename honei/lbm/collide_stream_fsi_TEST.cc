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

#include <honei/lbm/tags.hh>
#include <unittest/unittest.hh>
#include <honei/lbm/grid.hh>
#include <honei/lbm/collide_stream_fsi.hh>
#include <honei/lbm/collide_stream_grid.hh>
#include <honei/lbm/update_velocity_directions_grid.hh>
#include <honei/lbm/grid_packer.hh>

using namespace honei;
using namespace tests;
using namespace std;
using namespace lbm;
using namespace lbm_boundary_types;

template <typename Tag_, typename DataType_>
class CollideStreamFSITest :
    public TaggedTest<Tag_>
{
    public:
        CollideStreamFSITest(const std::string & type) :
            TaggedTest<Tag_>("collide_and_stream_test (FSI) <" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(12);
            unsigned long g_w(12);

            DenseMatrix<DataType_> h(g_h, g_w, DataType_(0.05));

            DenseMatrix<DataType_> u(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> v(g_h, g_w, DataType_(0.));

            DenseMatrix<DataType_> b(g_h, g_w, DataType_(0.));

            Grid<D2Q9, DataType_> grid;
            DenseMatrix<bool> obstacles_ref(g_h, g_w, false);
            obstacles_ref[5][5] = true;
            grid.obstacles = &obstacles_ref;
            grid.h = &h;
            grid.u = &u;
            grid.v = &v;
            grid.b = &b;

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
            UpdateVelocityDirectionsGrid<Tag_, NOSLIP>::value(info, data);

            data.f_temp_1->lock(lm_read_only);
            DenseVector<DataType_> ref_result(data.f_temp_8->copy());
            data.f_temp_1->unlock(lm_read_only);

            //------------------------------------------------------------------------------------------
            Grid<D2Q9, DataType_> grid_2;
            DenseMatrix<bool> obstacles(g_h, g_w, false);
            grid_2.obstacles = &obstacles;
            grid_2.h = &h;
            grid_2.u = &u;
            grid_2.v = &v;
            grid_2.b = &b;

            PackedGridData<D2Q9, DataType_>  data_2;
            PackedGridInfo<D2Q9> info_2;

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid_2, info_2, data_2);
            GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DataType_>::cuda_pack(info_2, data_2);
            PackedSolidData<D2Q9, DataType_> solids;

            DenseMatrix<bool> line(g_h, g_w, false);
            line[5][5] = true;
            DenseMatrix<bool> bound(g_h, g_w, false);
            DenseMatrix<bool> stf(g_h, g_w, false);
            DenseMatrix<bool> sol(g_h, g_w, false);
            GridPackerFSI<D2Q9, NOSLIP, DataType_>::allocate(data_2, solids);
            GridPackerFSI<D2Q9, NOSLIP, DataType_>::pack(grid_2, data_2, solids, line, bound, sol, stf, obstacles);

            for(unsigned long i(0); i < data_2.h->size(); i++)
            {
                (*data_2.f_eq_0)[i] = DataType_(1.234);
                (*data_2.f_eq_1)[i] = DataType_(1.234);
                (*data_2.f_eq_2)[i] = DataType_(1.234);
                (*data_2.f_eq_3)[i] = DataType_(1.234);
                (*data_2.f_eq_4)[i] = DataType_(1.234);
                (*data_2.f_eq_5)[i] = DataType_(1.234);
                (*data_2.f_eq_6)[i] = DataType_(1.234);
                (*data_2.f_eq_7)[i] = DataType_(1.234);
                (*data_2.f_eq_8)[i] = DataType_(1.234);
                (*data_2.f_0)[i] = DataType_(2.234);
                (*data_2.f_1)[i] = DataType_(2.234);
                (*data_2.f_2)[i] = DataType_(11.234);
                (*data_2.f_3)[i] = DataType_(2.234);
                (*data_2.f_4)[i] = DataType_(2.234);
                (*data_2.f_5)[i] = DataType_(2.234);
                (*data_2.f_6)[i] = DataType_(2.234);
                (*data_2.f_7)[i] = DataType_(2.234);
                (*data_2.f_8)[i] = DataType_(2.234);
                (*data_2.f_temp_1)[i] = DataType_(4711);
                (*data_2.f_temp_2)[i] = DataType_(4711);
                (*data_2.f_temp_3)[i] = DataType_(4711);
                (*data_2.f_temp_4)[i] = DataType_(4711);
                (*data_2.f_temp_5)[i] = DataType_(4711);
                (*data_2.f_temp_6)[i] = DataType_(4711);
                (*data_2.f_temp_7)[i] = DataType_(4711);
                (*data_2.f_temp_8)[i] = DataType_(4711);
            }
            CollideStreamGrid<Tag_, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>::
                value(info_2, data_2, tau);
            CollideStreamFSI<Tag_, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>::
                value(info_2, data_2, solids, DataType_(0.), DataType_(0.));
            UpdateVelocityDirectionsGrid<Tag_, NOSLIP>::value(info_2, data_2);

            //in matrix-form:
            DenseMatrix<DataType_> ref(h.rows(), h.columns());
            DenseMatrix<DataType_> res(h.rows(), h.columns());

            GridPacker<D2Q9, NOSLIP, DataType_>::deflate(grid, data, data.f_temp_1 , &ref);
            GridPacker<D2Q9, NOSLIP, DataType_>::deflate(grid_2, data_2, data_2.f_temp_1, &res);

            for(unsigned long i(0) ; i < ref.rows() ; ++i)
            {
                for(unsigned long j(0) ; j < ref.columns() ; ++j)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(res[i][j], ref[i][j], std::numeric_limits<DataType_>::epsilon() * 10);
                }
            }
            GridPacker<D2Q9, NOSLIP, DataType_>::deflate(grid, data, data.f_temp_2 , &ref);
            GridPacker<D2Q9, NOSLIP, DataType_>::deflate(grid_2, data_2, data_2.f_temp_2, &res);

            for(unsigned long i(0) ; i < ref.rows() ; ++i)
            {
                for(unsigned long j(0) ; j < ref.columns() ; ++j)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(res[i][j], ref[i][j], std::numeric_limits<DataType_>::epsilon() * 10);
                }
            }

            GridPacker<D2Q9, NOSLIP, DataType_>::deflate(grid, data, data.f_temp_3 , &ref);
            GridPacker<D2Q9, NOSLIP, DataType_>::deflate(grid_2, data_2, data_2.f_temp_3, &res);

            for(unsigned long i(0) ; i < ref.rows() ; ++i)
            {
                for(unsigned long j(0) ; j < ref.columns() ; ++j)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(res[i][j], ref[i][j], std::numeric_limits<DataType_>::epsilon() * 10);
                }
            }

            GridPacker<D2Q9, NOSLIP, DataType_>::deflate(grid, data, data.f_temp_4 , &ref);
            GridPacker<D2Q9, NOSLIP, DataType_>::deflate(grid_2, data_2, data_2.f_temp_4, &res);

            for(unsigned long i(0) ; i < ref.rows() ; ++i)
            {
                for(unsigned long j(0) ; j < ref.columns() ; ++j)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(res[i][j], ref[i][j], std::numeric_limits<DataType_>::epsilon() * 10);
                }
            }

            GridPacker<D2Q9, NOSLIP, DataType_>::deflate(grid, data, data.f_temp_5 , &ref);
            GridPacker<D2Q9, NOSLIP, DataType_>::deflate(grid_2, data_2, data_2.f_temp_5, &res);

            for(unsigned long i(0) ; i < ref.rows() ; ++i)
            {
                for(unsigned long j(0) ; j < ref.columns() ; ++j)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(res[i][j], ref[i][j], std::numeric_limits<DataType_>::epsilon() * 10);
                }
            }

            GridPacker<D2Q9, NOSLIP, DataType_>::deflate(grid, data, data.f_temp_6 , &ref);
            GridPacker<D2Q9, NOSLIP, DataType_>::deflate(grid_2, data_2, data_2.f_temp_6, &res);

            for(unsigned long i(0) ; i < ref.rows() ; ++i)
            {
                for(unsigned long j(0) ; j < ref.columns() ; ++j)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(res[i][j], ref[i][j], std::numeric_limits<DataType_>::epsilon() * 10);
                }
            }

            GridPacker<D2Q9, NOSLIP, DataType_>::deflate(grid, data, data.f_temp_7 , &ref);
            GridPacker<D2Q9, NOSLIP, DataType_>::deflate(grid_2, data_2, data_2.f_temp_7, &res);

            for(unsigned long i(0) ; i < ref.rows() ; ++i)
            {
                for(unsigned long j(0) ; j < ref.columns() ; ++j)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(res[i][j], ref[i][j], std::numeric_limits<DataType_>::epsilon() * 10);
                }
            }

            GridPacker<D2Q9, NOSLIP, DataType_>::deflate(grid, data, data.f_temp_8 , &ref);
            GridPacker<D2Q9, NOSLIP, DataType_>::deflate(grid_2, data_2, data_2.f_temp_8, &res);

            for(unsigned long i(0) ; i < ref.rows() ; ++i)
            {
                for(unsigned long j(0) ; j < ref.columns() ; ++j)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(res[i][j], ref[i][j], std::numeric_limits<DataType_>::epsilon() * 10);
                }
            }

        }
};
CollideStreamFSITest<tags::CPU, float> collidestream_grid_test_float("float");
