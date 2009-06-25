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
#include <honei/lbm/boundary_init_fsi.hh>
#include <honei/lbm/grid_packer.hh>

using namespace honei;
using namespace tests;
using namespace std;
using namespace lbm;
using namespace lbm_boundary_types;
using namespace lbm_lattice_types;

template <typename Tag_, typename DataType_>
class BoundaryInitFSITest :
    public TaggedTest<Tag_>
{
    public:
        BoundaryInitFSITest(const std::string & type) :
            TaggedTest<Tag_>("boundary_init_test (FSI) <" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long g_h(12);
            unsigned long g_w(12);

            DenseMatrix<DataType_> h(g_h, g_w, DataType_(0.05));
            h[5][1] = DataType_(0.4);
            h[5][2] = DataType_(0.3);
            h[5][3] = DataType_(0.2);

            DenseMatrix<DataType_> u(g_h, g_w, DataType_(0.));
            u[5][2] = DataType_(0);
            u[5][3] = DataType_(0);
            DenseMatrix<DataType_> v(g_h, g_w, DataType_(0.));
            v[5][2] = DataType_(0);
            v[5][3] = DataType_(0);

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

            std::vector<unsigned long> lines_i;
            std::vector<unsigned long> lines_j;
            PackedSolidData<D2Q9, DataType_> solids(lines_i, lines_j);

            DenseMatrix<bool> line(g_h, g_w, false);
            line[5][5] = true;
            DenseMatrix<bool> bound(g_h, g_w, false);
            bound[4][4] = true;
            bound[4][5] = true;
            bound[4][6] = true;
            bound[5][4] = true;
            bound[5][6] = true;
            bound[6][4] = true;
            bound[6][5] = true;
            bound[6][6] = true;

            DenseMatrix<bool> stf(g_h, g_w, false);
            stf[5][4] = true;

            DenseMatrix<bool> sol(g_h, g_w, false);
            sol[5][5] = true;

            GridPackerFSI<D2Q9, NOSLIP, DataType_>::allocate(data, solids);
            GridPackerFSI<D2Q9, NOSLIP, DataType_>::pack(grid, data, solids, line, bound, sol, stf, obstacles_ref);

            grid.d_x = DataType_(0.1);
            grid.d_y = DataType_(0.1);
            solids.current_u = DataType_(1);
            solids.current_v = DataType_(1);

            BoundaryInitFSI<Tag_, D2Q9::DIR_1>::value(grid, info, data, solids);

            //in matrix-form:
            DenseMatrix<DataType_> res_h(h.rows(), h.columns());
            DenseMatrix<DataType_> res_u(h.rows(), h.columns());
            DenseMatrix<DataType_> res_v(h.rows(), h.columns());

            GridPacker<D2Q9, NOSLIP, DataType_>::deflate(grid, data, data.h, &res_h);
            GridPacker<D2Q9, NOSLIP, DataType_>::deflate(grid, data, data.u, &res_u);
            GridPacker<D2Q9, NOSLIP, DataType_>::deflate(grid, data, data.v, &res_v);

            std::cout << res_h << std::endl;
            std::cout << res_u << std::endl;
            std::cout << res_v << std::endl;

        }
};
BoundaryInitFSITest<tags::CPU, float> collidestream_grid_test_float("float");
