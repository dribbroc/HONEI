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
#include <honei/lbm/scan_conversion_fsi.hh>
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
using namespace lbm::lbm_boundary_types;

template <typename Tag_, typename DataType_>
class ScanConversionFSITest :
    public TaggedTest<Tag_>
{
    public:
        ScanConversionFSITest(const std::string & type) :
            TaggedTest<Tag_>("scan_conversion_fsi_test<" + type + ">")
        {
        }

        virtual void run() const
        {

            unsigned long g_h(50);
            unsigned long g_w(50);

            Grid<D2Q9, DataType_> grid;

            ScenarioCollection::get_scenario(0, g_h, g_w, grid);

            PackedGridData<D2Q9, DataType_>  data;
            PackedSolidData<D2Q9, DataType_>  solids;
            PackedGridInfo<D2Q9> info;

            DenseMatrix<bool> line(g_h, g_w, false);
            DenseMatrix<bool> bound(g_h, g_w, false);
            DenseMatrix<bool> stf(g_h, g_w, false);
            DenseMatrix<bool> sol(g_h, g_w, false);

            GridPacker<D2Q9, NOSLIP, DataType_>::pack(grid, info, data);
            GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DataType_>::cuda_pack(info, data);
            GridPackerFSI<D2Q9, NOSLIP, DataType_>::allocate(data, solids);
            GridPackerFSI<D2Q9, NOSLIP, DataType_>::pack(grid, data, solids, line, bound, stf, sol, *grid.obstacles);

            //Directly dealing with omega-coordinates
            Line<DataType_, lbm_solid_dims::D2> line_1_0(DataType_(5) * grid.d_x, DataType_(20) * grid.d_y, DataType_(10)* grid.d_x, DataType_(20) * grid.d_y);
            Line<DataType_, lbm_solid_dims::D2> line_2_0(DataType_(10)* grid.d_x, DataType_(20) * grid.d_y, DataType_(10)* grid.d_x, DataType_(25) * grid.d_y);
            Line<DataType_, lbm_solid_dims::D2> line_3_0(DataType_(10)* grid.d_x, DataType_(25) * grid.d_y, DataType_(5)* grid.d_x, DataType_(25) * grid.d_y);
            Line<DataType_, lbm_solid_dims::D2> line_4_0(DataType_(5)* grid.d_x, DataType_(25) * grid.d_y, DataType_(5)* grid.d_x, DataType_(20) * grid.d_y);

            Polygon<DataType_, lbm_solid_dims::D2> tri_0(4);
            tri_0.add_line(line_1_0);
            tri_0.add_line(line_2_0);
            tri_0.add_line(line_3_0);
            tri_0.add_line(line_4_0);
            tri_0.value();

            ScanConversionFSI<Tag_>::value(grid, info, data, solids, tri_0, true);
            DenseMatrix<bool> line_res(g_h, g_w, false);
            GridPackerFSI<D2Q9, NOSLIP, DataType_>::deflate(grid, solids, solids.line_flags, &line_res);
            DenseMatrix<bool> boundary_res(g_h, g_w, false);
            GridPackerFSI<D2Q9, NOSLIP, DataType_>::deflate(grid, solids, solids.boundary_flags, &boundary_res);
            DenseMatrix<bool> solid_res(g_h, g_w, false);
            GridPackerFSI<D2Q9, NOSLIP, DataType_>::deflate(grid, solids, solids.solid_flags, &solid_res);

            std::cout << line_res << std::endl;
            std::cout << boundary_res << std::endl;
            std::cout << solid_res << std::endl;
        }


};
ScanConversionFSITest<tags::CPU, float> solver_test_float("float");
