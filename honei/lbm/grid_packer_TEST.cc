/* vim: set sw=4 sts=4 et foldmethod=syntax : */

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
#include <honei/lbm/grid_packer.hh>
#include <unittest/unittest.hh>
#include <iostream>

using namespace honei;
using namespace tests;
using namespace std;
using namespace lbm;
using namespace lbm_lattice_types;

template <typename Tag_, typename DataType_>
class GridPackerTest :
    public TaggedTest<Tag_>
{
    public:
        GridPackerTest(const std::string & type) :
            TaggedTest<Tag_>("grid_packer_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            DenseMatrix<DataType_> dummy(10, 11, DataType_(0));
            DenseMatrix<bool> obst(10, 11, false);
            PackedGridInfo<D2Q9> info;
            PackedGridData<D2Q9, DataType_> data;
            GridPacker<D2Q9, DataType_> packer;

            Grid<D2Q9, DataType_> grid;

            grid.h = &dummy.copy();
            grid.u = &dummy.copy();
            grid.v = &dummy.copy();
            grid.obstacles = &obst;

            packer.pack(grid, info, data);
            TEST_CHECK(true);
        }
};

GridPackerTest<tags::CPU, float> gptest_float("float");

