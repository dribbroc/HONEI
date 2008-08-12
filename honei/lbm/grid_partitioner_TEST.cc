/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
#include <honei/lbm/grid_partitioner.hh>
#include <honei/lbm/grid_packer.hh>
#include <unittest/unittest.hh>
#include <iostream>

using namespace honei;
using namespace tests;
using namespace std;
using namespace lbm;
using namespace lbm_lattice_types;

template <typename Tag_, typename DataType_>
class GridPartitionerTest :
    public TaggedTest<Tag_>
{
    public:
        GridPartitionerTest(const std::string & type) :
            TaggedTest<Tag_>("grid_partitioner_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            DenseMatrix<DataType_> dummy(10, 12, DataType_(2));
            dummy(2, 5) = DataType_(5);
            DenseMatrix<bool> obst(10, 12, false);
            /*
             +++++++++++*
             +++++++++++*
             ++++*+++++++
             +++X++++++++
             */
            //obst(0,5) = true;
            obst(0, 11) = true;
            obst(1, 11) = true;
            obst(2, 4) = true;

            /*obst(0,4) = true;
            obst(1,5) = true;
            obst(1,6) = true;
            obst(0,6) = true;
            for (MutableMatrix<bool>::ElementIterator i(obst.begin_elements()) ; i.index() < 62 ; ++i)
            {
                *i = true;
            }*/
            std::cout<<obst;
            PackedGridInfo<D2Q9> info;
            PackedGridData<D2Q9, DataType_> data;

            Grid<D2Q9, DataType_> grid;

            grid.h = new DenseMatrix<DataType_>(dummy.copy());
            grid.u = new DenseMatrix<DataType_>(dummy.copy());
            grid.v = new DenseMatrix<DataType_>(dummy.copy());
            grid.obstacles = new DenseMatrix<bool>(obst);

            std::list<PackedGridInfo<D2Q9> > info_list;
            std::list<PackedGridData<D2Q9, DataType_> > data_list;

            GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DataType_>::pack(grid, info, data);
            GridPartitioner<D2Q9, DataType_>::partition(2, info, data, info_list, data_list);

            TEST_CHECK(true);
        }
};

GridPartitionerTest<tags::CPU, float> gptest_float("float");

