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
            DenseMatrix<DataType_> dummy(10, 12, DataType_(0));
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

            GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DataType_>::pack(grid, info, data);
            TEST_CHECK(true);

            /*for (unsigned long i(0) ; i < info.limits->size() ; ++i)
            {
                std::cout<<(*info.limits)[i] <<std::endl;
            }*/
            unsigned long observe(19);
            std::cout<<(*info.dir_0)[observe] << " O: " << (*info.dir_1)[observe] << " N: " << (*info.dir_3)[observe] << " W: " << (*info.dir_5)[observe] << " S: " << (*info.dir_7)[observe] <<std::endl;
            std::cout<<" NO: " << (*info.dir_2)[observe] << " NW: " << (*info.dir_4)[observe] << " SW: " << (*info.dir_6)[observe] << " SO: " << (*info.dir_8)[observe] <<std::endl;
            ++observe;
            std::cout<<(*info.dir_0)[observe] << " O: " << (*info.dir_1)[observe] << " N: " << (*info.dir_3)[observe] << " W: " << (*info.dir_5)[observe] << " S: " << (*info.dir_7)[observe] <<std::endl;
            std::cout<<" NO: " << (*info.dir_2)[observe] << " NW: " << (*info.dir_4)[observe] << " SW: " << (*info.dir_6)[observe] << " SO: " << (*info.dir_8)[observe] <<std::endl;
            ++observe;
            std::cout<<(*info.dir_0)[observe] << " O: " << (*info.dir_1)[observe] << " N: " << (*info.dir_3)[observe] << " W: " << (*info.dir_5)[observe] << " S: " << (*info.dir_7)[observe] <<std::endl;
            std::cout<<" NO: " << (*info.dir_2)[observe] << " NW: " << (*info.dir_4)[observe] << " SW: " << (*info.dir_6)[observe] << " SO: " << (*info.dir_8)[observe] <<std::endl;

        }
};

GridPackerTest<tags::CPU, float> gptest_float("float");

