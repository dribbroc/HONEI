/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2011 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/util/unittest.hh>
#include <iostream>
#include <honei/woolb3/grid.hh>
#include <honei/woolb3/packed_grid.hh>
#include <honei/woolb3/collide_stream.hh>
#include <honei/util/time_stamp.hh>


using namespace honei;
using namespace tests;


template <typename Tag_, typename DataType_>
class CollideStreamTest :
    public TaggedTest<Tag_>
{
    public:
        CollideStreamTest(const std::string & type) :
            TaggedTest<Tag_>("collide_stream_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            //unsigned long size(2048);
            unsigned long size(256);
            DenseMatrix<bool> geometry(size, size, false);
            DenseMatrix<DataType_> h(size, size, 1);
            DenseMatrix<DataType_> b(size, size, 1);
            DenseMatrix<DataType_> u(size, size, 1);
            DenseMatrix<DataType_> v(size, size, 1);

            Grid<DataType_, 9> grid(geometry, h, b, u, v);
            PackedGrid<DataType_, 9> pgrid(grid);

            TimeStamp at, bt;
            at.take();
            CollideStream<Tag_>::value(pgrid, DataType_(1.23456));
            bt.take();
            std::cout<<"TOE: "<<bt.total()-at.total()<<std::endl;
        }

};
CollideStreamTest<tags::CPU, float> collide_stream_test_float("float");
