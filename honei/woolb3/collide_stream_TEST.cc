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
#include <honei/woolb3/grid3.hh>
#include <honei/woolb3/packed_grid3.hh>
#include <honei/woolb3/collide_stream.hh>


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
            unsigned long g_h(16);
            unsigned long g_w(16);

            DenseMatrix<DataType_> h(g_h, g_w, DataType_(0.05));

            DenseMatrix<DataType_> u(g_h, g_w, DataType_(0.));
            DenseMatrix<DataType_> v(g_h, g_w, DataType_(0.));

            DenseMatrix<DataType_> b(g_h, g_w, DataType_(0.));

            DenseMatrix<bool> geometry(g_h, g_w, false);
            Grid3<DataType_, 9> grid(geometry, h, b, u, v);
            PackedGrid3<DataType_, 9> pgrid(grid);

            DataType_ tau (0.5);

            for(unsigned long i(0); i < pgrid.h->size(); i++)
            {
                for (unsigned long j(0) ; j < 9 ; ++j)
                {
                    (*pgrid.f_eq[j])[i] = DataType_(1.234);
                    (*pgrid.f[j])[i] = DataType_(11.234);
                    (*pgrid.f_temp[j])[i] = DataType_(4711);
                }
            }
            CollideStream<Tag_>::value(pgrid, tau, false);
            CollideStream<Tag_>::value(pgrid, tau, true);
            //std::cout<< *pgrid.f_temp[2];
        }

};
CollideStreamTest<tags::CPU, float> collide_stream_test_float("float");
