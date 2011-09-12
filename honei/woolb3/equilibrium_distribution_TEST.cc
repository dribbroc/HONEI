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
#include <honei/woolb3/equilibrium_distribution.hh>


using namespace honei;
using namespace tests;


template <typename Tag_, typename DataType_>
class EquilibriumDistributionTest :
    public TaggedTest<Tag_>
{
    public:
        EquilibriumDistributionTest(const std::string & type) :
            TaggedTest<Tag_>("equilibrium_distribution_test<" + type + ">")
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
            Grid<DataType_, 9> grid(geometry, h, b, u, v);
            PackedGrid<DataType_, 9> pgrid(grid);

            DataType_ e (0.5);
            DataType_ g (9.81);

            EquilibriumDistribution<Tag_>::value(pgrid, e, g);
            //std::cout<< *pgrid.f_temp[2];
        }

};
EquilibriumDistributionTest<tags::CPU, float> equilibrium_distribution_test_float("float");
