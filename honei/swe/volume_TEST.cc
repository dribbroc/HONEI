/* vim: set number sw=4 sts=4 et foldmethod=syntax : */

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

#include <volume.hh>
#include <honei/util/unittest.hh>
#include <iostream>
#include <string>

using namespace honei;
using namespace tests;

template <typename DataType_>
class VolumeQuickTest :
    public QuickTest
{
    public:
        VolumeQuickTest(const std::string & type) :
            QuickTest("volume_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            DenseMatrix<DataType_> height(30 ,30, DataType_(1));
            Cylinder<DataType_> c1(height, DataType_(5), 0, 15);
            Cylinder<DataType_> c2(height, DataType_(5), 15, 0);
            Cylinder<DataType_> c3(height, DataType_(5), 29, 15);
            Cylinder<DataType_> c4(height, DataType_(5), 15, 29);

            c1.value();
            c2.value();
            c3.value();
            c4.value();

            std::cout << height << std::endl;

            DenseMatrix<DataType_> height2(30 ,30, DataType_(1));
            Cuboid<DataType_> q1(height2, 15, 10, DataType_(5), 4, 4);
            q1.value();

            Cuboid<DataType_> q2(height2, 15, 10, DataType_(5), 28, 4);
            q2.value();

            std::cout << height2 << std::endl;

            DenseMatrix<DataType_> height3(30 ,30, DataType_(1));
            Cuboid<DataType_> q3(height3, 15, 10, DataType_(6), 5, 5);
            Cylinder<DataType_> c5(height3, DataType_(5), 15, 29);

            VolumeList v;
            v.insert(&q3);
            v.insert(&c5);
            v.convex_hull();

            std::cout << height3 << std::endl;

            TEST_CHECK(true);
        }
};
VolumeQuickTest<float> ws_quick_test_float("float");
VolumeQuickTest<double> ws_quick_test_double("double");
