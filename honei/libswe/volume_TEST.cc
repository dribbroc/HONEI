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
#include <unittest/unittest.hh>
#include <iostream>
#include <string>

using namespace honei;
using namespace tests;
using namespace volume_types;

template <typename DataType_>
class VolumeQuickTest :
    public QuickTest
{
    public:
        VolumeQuickTest(const std::string & type) :
            QuickTest("water_shot_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            DenseMatrix<DataType_> height(30 ,30, DataType_(1));
            Volume<CYLINDRIC::STENCIL>::value(height, DataType_(5), 0, 15);
            Volume<CYLINDRIC::STENCIL>::value(height, DataType_(5), 15, 0);
            Volume<CYLINDRIC::STENCIL>::value(height, DataType_(5), 29, 15);
            Volume<CYLINDRIC::STENCIL>::value(height, DataType_(5), 15, 29);
            std::cout << height << std::endl;
            TEST_CHECK(true);
        }
};
VolumeQuickTest<float> ws_quick_test_float("float");
VolumeQuickTest<double> ws_quick_test_double("double");
