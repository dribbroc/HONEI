/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the Math C++ library. LibMath is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibMath is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <libmath/interpolation.hh>
#include <unittest/unittest.hh>
#include <libutil/stringify.hh>
#include <iostream>

using namespace honei;
using namespace tests;
using namespace std;
using namespace interpolation_methods;
template <typename Tag_, typename DT1_>
class InterpolationTest:
    public BaseTest
{
    public:
        InterpolationTest(const std::string & tag) :
            BaseTest("Interpolation test<" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            DenseMatrix<DT1_> height(10,10, DT1_(0));
            height[5][5] = DT1_(20);
            height[6][5] = DT1_(20);
            height[6][6] = DT1_(20);
            height[5][6] = DT1_(20);

            DT1_ delta_x(1);
            DT1_ delta_y(1);
            DT1_ result = Interpolation<Tag_, LINEAR>::value(delta_x, delta_y, height, DT1_(5.1234), DT1_(5.00001));
            TEST_CHECK_EQUAL_WITHIN_EPS(DT1_(20), result, std::numeric_limits<DT1_>::epsilon());

            DT1_ result2 = Interpolation<Tag_, NN>::value(delta_x, delta_y, height, DT1_(5.1234), DT1_(5.00001));
            TEST_CHECK_EQUAL_WITHIN_EPS(DT1_(20), result, std::numeric_limits<DT1_>::epsilon());
        }
};
InterpolationTest<tags::CPU, double> interpolation_test_double("double");
InterpolationTest<tags::CPU, float> interpolation_test_float("float");
#ifdef HONEI_SSE
InterpolationTest<tags::CPU::SSE, double> interpolation_test_double_sse("SSE double");
InterpolationTest<tags::CPU::SSE, float> interpolation_test_float_sse("SSE float");
#endif
