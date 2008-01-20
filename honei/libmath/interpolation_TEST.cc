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

#include <honei/libmath/interpolation.hh>
#include <unittest/unittest.hh>
#include <honei/libutil/stringify.hh>
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

template <typename Tag_, typename DT1_>
class InterpolationBoundaryTest:
    public BaseTest
{
    public:
        InterpolationBoundaryTest(const std::string & tag) :
            BaseTest("Interpolation boundary test<" + tag + ">")
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
            DT1_ result = Interpolation<Tag_, LINEAR>::value(delta_x, delta_y, height, DT1_(345), DT1_(123));
            TEST_CHECK_EQUAL_WITHIN_EPS(DT1_(0), result, std::numeric_limits<DT1_>::epsilon());

            DT1_ result2 = Interpolation<Tag_, NN>::value(delta_x, delta_y, height, DT1_(423), DT1_(12212));
            TEST_CHECK_EQUAL_WITHIN_EPS(DT1_(0), result, std::numeric_limits<DT1_>::epsilon());
        }
};

template <typename Tag_, typename DT1_>
class InterpolationRimTest:
    public BaseTest
{
    public:
        InterpolationRimTest(const std::string & tag) :
            BaseTest("Interpolation at jump rim test<" + tag + ">")
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
            DT1_ result = Interpolation<Tag_, LINEAR>::value(delta_x, delta_y, height, DT1_(6.5), DT1_(5.5));

            DT1_ result2 = Interpolation<Tag_, NN>::value(delta_x, delta_y, height, DT1_(6.1), DT1_(6.1));
            TEST_CHECK_EQUAL_WITHIN_EPS(result, DT1_(10), std::numeric_limits<DT1_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(result2, DT1_(20), std::numeric_limits<DT1_>::epsilon());
        }
};

InterpolationTest<tags::CPU, double> interpolation_test_double("double");
InterpolationTest<tags::CPU, float> interpolation_test_float("float");
#ifdef HONEI_SSE
InterpolationTest<tags::CPU::SSE, double> interpolation_test_double_sse("SSE double");
InterpolationTest<tags::CPU::SSE, float> interpolation_test_float_sse("SSE float");
#endif

InterpolationBoundaryTest<tags::CPU, double> interpolation_b_test_double("double");
InterpolationBoundaryTest<tags::CPU, float> interpolation_b_test_float("float");
#ifdef HONEI_SSE
InterpolationBoundaryTest<tags::CPU::SSE, double> interpolation_b_test_double_sse("SSE double");
InterpolationBoundaryTest<tags::CPU::SSE, float> interpolation_b_test_float_sse("SSE float");
#endif

InterpolationRimTest<tags::CPU, double> interpolation_r_test_double("double");
InterpolationRimTest<tags::CPU, float> interpolation_r_test_float("float");
#ifdef HONEI_SSE
InterpolationRimTest<tags::CPU::SSE, double> interpolation_r_test_double_sse("SSE double");
InterpolationRimTest<tags::CPU::SSE, float> interpolation_r_test_float_sse("SSE float");
#endif
