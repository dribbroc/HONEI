/* vim: set number sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Volker Jung <volker.m.jung@t-online.de>
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <limiter.hh>
#include <unittest/unittest.hh>

#include <string>

using namespace pg512;
using namespace tests;

template <typename DataType_>
class MinModLimiterQuickTest :
    public QuickTest
{
    public:
        MinModLimiterQuickTest(const std::string & type) :
            QuickTest("mm_limiter_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            DataType_ s1(-1.23456);
            DataType_ s2(1.23456);
            DataType_ s3(0.23456);
            min_mod_limiter(s1);
            TEST_CHECK_EQUAL_WITHIN_EPS(s1, DataType_(0), std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(min_mod_limiter(s2), DataType_(1), std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(min_mod_limiter(s3), DataType_(0.23456), std::numeric_limits<DataType_>::epsilon());
        }
};
MinModLimiterQuickTest<float> mm_limiter_quick_test_float("float");
MinModLimiterQuickTest<double> mm_limiter_quick_test_double("double");

template <typename DataType_>
class SuperBeeLimiterQuickTest :
    public QuickTest
{
    public:
        SuperBeeLimiterQuickTest(const std::string & type) :
            QuickTest("sb_limiter_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            DataType_ s1(-1.23456);
            DataType_ s2(2.23456);
            DataType_ s3(0.23456);
            DataType_ s4(0.53456);
            DataType_ s5(1.53456);
            TEST_CHECK_EQUAL_WITHIN_EPS(super_bee_limiter(s1), DataType_(0), std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(super_bee_limiter(s2), DataType_(2), std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(super_bee_limiter(s3), DataType_(2*0.23456), std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(super_bee_limiter(s4), DataType_(1), std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(super_bee_limiter(s5), DataType_(1.53456), std::numeric_limits<DataType_>::epsilon());
        }
};
SuperBeeLimiterQuickTest<float> sb_limiter_quick_test_float("float");
SuperBeeLimiterQuickTest<double> sb_limiter_quick_test_double("double");

template <typename DataType_>
class MonotonizedCentralLimiterQuickTest :
    public QuickTest
{
    public:
        MonotonizedCentralLimiterQuickTest(const std::string & type) :
            QuickTest("mc_limiter_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            DataType_ s1(-1.23456);
            DataType_ s2(0.33456);
            DataType_ s3(0.23456);
            DataType_ s4(3.13456);
            TEST_CHECK_EQUAL_WITHIN_EPS(monotonized_central_limiter(s1), DataType_(0), std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(monotonized_central_limiter(s2), DataType_((1 + s2) / 2), std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(monotonized_central_limiter(s3), DataType_(2 * 0.23456), std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(monotonized_central_limiter(s4), DataType_(2), std::numeric_limits<DataType_>::epsilon());
        }
};
MonotonizedCentralLimiterQuickTest<float> mc_limiter_quick_test_float("float");
MonotonizedCentralLimiterQuickTest<double> mc_limiter_quick_test_double("double");

template <typename DataType_>
class VanLeerLimiterQuickTest :
    public QuickTest
{
    public:
        VanLeerLimiterQuickTest(const std::string & type) :
            QuickTest("vl_limiter_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            DataType_ s1(-1.23456);
            DataType_ s2(3.33456);
            DataType_ s3(0.23456);
            DataType_ s4(999999.13456);
            TEST_CHECK_EQUAL_WITHIN_EPS(van_leer_limiter(s1), DataType_(0), std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(van_leer_limiter(s2), DataType_((s2 + s2) / (s2 + 1)), std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(van_leer_limiter(s3), DataType_((s3 + s3) / (s3 + 1)), std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK(van_leer_limiter(s4) < DataType_(2));
        }
};
VanLeerLimiterQuickTest<float> vl_limiter_quick_test_float("float");
VanLeerLimiterQuickTest<double> vl_limiter_quick_test_double("double");
