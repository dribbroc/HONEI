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

/**************************************************************
 * Testcases for each derived limiter class 
 **************************************************************/

// Testcases for class MM_Limiter

template <typename DataType_>
class MM_LimiterQuickTest :
    public QuickTest
{
    public:
        MM_LimiterQuickTest(const std::string & type) :
            QuickTest("mm_limiter_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            DataType_ s1(-1.23456);
            DataType_ s2(1.23456);
            DataType_ s3(0.23456);
            TEST_CHECK_EQUAL(MM_Limiter<>::value(s1), static_cast<DataType_>(0));
            TEST_CHECK_EQUAL(MM_Limiter<>::value(s2), static_cast<DataType_>(1));
            TEST_CHECK_EQUAL(MM_Limiter<>::value(s3), static_cast<DataType_>(0.23456));           
        }
 }; 

MM_LimiterQuickTest<float> mm_limiter_quick_test_float("float");
MM_LimiterQuickTest<double> mm_limiter_quick_test_double("double");


// Testcases for class SB_Limiter

template <typename DataType_>
class SB_LimiterQuickTest :
    public QuickTest
{
    public:
        SB_LimiterQuickTest(const std::string & type) :
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
            TEST_CHECK_EQUAL(SB_Limiter<>::value(s1), static_cast<DataType_>(0));
            TEST_CHECK_EQUAL(SB_Limiter<>::value(s2), static_cast<DataType_>(2));
            TEST_CHECK_EQUAL(SB_Limiter<>::value(s3), static_cast<DataType_>(2*0.23456));
            TEST_CHECK_EQUAL(SB_Limiter<>::value(s4), static_cast<DataType_>(1));            
            TEST_CHECK_EQUAL(SB_Limiter<>::value(s5), static_cast<DataType_>(1.53456)); 
        }
 }; 

SB_LimiterQuickTest<float> sb_limiter_quick_test_float("float");
SB_LimiterQuickTest<double> sb_limiter_quick_test_double("double");

// Testcases for class MC_Limiter

template <typename DataType_>
class MC_LimiterQuickTest :
    public QuickTest
{
    public:
        MC_LimiterQuickTest(const std::string & type) :
            QuickTest("mc_limiter_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            DataType_ s1(-1.23456);
            DataType_ s2(0.33456);
            DataType_ s3(0.23456);
            DataType_ s4(3.13456);
            TEST_CHECK_EQUAL(MC_Limiter<>::value(s1), static_cast<DataType_>(0));
            TEST_CHECK_EQUAL(MC_Limiter<>::value(s2), static_cast<DataType_>((1+s2)/2));
            TEST_CHECK_EQUAL(MC_Limiter<>::value(s3), static_cast<DataType_>(2*0.23456));
            TEST_CHECK_EQUAL(MC_Limiter<>::value(s4), static_cast<DataType_>(2));            
        }
 }; 

MC_LimiterQuickTest<float> mc_limiter_quick_test_float("float");
MC_LimiterQuickTest<double> mc_limiter_quick_test_double("double");

// Testcases for class VL_Limiter

template <typename DataType_>
class VL_LimiterQuickTest :
    public QuickTest
{
    public:
        VL_LimiterQuickTest(const std::string & type) :
            QuickTest("vl_limiter_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            DataType_ s1(-1.23456);
            DataType_ s2(3.33456);
            DataType_ s3(0.23456);
            DataType_ s4(999999.13456);
            TEST_CHECK_EQUAL(VL_Limiter<>::value(s1), static_cast<DataType_>(0));
            TEST_CHECK_EQUAL(VL_Limiter<>::value(s2), static_cast<DataType_>((s2+s2)/(s2+1)));
            TEST_CHECK_EQUAL(VL_Limiter<>::value(s3), static_cast<DataType_>((s3+s3)/(s3+1)));
            TEST_CHECK(VL_Limiter<>::value(s4) < static_cast<DataType_>(2));            
        }
 }; 

VL_LimiterQuickTest<float> vl_limiter_quick_test_float("float");
VL_LimiterQuickTest<double> vl_limiter_quick_test_double("double");
