/* vim: set number sw=4 sts=4 et foldmethod=syntax : */

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
class MM_LimiterTest :
    public BaseTest
{
    public:
        MM_LimiterTest(const std::string & type) :
            BaseTest("mm_limiter_test<" + type + ">")
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

MM_LimiterTest<float> mm_limiter_test_float("float");
MM_LimiterTest<double> mm_limiter_test_double("double");


// Testcases for class SB_Limiter

template <typename DataType_>
class SB_LimiterTest :
    public BaseTest
{
    public:
        SB_LimiterTest(const std::string & type) :
            BaseTest("sb_test<" + type + ">")
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

SB_LimiterTest<float> sb_limiter_test_float("float");
SB_LimiterTest<double> sb_limiter_test_double("double");

// Testcases for class MC_Limiter

template <typename DataType_>
class MC_LimiterTest :
    public BaseTest
{
    public:
        MC_LimiterTest(const std::string & type) :
            BaseTest("mc_test<" + type + ">")
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

MC_LimiterTest<float> mc_limiter_test_float("float");
MC_LimiterTest<double> mc_limiter_test_double("double");

// Testcases for class VL_Limiter

template <typename DataType_>
class VL_LimiterTest :
    public BaseTest
{
    public:
        VL_LimiterTest(const std::string & type) :
            BaseTest("vl_test<" + type + ">")
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

VL_LimiterTest<float> vl_limiter_test_float("float");
VL_LimiterTest<double> vl_limiter_test_double("double");
