/* vim: set number sw=4 sts=4 et foldmethod=syntax : */

#include <limiter.hh>
#include <unittest/unittest.hh>

#include <string>

using namespace pg512;
using namespace tests;

/*** ***********************************************************
 * 1. Testcases for each derived limiter class 
 **************************************************************/

// Testcases for class MM_Limiter

template <typename Tag_ = tags::CPU>
class MM_LimiterCreationTest :
    public BaseTest
{
    public:
        MM_LimiterCreationTest(const std::string & type) :
            BaseTest("mm_limiter_creation_test<" + type + ">")
        {
         }

        virtual void run() const
        {
            std::tr1::shared_ptr<MM_Limiter<> > mm_limiter(new MM_Limiter<>() );
            TEST_CHECK(true);
         }
}; 

MM_LimiterCreationTest<float> mm_limiter_test_float("float");
MM_LimiterCreationTest<double> mm_limiter_test_double("double");
MM_LimiterCreationTest<int> mm_limiter_creation_test("int");


// Testcases for class SB_Limiter

template <typename Tag_ = tags::CPU>
class SB_LimiterCreationTest :
    public BaseTest
{ 
    public:
        SB_LimiterCreationTest(const std::string & type) :
            BaseTest("sb_limiter_creation_test<" + type + ">")
         {
        }

        virtual void run() const
         {
            std::tr1::shared_ptr<SB_Limiter<> > sb_limiter(new SB_Limiter<>() );
            TEST_CHECK(true);
        }
};

SB_LimiterCreationTest<float> sb_limiter_test_float("float");
SB_LimiterCreationTest<double> sb_limiter_test_double("double");
SB_LimiterCreationTest<int> sb_limiter_creation_test("int");



// Testcases for class MC_Limiter

template <typename Tag_ = tags::CPU>
class MC_LimiterCreationTest :
    public BaseTest
{
    public:
        MC_LimiterCreationTest(const std::string & type) :
            BaseTest("mc_limiter_creation_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            std::tr1::shared_ptr<MC_Limiter<> > mc_limiter(new MC_Limiter<>() );
            TEST_CHECK(true);
         }
};

MC_LimiterCreationTest<float> mc_limiter_test_float("float");
MC_LimiterCreationTest<double> mc_limiter_test_double("double");
MC_LimiterCreationTest<int> mc_limiter_creation_test("int");



// Testcases for class VL_Limiter

template <typename Tag_ = tags::CPU>
class VL_LimiterCreationTest :
    public BaseTest
{
    public:
        VL_LimiterCreationTest(const std::string & type) :
            BaseTest("VL_limiter_creation_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            std::tr1::shared_ptr<VL_Limiter<> > vl_limiter(new VL_Limiter<>() );
            TEST_CHECK(true);
         }
};

VL_LimiterCreationTest<float> vl_limiter_test_float("float");
VL_LimiterCreationTest<double> vl_limiter_test_double("double");
VL_LimiterCreationTest<int> vl_limiter_creation_test("int");
