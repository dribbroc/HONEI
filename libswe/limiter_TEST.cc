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

MM_LimiterCreationTest<> mm_limiter_test_float("");


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

SB_LimiterCreationTest<> sb_limiter_test_float("");


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

MC_LimiterCreationTest<> mc_limiter_test_float("");


// Testcases for class VL_Limiter

template <typename Tag_ = tags::CPU>
class VL_LimiterCreationTest :
    public BaseTest
{
    public:
        VL_LimiterCreationTest(const std::string & type) :
            BaseTest("vl_limiter_creation_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            std::tr1::shared_ptr<VL_Limiter<> > vl_limiter(new VL_Limiter<>() );
            TEST_CHECK(true);
        }
}; 

VL_LimiterCreationTest<> vl_limiter_test_float("");
