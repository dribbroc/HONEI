/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <limiter.hh>
#include <unittest/unittest.hh>

#include <string>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
class LimiterTest :
    public BaseTest
{
    public:
        LimiterTest(const std::string & type) :
            BaseTest("limiter_test<" + type + ">")
        {
        }

        virtual void run() const
        {
                TEST_CHECK(true);
            }
        }
};

LimiterTest<float> limiter_test_float("float");
LimiterTest<double> limiter_test_double("double");
