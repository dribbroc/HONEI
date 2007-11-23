/* vim: set sw=4 sts=4 et foldmethod=syntax : */

// This test case needs debug support.
#ifndef DEBUG
#define DEBUG 1
#endif

#include <libutil/time_stamp.hh>
#include <unittest/unittest.hh>

using namespace honei;
using namespace tests;

class TimeStampTest :
    public QuickTest
{
    public:
        TimeStampTest() :
            QuickTest("time_stamp_test")
        {
        }

        virtual void run() const
        {
            TimeStamp infinity, a, b, c, d;

            TEST_CHECK(! (a < infinity));
            TEST_CHECK(! (b < infinity));
            TEST_CHECK(! (c < infinity));
            TEST_CHECK(! (d < infinity));

            a.take();
            sleep(1);

            b.take();
            sleep(1);

            c.take();
            sleep(1);

            d.take();

            TEST_CHECK(a < b);
            TEST_CHECK(b < c);
            TEST_CHECK(c < d);
            TEST_CHECK(d < infinity);
        }
} time_stamp_test;
