/* vim: set sw=4 sts=4 et foldmethod=syntax : */


#include <honei/util/time_stamp.hh>
#include <unittest/unittest.hh>
#include <honei/backends/memory/memory_arbiter.hh>

using namespace honei;
using namespace tests;

class MemoryArbiterQuickTest :
    public QuickTest
{
    public:
        MemoryArbiterQuickTest() :
            QuickTest("memory_arbiter_test")
        {
        }

        virtual void run() const
        {
            MemoryArbiter::instance()->read<tags::CPU>(0, 1);
            MemoryArbiter::instance()->release_read<tags::CPU>(0, 1);
            MemoryArbiter::instance()->write<tags::CPU>(0, 1);
            MemoryArbiter::instance()->release_write<tags::CPU>(0, 1);

            TEST_CHECK_THROWS(MemoryArbiter::instance()->release_read<tags::CPU>(25, 50), InternalError);
            TEST_CHECK_THROWS(MemoryArbiter::instance()->release_read<tags::CPU>(0, 1), InternalError);
        }
} memory_arbiter_quick_test;
