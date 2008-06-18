/* vim: set sw=4 sts=4 et foldmethod=syntax : */


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
            MemoryArbiter::instance()->read<tags::CPU>(0);
            MemoryArbiter::instance()->read<tags::CPU>(1);
            MemoryArbiter::instance()->release_read<tags::CPU>(1);
            MemoryArbiter::instance()->release_read<tags::CPU>(0);
            MemoryArbiter::instance()->write<tags::CPU>(0);
            MemoryArbiter::instance()->release_write<tags::CPU>(0);

            TEST_CHECK_THROWS(MemoryArbiter::instance()->release_read<tags::CPU>(25), InternalError);
            TEST_CHECK_THROWS(MemoryArbiter::instance()->release_read<tags::CPU>(0), InternalError);
        }
} memory_arbiter_quick_test;


#ifdef HONEI_CUDA
class CUDAMemoryArbiterQuickTest :
    public QuickTest
{
    public:
        CUDAMemoryArbiterQuickTest() :
            QuickTest("cuda_memory_arbiter_test")
        {
        }

        virtual void run() const
        {
            MemoryArbiter::instance()->read<tags::CPU>(0);
            MemoryArbiter::instance()->read<tags::GPU::CUDA>(0);
            MemoryArbiter::instance()->read<tags::CPU>(1);
            MemoryArbiter::instance()->release_read<tags::CPU>(1);
            MemoryArbiter::instance()->release_read<tags::CPU>(0);
            MemoryArbiter::instance()->release_read<tags::GPU::CUDA>(0);
            MemoryArbiter::instance()->write<tags::CPU>(0);
            MemoryArbiter::instance()->release_write<tags::CPU>(0);

            TEST_CHECK_THROWS(MemoryArbiter::instance()->release_read<tags::GPU::CUDA>(25), InternalError);
            TEST_CHECK_THROWS(MemoryArbiter::instance()->release_read<tags::GPU::CUDA>(0), InternalError);
        }
} cuda_memory_arbiter_quick_test;
#endif
