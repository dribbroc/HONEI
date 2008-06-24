/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <unittest/unittest.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/backends/memory/memory_backend.hh>

using namespace honei;
using namespace tests;

template <typename Tag_>
class MemoryBackendQuickTest :
    public QuickTaggedTest<Tag_>
{
    public:
        MemoryBackendQuickTest() :
            QuickTaggedTest<Tag_>("memory_backend_test")
        {
        }

        virtual void run() const
        {
            int data_array [10];
            for (int i(0) ; i < 10 ; ++i)
            {
                data_array[i] = i;
            }
            void * data = data_array;

            void * device(MemoryBackend<Tag_>::instance()->upload(0, data, 5 * sizeof(int)));
            if (device != data)
            {
                data_array[3] = -50;
            }
            MemoryBackend<Tag_>::instance()->download(0, data, 5 * sizeof(int));
            MemoryBackend<Tag_>::instance()->free(0);
            TEST_CHECK_EQUAL(data_array[3], 3);
        }
};
MemoryBackendQuickTest<tags::CPU> memory_backend_quick_test;
#ifdef HONEI_CUDA
MemoryBackendQuickTest<tags::GPU::CUDA> cuda_memory_backend_quick_test;
#endif

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
            char  data_array1 [10];
            void * data1 = data_array1;
            unsigned long mem1((unsigned long)data1);
            char  data_array2 [10];
            void * data2 = data_array2;
            unsigned long mem2((unsigned long)data2);

            MemoryArbiter::instance()->add_memblock(mem2);
            MemoryArbiter::instance()->add_memblock(mem1);
            MemoryArbiter::instance()->read(tags::CPU::memory_value, mem1, data1, 2);
            MemoryArbiter::instance()->read(tags::CPU::memory_value, mem2, data2, 1);
            MemoryArbiter::instance()->release_read(mem2);
            MemoryArbiter::instance()->release_read(mem1);
            MemoryArbiter::instance()->write(tags::CPU::memory_value, mem1, data1, 2);
            MemoryArbiter::instance()->release_write(mem1);

            TEST_CHECK_THROWS(MemoryArbiter::instance()->release_read(25), InternalError);
            TEST_CHECK_THROWS(MemoryArbiter::instance()->release_read(mem1), InternalError);

            MemoryArbiter::instance()->remove_memblock(mem2);
            MemoryArbiter::instance()->remove_memblock(mem1);
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
            char  data_array1 [10];
            void * data1 = data_array1;
            unsigned long mem1((unsigned long)data1);
            char  data_array2 [10];
            void * data2 = data_array2;
            unsigned long mem2((unsigned long)data2);
            MemoryArbiter::instance()->add_memblock(mem1);
            MemoryArbiter::instance()->read(tags::CPU::memory_value, mem1, data1, 1);
            MemoryArbiter::instance()->add_memblock(mem2);
            MemoryArbiter::instance()->read(tags::GPU::CUDA::memory_value, mem2, data2, 1);
            MemoryArbiter::instance()->read(tags::CPU::memory_value, mem2, data2, 1);
            MemoryArbiter::instance()->release_read(mem1);
            MemoryArbiter::instance()->release_read(mem2);
            MemoryArbiter::instance()->release_read(mem2);
            MemoryArbiter::instance()->write(tags::CPU::memory_value, mem1, data1, 1);
            MemoryArbiter::instance()->write(tags::CPU::memory_value, mem2, data2, 1);
            MemoryArbiter::instance()->release_write(mem1);
            MemoryArbiter::instance()->release_write(mem2);

            TEST_CHECK_THROWS(MemoryArbiter::instance()->release_read(25), InternalError);
            TEST_CHECK_THROWS(MemoryArbiter::instance()->release_read(mem1), InternalError);

            MemoryArbiter::instance()->remove_memblock(mem2);
            MemoryArbiter::instance()->remove_memblock(mem1);
        }
} cuda_memory_arbiter_quick_test;
#endif
