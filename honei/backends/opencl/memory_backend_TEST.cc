/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/util/unittest.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/util/memory_backend.hh>
#include <honei/backends/opencl/opencl_backend.hh>
#include <honei/backends/opencl/operations.hh>

using namespace honei;
using namespace tests;

template <typename Tag_>
class MemoryBackendQuickTest :
    public QuickTaggedTest<Tag_>
{
    public:
        MemoryBackendQuickTest() :
            QuickTaggedTest<Tag_>("opencl_memory_backend_test")
        {
        }

        virtual void run() const
        {
            float data_array [10];
            for (int i(0) ; i < 10 ; ++i)
            {
                data_array[i] = i;
            }
            void * data = data_array;

            MemoryBackend<Tag_>::instance()->upload(0, data, 5 * sizeof(float));

            if (Tag_::memory_value == honei::tags::tv_opencl_gpu)
                data_array[3] = -50;
            MemoryBackend<Tag_>::instance()->download(0, data, 5 * sizeof(float));
            MemoryBackend<Tag_>::instance()->free(0);
            TEST_CHECK_EQUAL(data_array[3], 3);
        }
};
MemoryBackendQuickTest<tags::OpenCL::CPU> ocl_cpu_memory_backend_quick_test;
#ifdef HONEI_OPENCL_GPU
MemoryBackendQuickTest<tags::OpenCL::GPU> ocl_gpu_memory_backend_quick_test;
#endif
#ifdef HONEI_OPENCL_ACC
MemoryBackendQuickTest<tags::OpenCL::Accelerator> ocl_accelerator_memory_backend_quick_test;
#endif

template <typename Tag_>
class OpenCLMemoryArbiterQuickTest :
    public QuickTaggedTest<Tag_>
{
    public:
        OpenCLMemoryArbiterQuickTest() :
            QuickTaggedTest<Tag_>("opencl_memory_arbiter_test")
        {
        }

        virtual void run() const
        {
            char  data_array1 [10];
            void * data1 = data_array1;
            void * mem1(data1);
            char  data_array2 [10];
            void * data2 = data_array2;
            void * mem2(data2);
            MemoryArbiter::instance()->register_address(mem1);
            MemoryArbiter::instance()->lock(lm_read_only, tags::CPU::memory_value, mem1, data1, 1);
            MemoryArbiter::instance()->register_address(mem2);
            MemoryArbiter::instance()->lock(lm_read_only, Tag_::memory_value, mem2, data2, 1);
            MemoryArbiter::instance()->lock(lm_read_only, tags::CPU::memory_value, mem2, data2, 1);
            MemoryArbiter::instance()->unlock(lm_read_only, mem1);
            MemoryArbiter::instance()->unlock(lm_read_only, mem2);
            MemoryArbiter::instance()->unlock(lm_read_only, mem2);
            MemoryArbiter::instance()->lock(lm_read_and_write, Tag_::memory_value, mem2, data2, 1);
            MemoryArbiter::instance()->unlock(lm_read_and_write, mem2);
            MemoryArbiter::instance()->lock(lm_read_only, tags::CPU::memory_value, mem2, data2, 1);
            MemoryArbiter::instance()->unlock(lm_read_only, mem2);
            MemoryArbiter::instance()->lock(lm_read_and_write, tags::CPU::memory_value, mem1, data1, 1);
            MemoryArbiter::instance()->lock(lm_read_and_write, tags::CPU::memory_value, mem2, data2, 1);
            MemoryArbiter::instance()->unlock(lm_read_and_write, mem1);
            MemoryArbiter::instance()->unlock(lm_read_and_write, mem2);

            TEST_CHECK_THROWS(MemoryArbiter::instance()->unlock(lm_read_only, (void *)25), InternalError);
            TEST_CHECK_THROWS(MemoryArbiter::instance()->unlock(lm_read_only, mem1), InternalError);

            MemoryArbiter::instance()->remove_address(mem2);
            MemoryArbiter::instance()->remove_address(mem1);
        }
};
OpenCLMemoryArbiterQuickTest<tags::OpenCL::CPU> ocl_cpu_memory_arbiter_quick_test;
#ifdef HONEI_OPENCL_GPU
OpenCLMemoryArbiterQuickTest<tags::OpenCL::GPU> ocl_gpu_memory_arbiter_quick_test;
#endif
#ifdef HONEI_OPENCL_ACC
OpenCLMemoryArbiterQuickTest<tags::OpenCL::Accelerator> ocl_accelerator_memory_arbiter_quick_test;
#endif
