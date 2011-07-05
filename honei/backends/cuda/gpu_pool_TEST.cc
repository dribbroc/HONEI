/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@tu-dortmund.de>
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

#ifndef DEBUG
#define DEBUG 1
#endif

#include <honei/backends/cuda/gpu_pool.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/backends/cuda/multi_gpu.hh>
#include <honei/backends/cuda/transfer.hh>
#include <honei/util/memory_backend.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/util/lock.hh>
#include <honei/util/stringify.hh>
#include <honei/util/unittest.hh>

#include <iostream>

using namespace honei::cuda;
using namespace tests;

namespace
{
    class InfoTask
    {
        private:
            int thread;
        public:
            InfoTask(int thread) :
                thread(thread)
            {
            }

            void operator() ()
            {
                std::cout<<"GPU Thread: " << thread << std::endl;
                int device = cuda_get_device();
                std::cout<<"Device: " << device << std::endl;
                cuda_print_device_name(device);
            }
    };

    class TransferTask
    {
        private:
            int * data;
        public:
            TransferTask(int * d) :
                data(d)
            {
            }

            void operator() ()
            {
                MemoryBackend<tags::GPU::CUDA>::instance()->upload((void *) data, (void *)data, 5 * sizeof(int));
                for (int i(0) ; i < 5 ; ++i)
                {
                    data[i] = -4712;
                }
                MemoryBackend<tags::GPU::CUDA>::instance()->download((void *) data, (void *)data, 5 * sizeof(int));
                MemoryBackend<tags::GPU::CUDA>::instance()->free((void *) data);
            }
    };

    class LockTask
    {
        private:
            void * data;
        public:
            LockTask(void * d) :
                data(d)
            {
            }

            void operator() ()
            {
                //void * device(MemoryArbiter::instance()->lock(lm_read_and_write, tags::GPU::CUDA::memory_value, data, data, 1));
                MemoryArbiter::instance()->lock(lm_read_and_write, tags::GPU::CUDA::memory_value, data, data, 1);
                ((char*)data)[0]='b';
                //cuda_fill_zero(device, sizeof(char));
                MemoryArbiter::instance()->unlock(lm_read_and_write, data);
            }
    };
}

class GPUPoolQuickTest :
    public QuickTest
{
    public:
        GPUPoolQuickTest() :
            QuickTest("gpu_pool_quick_test")
        {
        }

        virtual void run() const
        {
            TicketVector tickets;

            int data [15];
            TransferTask trans(data);
            TransferTask trans2(data+5);
            TransferTask trans3(data+10);
            for (int i(0) ; i < 15 ; ++i)
            {
                data[i] = i;
            }
            tickets.push_back(GPUPool::instance()->enqueue(trans,0));
            tickets.wait();
            tickets.push_back(GPUPool::instance()->enqueue(trans2,1));
            tickets.wait();
            tickets.push_back(GPUPool::instance()->enqueue(trans3,2));
            tickets.wait();
            for (int i(0) ; i < 15 ; ++i)
            {
                TEST_CHECK_EQUAL(data[i], i);
            }
            TEST_CHECK(GPUPool::instance()->idle());
            std::cout<<"Used GPU's: "<<GPUPool::instance()->get_num_gpus()<<std::endl;

            for (unsigned long device(0) ; device < GPUPool::instance()->get_num_gpus() ; ++device)
            {
                InfoTask info(device);
                GPUPool::instance()->enqueue(info, device).wait();
            }

            GPUPool::instance()->single_start(0);
            TEST_CHECK_EQUAL(GPUPool::instance()->get_num_gpus(), 1ul);
        }
} gpu_pool_quick_test;

class GPUPoolArbiterQuickTest :
    public QuickTest
{
    public:
        GPUPoolArbiterQuickTest() :
            QuickTest("gpu_pool_arbiter_quick_test")
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
            char  data_array3 [10];
            void * data3 = data_array3;
            void * mem3(data3);
            data_array1[0]='a';
            data_array2[0]='a';
            data_array3[0]='a';
            MemoryArbiter::instance()->register_address(mem1);
            MemoryArbiter::instance()->register_address(mem2);
            MemoryArbiter::instance()->register_address(mem3);
            TicketVector tickets;
            LockTask lt1(data1);
            LockTask lt2(data2);
            LockTask lt3(data3);

            tickets.push_back(GPUPool::instance()->enqueue(lt1,0));
            tickets.push_back(GPUPool::instance()->enqueue(lt2,1));
            tickets.push_back(GPUPool::instance()->enqueue(lt3,2));
            tickets.wait();

            MemoryArbiter::instance()->lock(lm_read_only, tags::CPU::memory_value, mem1, data1, 1);
            MemoryArbiter::instance()->lock(lm_read_only, tags::CPU::memory_value, mem2, data2, 1);
            MemoryArbiter::instance()->lock(lm_read_only, tags::CPU::memory_value, mem3, data3, 1);
            TEST_CHECK_EQUAL(data_array1[0], 'a');
            TEST_CHECK_EQUAL(data_array2[0], 'a');
            TEST_CHECK_EQUAL(data_array3[0], 'a');
            MemoryArbiter::instance()->unlock(lm_read_only, mem1);
            MemoryArbiter::instance()->unlock(lm_read_only, mem2);
            MemoryArbiter::instance()->unlock(lm_read_only, mem3);

            MemoryArbiter::instance()->remove_address(mem2);
            MemoryArbiter::instance()->remove_address(mem1);
            MemoryArbiter::instance()->remove_address(mem3);
        }
} gpu_pool_arbiter_quick_test;
