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
#include <honei/backends/cuda/multi_gpu.hh>
#include <honei/util/lock.hh>
#include <honei/util/stringify.hh>
#include <unittest/unittest.hh>
#include <iostream>

using namespace honei::cuda;
using namespace tests;

namespace
{
    class TestTask
    {
        private:
            Mutex * const _mutex;

        public:
            TestTask() :
                _mutex(new Mutex)
            {
            }

            void operator() ()
            {
                Lock l(*_mutex);
                std::cout<<cuda_get_device()<<std::endl;
            }
    };
}
/*
class ThreadPoolTest :
    public BaseTest
{
    public:
        ThreadPoolTest() :
            BaseTest("thread_pool_test")
        {
        }

        virtual void run() const
        {
            unsigned v(34), w(34);
            TestTask t(v);
            TestTask u(w);

            TicketVector tickets;

            for (unsigned i(0) ; i < 1250 ; ++i)
            {
                tickets.push_back(ThreadPool::instance()->enqueue(t));
            }

            Ticket<tags::CPU::MultiCore> * first_t(static_cast<Ticket<tags::CPU::MultiCore> *>(tickets[0]));

            for (unsigned i(1250) ; i < 2500 ; ++i)
            {
                tickets.push_back(ThreadPool::instance()->enqueue(t, DispatchPolicy::on_core(2 * sysconf(_SC_NPROCESSORS_CONF))));
            }

            for (unsigned i(2500) ; i < 5000 ; ++i)
            {
                tickets.push_back(ThreadPool::instance()->enqueue(t, DispatchPolicy::same_core_as(first_t)));
            }

            tickets.wait();

            TEST_CHECK_EQUAL(v, 5034u);

            for (unsigned i(0) ; i < 1250 ; ++i)
            {
                tickets.push_back(ThreadPool::instance()->enqueue(u));
            }

            first_t = static_cast<Ticket<tags::CPU::MultiCore> *>(tickets[251]);

            ThreadPool::instance()->add_threads(3);

            for (unsigned i(1250) ; i < 2500 ; ++i)
            {
                tickets.push_back(ThreadPool::instance()->enqueue(u, DispatchPolicy::on_core(2 * sysconf(_SC_NPROCESSORS_CONF))));
            }

            ThreadPool::instance()->delete_threads(3);

            for (unsigned i(2500) ; i < 5000 ; ++i)
            {
                tickets.push_back(ThreadPool::instance()->enqueue(u, DispatchPolicy::same_core_as(first_t)));
            }

            tickets.wait();

            TEST_CHECK_EQUAL(w, 5034u);
        }
} thread_pool_test;
*/

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
            TestTask t;

            TicketVector tickets;

            tickets.push_back(GPUPool::instance()->enqueue(t,1));
            tickets.wait();
            tickets.push_back(GPUPool::instance()->enqueue(t,0));
            tickets.wait();
            tickets.push_back(GPUPool::instance()->enqueue(t,1));
            tickets.push_back(GPUPool::instance()->enqueue(t,0));
            tickets.wait();


            //TEST_CHECK_EQUAL(v, 534u);
        }
} gpu_pool_quick_test;
