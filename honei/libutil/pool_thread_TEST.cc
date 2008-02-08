/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Volker Jung <volker.jung@uni-dortmund.de>
 *
 * This file is part of the Utility C++ library. LibUtil is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibUtil is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

// This test case needs debug support.
#ifndef DEBUG
#define DEBUG 1
#endif

#include <honei/libutil/pool_task.hh>
#include <honei/libutil/pool_thread.hh>
#include <honei/libutil/thread_pool.hh>
#include <unittest/unittest.hh>

using namespace honei;
using namespace tests;

namespace
{
    class TestTask
    {
        private:
            unsigned & _v;

        public:
            TestTask(unsigned & v) :
                _v(v)
            {
            }

            virtual void operator() ()
            {
                ++_v;
            }
    };
}

class PoolThreadTest :
    public QuickTest
{
    public:
        PoolThreadTest() :
            QuickTest("pool_thread_test")
        {
        }

        virtual void run() const
        {
            unsigned v(34);
            TestTask t(v);
            WorkerTask wt(t);
            PoolThread * poolthread(new PoolThread(ThreadPool::instance()));
            wt();
            TEST_CHECK_EQUAL(v, 35);
            PoolTask * pt(new PoolTask(wt));
            poolthread->run(pt);
            pt->wait_on();
            TEST_CHECK_EQUAL(v, 36);
       }
} pool_thread_test;
