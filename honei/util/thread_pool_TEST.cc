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

#include <honei/util/lock.hh>
#include <honei/util/pool_task.hh>
#include <honei/util/pool_thread.hh>
#include <honei/util/stringify.hh>
#include <honei/util/thread_pool.hh>
#include <unittest/unittest.hh>

#include <iostream>

using namespace honei;
using namespace tests;

namespace
{
    class TestTask
    {
        private:
            unsigned & _v;
            Mutex * const _mutex;

        public:
            TestTask(unsigned & v) :
                _v(v),
                _mutex(new Mutex)
            {
            }

            virtual void operator() ()
            {
                Lock l(*_mutex);
                ++_v;
            }
    };
}

class ThreadPoolTest :
    public QuickTest
{
    public:
        ThreadPoolTest() :
            QuickTest("thread_pool_test")
        {
        }

        virtual void run() const
        {
            unsigned v(34);
            TestTask t(v);
            WorkerTask wt(t);
            PoolTask * pt[500];
            for (unsigned i(0) ; i < 500 ; ++i)
            {
                pt[i] = ThreadPool::instance()->dispatch(wt);
            }
            for (unsigned i(0) ; i < 500 ; ++i)
            {
                pt[i]->wait_on();
            }
            TEST_CHECK_EQUAL(v, 534);
        }
} thread_pool_test;
