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

#include <honei/util/pool_task.hh>
#include <honei/util/thread_pool.hh>
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

class PoolTaskTest :
    public QuickTest
{
    public:
        PoolTaskTest() :
            QuickTest("pool_task_test")
        {
        }

        virtual void run() const
        {
            unsigned v(34);
            TestTask t(v);
//            ThreadPool * pool = ThreadPool.getInstance(0);
            WorkerTask wt(t);
            wt();
            TEST_CHECK_EQUAL(v, 35);
            PoolTask pt(wt);
            pt();
            TEST_CHECK_EQUAL(v, 36);
        }
} pool_task_test;
