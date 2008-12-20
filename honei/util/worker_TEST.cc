/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#include <honei/util/worker.hh>
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

            virtual ~TestTask(){}

            virtual void operator() ()
            {
                ++_v;
            }
    };
}

class WorkerQueueTest :
    public QuickTest
{
    public:
        WorkerQueueTest() :
            QuickTest("worker_queue_test")
        {
        }

        virtual void run() const
        {
            unsigned v(34);
            TestTask t(v);
            WorkerThread thread;
            WorkerTask wt(t);

            for (unsigned i(0) ; i < 8 ; ++i)
            {
                thread.enqueue(wt);
            }

            while (! thread.idle())
                sleep(1);

            TEST_CHECK_EQUAL(v, 42ul);
        }
} worker_queue_test;
