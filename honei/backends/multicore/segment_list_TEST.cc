/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Sven Mallach <mallach@honei.org>
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

#include <honei/backends/multicore/segment_list.hh>
#include <honei/backends/multicore/thread_pool.hh>

#include <honei/util/lock.hh>

#include <honei/util/stringify.hh>
#include <unittest/unittest.hh>

#include <list>

using namespace honei::mc;
using namespace tests;

namespace
{
    class TestTask
    {
        private:
            Mutex * const _mutex;

        public:
            int exec_ctr;

            TestTask() :
                exec_ctr(0),
                _mutex(new Mutex)
            {
            }

            void operator() ()
            {
                Lock l(*_mutex);
                ++exec_ctr;
            }
    };

    class EnqueuerTask
    {
        private:

            TicketVector * tv;
            std::list<TestTask *> * l;

        public:

            EnqueuerTask() :
                tv(new TicketVector),
                l(new std::list<TestTask *>)

            {
            }

            EnqueuerTask(const EnqueuerTask & other) :
                l(other.l),
                tv(other.tv)
            {
            }

            bool check()
            {
                tv->wait();

                while (! l->empty())
                {
                    TestTask * t = l->front();

                    if (t->exec_ctr != 1)
                        return false;

                    l->pop_front();

                    delete t;
                }

                delete tv;
                delete l;

                return true;
            }

            void operator() ()
            {
                for (int i(0) ; i < 1 ; ++i)
                {
                    TestTask * t = new TestTask;
                    l->push_back(t);

                    (*t)();
                    tv->push_back(ThreadPool::instance()->enqueue(*t));
                }
            }
    };
}

class SegmentListTest :
    public BaseTest
{
    public:
        SegmentListTest() :
            BaseTest("segment_list_test")
        {
        }

        virtual void run() const
        {
            std::list<EnqueuerTask *> l;

            TicketVector tickets;

            for (unsigned i(0) ; i < 50 ; ++i)
            {
                EnqueuerTask * et = new EnqueuerTask;
                l.push_back(et);

                tickets.push_back(ThreadPool::instance()->enqueue((*et)));
            }

            Topology * top = Topology::instance();

            for (unsigned i(50) ; i < 100 ; ++i)
            {
                EnqueuerTask * et = new EnqueuerTask;
                l.push_back(et);

                tickets.push_back(ThreadPool::instance()->enqueue(*et, DispatchPolicy::on_node(i % top->num_nodes())));
            }

            for (unsigned i(100) ; i < 150 ; ++i)
            {
                EnqueuerTask * et = new EnqueuerTask;
                l.push_back(et);

                tickets.push_back(ThreadPool::instance()->enqueue(*et, DispatchPolicy::on_core(i % top->num_lpus())));
            }

            tickets.wait();

            while (! l.empty())
            {
                EnqueuerTask * et = l.front();
                if (! et->check())
                {
                    TEST_CHECK(false);
                }

                l.pop_front();
                delete et;
            }

            TEST_CHECK(true);
        }
} segment_list_test;
