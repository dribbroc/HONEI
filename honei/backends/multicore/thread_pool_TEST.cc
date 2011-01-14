/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008, 2009, 2010, 2011 Sven Mallach <mallach@honei.org>
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

#include <honei/backends/multicore/thread_pool.hh>
#include <honei/util/lock.hh>
#include <honei/util/stringify.hh>
#include <honei/util/unittest.hh>

using namespace honei::mc;
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

            void operator() ()
            {
                Lock l(*_mutex);
                ++_v;
            }
    };
}

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

            for (unsigned i(0) ; i < 500 ; ++i)
            {
                tickets.push_back(ThreadPool::instance()->enqueue(t));
            }
            for (unsigned i(500) ; i < 750 ; ++i)
            {
                tickets.push_back(ThreadPool::instance()->enqueue(t, DispatchPolicy::on_core(2 * sysconf(_SC_NPROCESSORS_CONF))));
            }

            for (unsigned i(750) ; i < 1000 ; ++i)
            {
                tickets.push_back(ThreadPool::instance()->enqueue(t, DispatchPolicy::on_core(i % (sysconf(_SC_NPROCESSORS_CONF) - 1))));
            }

            for (unsigned i(1000) ; i < 1250 ; ++i)
            {
                Ticket<tags::CPU::MultiCore> * tick(static_cast<Ticket<tags::CPU::MultiCore> *>(tickets[750 + (i % 200)]));
                tickets.push_back(ThreadPool::instance()->enqueue(t, DispatchPolicy::same_core_as(tick)));
            }

            Topology * top = Topology::instance();

            for (unsigned i(1250) ; i < 1500 ; ++i)
            {
                tickets.push_back(ThreadPool::instance()->enqueue(t, DispatchPolicy::on_node(i % top->num_nodes())));
            }

            for (unsigned i(1500) ; i < 1750 ; ++i)
            {
                Ticket<tags::CPU::MultiCore> * tick(static_cast<Ticket<tags::CPU::MultiCore> *>(tickets[1250 + (i % 200)]));
                tickets.push_back(ThreadPool::instance()->enqueue(t, DispatchPolicy::same_node_as(tick)));
            }

            tickets.wait();

            TEST_CHECK_EQUAL(v, 1784u);
        }
} thread_pool_test;

class ThreadPoolQuickTest :
    public QuickTest
{
    public:
        ThreadPoolQuickTest() :
            QuickTest("thread_pool_quick_test")
        {
        }

        virtual void run() const
        {
            unsigned v(34);
            TestTask t(v);

            TicketVector tickets;

            for (unsigned i(0) ; i < 125 ; ++i)
            {
                tickets.push_back(ThreadPool::instance()->enqueue(t));
            }

            Ticket<tags::CPU::MultiCore> * first_t(static_cast<Ticket<tags::CPU::MultiCore> *>(tickets[0]));

            for (unsigned i(125) ; i < 250 ; ++i)
            {
                tickets.push_back(ThreadPool::instance()->enqueue(t, DispatchPolicy::on_core(2 * sysconf(_SC_NPROCESSORS_CONF))));
            }

            for (unsigned i(250) ; i < 500 ; ++i)
            {
                tickets.push_back(ThreadPool::instance()->enqueue(t, DispatchPolicy::same_core_as(first_t)));
            }

            tickets.wait();

            TEST_CHECK_EQUAL(v, 534u);
        }
} thread_pool_quick_test;
