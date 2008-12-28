/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Sven Mallach <sven.mallach@cs.tu-dortmund.de>
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
#include <unittest/unittest.hh>

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

            virtual void operator() ()
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

            TicketList tickets;

            std::tr1::shared_ptr<Ticket<tags::CPU::MultiCore> > first_t(ThreadPool::instance(4 * sysconf(_SC_NPROCESSORS_CONF), false)->enqueue(t));

            tickets.push_back(first_t);

            for (unsigned i(1) ; i < 1250 ; ++i)
            {
                tickets.push_back(ThreadPool::instance()->enqueue(t));
            }

            for (unsigned i(1250) ; i < 2500 ; ++i)
            {
                tickets.push_back(ThreadPool::instance()->enqueue(t, DispatchPolicy::on_core(2 * sysconf(_SC_NPROCESSORS_CONF))));
            }

            for (unsigned i(2500) ; i < 5000 ; ++i)
            {
                tickets.push_back(ThreadPool::instance()->enqueue(t, DispatchPolicy::same_core_as(first_t)));
            }

            tickets.wait();

            ThreadPool::instance()->destroy();

            TEST_CHECK_EQUAL(v, 5034);

            first_t = ThreadPool::instance(4 * sysconf(_SC_NPROCESSORS_CONF), true)->enqueue(u);
            tickets.push_back(first_t);

            for (unsigned i(1) ; i < 1250 ; ++i)
            {
                tickets.push_back(ThreadPool::instance()->enqueue(u));
            }

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

            ThreadPool::instance()->destroy();

            TEST_CHECK_EQUAL(w, 5034);
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

            TicketList tickets;

            std::tr1::shared_ptr<Ticket<tags::CPU::MultiCore> > first_t(ThreadPool::instance(4 * sysconf(_SC_NPROCESSORS_CONF), true)->enqueue(t));

            tickets.push_back(first_t);

            for (unsigned i(1) ; i < 125 ; ++i)
            {
                tickets.push_back(ThreadPool::instance()->enqueue(t));
            }

            for (unsigned i(125) ; i < 250 ; ++i)
            {
                tickets.push_back(ThreadPool::instance()->enqueue(t, DispatchPolicy::on_core(2 * sysconf(_SC_NPROCESSORS_CONF))));
            }

            for (unsigned i(250) ; i < 500 ; ++i)
            {
                tickets.push_back(ThreadPool::instance()->enqueue(t, DispatchPolicy::same_core_as(first_t)));
            }

            tickets.wait();

            ThreadPool::instance()->destroy();

            TEST_CHECK_EQUAL(v, 534);
        }
} thread_pool_quick_test;
