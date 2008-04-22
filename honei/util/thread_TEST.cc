/* vim: set number sw=4 sts=4 et foldmethod=syntax : */

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

// This test needs DEBUG defined.
#ifndef DEBUG
#define DEBUG 1
#endif

#include <honei/util/lock.hh>
#include <honei/util/mutex.hh>
#include <honei/util/thread.hh>
#include <unittest/unittest.hh>

#include <string>
#include <vector>

using namespace honei;
using namespace tests;

namespace
{
    class Counter
    {
        private:
            std::tr1::shared_ptr<unsigned> _value;

            std::tr1::shared_ptr<Mutex> _mutex;

        public:
            Counter() :
                _value(new unsigned(0)),
                _mutex(new Mutex)
            {
            }

            void operator() ()
            {
                Lock l(*_mutex);

                ++(*_value);
            }

            unsigned value()
            {
                return *_value;
            }
    };
}

class ThreadTest :
    public QuickTest
{
    public:
        ThreadTest() :
            QuickTest("thread_test")
        {
        }

        virtual void run() const
        {
            Counter c;
            Thread::Function f(c);
            std::vector<Thread *> threads(16);

            for (std::vector<Thread *>::iterator i(threads.begin()), i_end(threads.end()) ;
                    i != i_end ; ++i)
            {
                *i = new Thread(f);
            }

            bool done;
            do
            {
                done = true;
                for (std::vector<Thread *>::iterator i(threads.begin()), i_end(threads.end()) ;
                        i != i_end ; ++i)
                {
                    done &= (*i)->completed();
                }
            }
            while (! done);

            TEST_CHECK_EQUAL(c.value(), threads.size());
        }
} thread_test;
