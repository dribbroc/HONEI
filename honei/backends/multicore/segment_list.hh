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

#ifndef MULTICORE_GUARD_ATOMICLIST_HH
#define MULTICORE_GUARD_ATOMICLIST_HH 1

#include <honei/backends/multicore/thread_task.hh>
#include <honei/util/mutex.hh>
#include <tr1/functional>

#include <iostream>

namespace honei
{
    struct SegmentEntry
    {
        mc::ThreadTask * const data;
        SegmentEntry * next;

        SegmentEntry(mc::ThreadTask * const t) :
            data(t),
            next(NULL)
        {
        }
    };

    struct Segment
    {
        SegmentEntry * head;
        unsigned count;
        Mutex * mutex;

        Segment() :
            head(NULL),
            count(0),
            mutex(new Mutex)
        {
        }

        ~Segment()
        {
            delete mutex;
        }

        void lock()
        {
            pthread_mutex_lock(mutex->mutex());
        }

        void unlock()
        {
            pthread_mutex_unlock(mutex->mutex());
        }

        mc::ThreadTask * find_and_extract(std::tr1::function<bool (mc::ThreadTask *)> & comp)
        {
            SegmentEntry * prev(NULL);
            SegmentEntry * it(head);

            mc::ThreadTask * res(NULL);

            while (it != NULL)
            {
                if (comp(it->data))
                {
                    if (prev != NULL)
                        prev->next = it->next;
                    else
                        head = it->next;

                    --count;
                    res = it->data;
                    delete it;
                    break;
                }

                prev = it;
                it = it->next;
            }

            return res;
        }

        void push(mc::ThreadTask * const t)
        {
            if (count == 0)
            {
//                std::cout << "push. count = 0, head = " << head << std::endl;
                head = new SegmentEntry(t);
                ++count;
                return;
            }

            SegmentEntry * h(head);

            for (unsigned i(1) ; i < count; ++i)
                h = h->next;

//            std::cout << "push. count = " << count << ",  using h = " << h << std::endl;

            h->next = new SegmentEntry(t);
            ++count;
        }
    };

    class SegmentList
    {
        private:
            /// \name Private members
            /// \{

                const unsigned _num_segs;
                volatile unsigned _size;

                // A fixed size array of segments
                Segment ** _segs;

                Mutex * _global_mutex;

            public:

                SegmentList(unsigned procs);

                ~SegmentList();

                void push_back(mc::ThreadTask * elem);

                void consistency_check();

                mc::ThreadTask * extract(std::tr1::function<bool (mc::ThreadTask *)> & comp);
    };
}
#endif
