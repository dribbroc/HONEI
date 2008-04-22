/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#ifndef LIBUTIL_GUARD_SYNC_POINT_HH
#define LIBUTIL_GUARD_SYNC_POINT_HH 1

#include <honei/util/condition_variable.hh>
#include <honei/util/instantiation_policy.hh>
#include <honei/util/lock.hh>
#include <honei/util/mutex.hh>

namespace honei
{
    class SyncPoint :
        public InstantiationPolicy<SyncPoint, NonCopyable>
    {
        private:
            /// Condition variable for synchronisation state.
            ConditionVariable * const _synchronised;

            /// Our mutex.
            Mutex * const _mutex;

            /// Our counter.
            unsigned _counter;

            /// Our sync-count.
            unsigned _sync_count;

        public:
            /// Constructor.
            SyncPoint(unsigned sync_count) :
                _synchronised(new ConditionVariable),
                _mutex(new Mutex),
                _counter(0),
                _sync_count(sync_count)
            {
            }

            /// Destructor.
            ~SyncPoint()
            {
                delete _mutex;
                delete _synchronised;
            }

            /// Increment our counter and return it.
            unsigned inline increment()
            {
                Lock l(*_mutex);

                return ++_counter;
            }

            /// Return our counter.
            unsigned inline counter()
            {
                Lock l(*_mutex);

                return _counter;
            }

            /// Signal the reach of this sync point in one of the critical sections.
            void signal_and_wait(Mutex * sync_mutex)
            {
                if (increment() == _sync_count)
                {
                    _synchronised->broadcast();
                }
                else
                {
                    do
                    {
                        _synchronised->wait(*sync_mutex);
                    } while (counter() != _sync_count);
                }
            }
    };
}

#endif
