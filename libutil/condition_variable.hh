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

#ifndef LIBUTIL_GUARD_CONDITION_VARIABLE_HH
#define LIBUTIL_GUARD_CONDITION_VARIABLE_HH 1

#include <libutil/mutex.hh>

#include <pthread.h>

namespace honei
{
    class ConditionVariable
    {
        private:
            /// Our pthread condition variable.
            pthread_cond_t * const _cond;

            /// Unwanted copy-constructor: Do not implement. See EffCpp, Item 27.
            ConditionVariable(const ConditionVariable &);

            /// Unwanted assignment operator: Do not implement. See EffCpp, Item 27.
            ConditionVariable & operator= (const ConditionVariable &);

        public:
            /// Constructor.
            ConditionVariable();

            /// Destructor.
            ~ConditionVariable();

            /// Broadcast a wake-up to all waiting threads.
            void broadcast();

            /// Sign a wake-up to the waiting thread.
            void signal();

            /// Acquire a lock for signaling.
            void acquire_then_signal(Mutex &);

            /// Wait until condition is met.
            void wait(Mutex &);
    };
}

#endif
