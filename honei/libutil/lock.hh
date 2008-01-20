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

#ifndef LIBUTIL_GUARD_LOCK_HH
#define LIBUTIL_GUARD_LOCK_HH 1

#include <honei/libutil/mutex.hh>

namespace honei
{
    class Lock
    {
        private:
            /// Our mutex.
            Mutex * const _mutex;

            /// Unwanted copy-constructor: Do not implement. See EffCpp, Item 27.
            Lock(const Lock &);

            /// Unwanted assignment operator: Do not implement. See EffCpp, Item 27.
            Lock & operator= (const Lock &);

        public:
            /// (Explicit) constructor.
            explicit Lock(Mutex &);

            /// Destructor.
            ~Lock();
    };

    class TryLock
    {
        private:
            /// Our mutex.
            Mutex * _mutex;

            /// Unwanted copy-constructor: Do not implement. See EffCpp, Item 27.
            TryLock(const TryLock &);

            /// Unwanted assignment operator: Do not implement. See EffCpp, Item 27.
            TryLock & operator= (const TryLock &);

        public:
            /// (Explicit) constructor.
            explicit TryLock(Mutex & e);

            /// Destructor.
            ~TryLock();

            /// Return true if the lock worked.
            bool operator() () const
            {
                return 0 != _mutex;
            }
    };
}

#endif
