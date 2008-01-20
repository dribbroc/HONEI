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

#ifndef LIBUTIL_GUARD_TIME_STAMP_HH
#define LIBUTIL_GUARD_TIME_STAMP_HH 1

#ifdef linux
#  include <sys/time.h>
#else
#  error "TimeStamp is not yet supported on your platform!"
#endif

#include <limits>

namespace honei
{
    /**
     * \ingroup grplibutil
     */
    class TimeStamp
    {
        private:
            /// Our time-stamp.
            timeval _time;

        public:
            /// Constructor.
            TimeStamp()
            {
                _time.tv_sec = std::numeric_limits<typeof(_time.tv_sec)>::max();
                _time.tv_usec = std::numeric_limits<typeof(_time.tv_usec)>::max();
            }

            /// Take a new time stamp.
            void take()
            {
                gettimeofday(&_time, 0);
            }

            /// Returns our seconds.
            unsigned long sec()
            {
                return _time.tv_sec;
            }

            /// Returns our useconds.
            long usec()
            {
                return _time.tv_usec;
            }

            /**
             * Our comparison operator.
             *
             * Return true if our time-stamp has been taken earlier than the
             * other.
             *
             * \param other Another time-stamp.
             */
            bool operator< (const TimeStamp & other)
            {
                if (_time.tv_sec < other._time.tv_sec)
                {
                    return true;
                }
                else if (_time.tv_sec == other._time.tv_sec)
                {
                    return _time.tv_usec < other._time.tv_usec;
                }
                else
                {
                    return false;
                }
            }
    };
}

#endif
