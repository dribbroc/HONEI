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

#ifndef LIBUTIL_GUARD_LOG_HH
#define LIBUTIL_GUARD_LOG_HH 1

#include <libutil/mutex.hh>

#include <string>
#include <list>

namespace honei
{
    enum LogLevel
    {
        ll_full, ///< Output everything, even structure contents.
        ll_stubs, ///< Output important data alongside the usual message.
        ll_minimal ///< Output only minimal data.
    };

    class Log
    {
        private:
            class LogOutput;

            /// Our list of outputs.
            std::list<LogOutput> _outputs;

            /// Our mutex.
            Mutex * const _mutex;

            /// Constructor.
            Log();

        public:
            /// Return the singleton instance of Log.
            static Log * instance();

            /**
             * Log a message.
             *
             * \param level Log-level of the message.
             * \param msg Message to be logged.
             **/
            void message(const LogLevel level, const std::string & msg);
    };
}

#endif
