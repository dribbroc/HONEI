/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * Base in part upon 'log.hh' from Paludis, which is:
 *     Copyright (c) 2006, 2007, 2008 Ciaran McCreesh
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

#include <honei/libutil/instantiation_policy.hh>

#include <string>

namespace honei
{
    struct LogQueue;

    /**
     * LogLevel governs the severity of a given LogMessage.
     */
    enum LogLevel
    {
        ll_full, ///< Output everything, even structure contents.
        ll_stubs, ///< Output important data alongside the usual message.
        ll_minimal ///< Output only minimal data.
    };

    /**
     * LogMessage enqueues a log message with the LogMessageQueue.
     */
    class LogMessage :
        public InstantiationPolicy<LogMessage, NonCopyable>
    {
        public:
            friend struct LogQueue;

            /**
             * Constructor.
             *
             * \param level Log-level of the message.
             * \param messag Message to be logged.
             */
            LogMessage(const LogLevel level, const std::string & message);

            /// Destructor.
            ~LogMessage();
    };

    /// Log a message.
    static inline void log(LogLevel level, const std::string & message)
    {
        LogMessage(level, message);
    }

/**
 * \def LOGMESSAGE
 *
 * \brief Convenience definition that provides a way to declare uniquely-named
 * instances of class LogMessage.
 *
 * The created LogMessage will be automatically enqueued with the LogQueue and written.
 *
 * \param l Log level of the message.
 * \param m The message.
 *
 * \warning Will only be compiled in when debug support is enabled.
 *
 * \ingroup grpdebug
 */
#if defined (DEBUG)
// C preprocessor abomination following...
#define LOGMESSAGE(l, m) log(l, m)
#else
#define LOGMESSAGE(l, m)
#endif
}

#endif
