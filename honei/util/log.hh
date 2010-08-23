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

#pragma once
#ifndef LIBUTIL_GUARD_LOG_HH
#define LIBUTIL_GUARD_LOG_HH 1

#include <honei/util/instantiation_policy.hh>

#include <string>

namespace honei
{
    struct LogQueue;

    /**
     * LogCategory governs the severity of a given LogMessage.
     */
    enum LogCategory
    {
        ll_minimal, ///< Deprecated log category
        lc_transfer, ///< Data transfer specific messages
        lc_backend, ///< Backend specific messages
        lc_application, ///<Application specific messages
        lc_none ///< Miscellaneous messages, not fitting in any other category
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
             * \param category Log-category of the message.
             * \param messag Message to be logged.
             */
            LogMessage(const LogCategory category, const std::string & message);

            /// Destructor.
            ~LogMessage();
    };

    /// Log a message.
    static inline void log(LogCategory category, const std::string & message)
    {
        LogMessage(category, message);
    }

/**
 * \def LOGMESSAGE
 *
 * \brief Convenience definition that provides a way to declare uniquely-named
 * instances of class LogMessage.
 *
 * The created LogMessage will be automatically enqueued with the LogQueue and written.
 *
 * \param l Log category of the message.
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
