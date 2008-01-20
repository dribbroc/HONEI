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

#include <honei/libutil/condition_variable.hh>
#include <honei/libutil/log.hh>
#include <honei/libutil/lock.hh>
#include <honei/libutil/exception.hh>
#include <honei/libutil/stringify.hh>
#include <honei/libutil/thread.hh>
#include <honei/libutil/assertion.hh>

#include <iostream>
#include <list>
#include <string>
#include <tr1/functional>

#include <syscall.h>

namespace honei
{
    namespace intern
    {
        struct LogData
        {
            /// Our finish-the-thread flag.
            bool finish_message;

            /// Our context.
            std::string context;

            /// Our log level.
            LogLevel level;

            /// Our message.
            std::string message;

            /// Constructor.
            LogData(const std::string & c, LogLevel l, const std::string & m) :
                finish_message(false),
                context(c),
                level(l),
                message(m)
            {
            }

            /// Constructor.
            LogData() :
                finish_message(true),
                context(""),
                level(ll_minimal),
                message("")
            {
            }
        };
    }

    /**
     * LogQueue prints LogMessages using a single thread.
     */
    struct LogQueue
    {
        /// Our mutex.
        Mutex * const mutex;

        /// Our work-pending condition variable.
        ConditionVariable * const work_pending;

        /// Our logging thread.
        Thread * thread;

        /// Our list of log messages.
        std::list<intern::LogData *> messages;

        /// Our previous context.
        std::string previous_context;

        /// Write out our log messages.
        void log_function()
        {
            bool finish(false);

            while (true)
            {
                intern::LogData * data(0);

                {
                    Lock l(*mutex);

                    if (messages.empty())
                    {
                        if (finish)
                            break;

                        work_pending->wait(*mutex);
                        continue;
                    }
                    else
                    {
                        ASSERT(! messages.empty(), "messages should not be empty!");

                        data = messages.front();
                        messages.pop_front();

                        finish |= data->finish_message;
                        if (data->finish_message)
                        {
                            delete data;
                            data = 0;
                            continue;
                        }
                    }
                }

                if (previous_context == data->context)
                    std::cerr << "(same context) " << data->message << std::endl;
                else
                    std::cerr << data->context << data->message << std::endl;
                previous_context = data->context;

                delete data;
            }
        }

        LogQueue() :
            mutex(new Mutex),
            work_pending(new ConditionVariable),
            previous_context("(none)")
        {
            thread = new Thread(std::tr1::bind(std::tr1::mem_fn(&LogQueue::log_function), this));
        }

        ~LogQueue()
        {
            enqueue(new intern::LogData);

            delete thread;
            delete work_pending;
            delete mutex;
        }

        void enqueue(intern::LogData * data)
        {
            Lock l(*mutex);

            messages.push_back(data);
            work_pending->signal();
        }
    };

    namespace intern
    {
        LogQueue log_queue;
    }

    LogMessage::LogMessage(const LogLevel level, const std::string & message)
    {
        intern::log_queue.enqueue(new intern::LogData(
                    "In thread ID '" + stringify(syscall(SYS_gettid)) + "':\n ... " + Context::backtrace("\n ... "),
                    level, message));
    }
}
