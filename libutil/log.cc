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

#include <libutil/condition_variable.hh>
#include <libutil/log.hh>
#include <libutil/lock.hh>
#include <libutil/exception.hh>
#include <libutil/stringify.hh>
#include <libutil/thread.hh>

#include <iostream>
#include <list>
#include <string>
#include <tr1/functional>

#include <syscall.h>

namespace honei
{
    /**
     * LogQueue prints LogMessages using a single thread.
     */
    struct LogQueue
    {
        /// Our mutex.
        Mutex * const mutex;

        /// Our work-pending condition variable.
        ConditionVariable * const work_pending;

        /// Our finish-the-thread flag.
        bool finish;

        /// Our logging thread.
        Thread * thread;

        /// Our list of log messages.
        std::list<LogMessage> messages;

        /// Write out our log messages.
        void log_function()
        {
            static std::string previous_context;

            while (true)
            {
                LogMessage message;

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
                        message = messages.front();
                        messages.pop_front();
                    }
                }

                if (previous_context == message.context())
                    std::cerr << "(same context) " << message.message() << std::endl;
                else
                    std::cerr << message.context() << message.message() << std::endl;
                previous_context = message.context();
            }
        }

        LogQueue() :
            mutex(new Mutex),
            work_pending(new ConditionVariable),
            finish(false),
            thread(new Thread(std::tr1::bind(std::tr1::mem_fn(&LogQueue::log_function), this)))
        {
        }

        ~LogQueue()
        {
            {
                Lock l(*mutex);

                finish = true;
                work_pending->signal();
            }

            delete thread;
            delete work_pending;
            delete mutex;
        }

        void enqueue(const LogMessage & message)
        {
            Lock l(*mutex);

            messages.push_back(message);
            work_pending->signal();
        }
    };

    namespace intern
    {
        static LogQueue log_queue;
    }

    LogMessage::LogMessage(const LogLevel level, const std::string & message) :
        _context("In thread ID '" + stringify(syscall(SYS_gettid)) + "':\n ... " + Context::backtrace("\n ... ")),
        _level(level),
        _message(message)
    {
        intern::log_queue.enqueue(*this);
    }

    LogMessage::LogMessage()
    {
    }
}
