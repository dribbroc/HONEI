/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * Based in part upon 'log.cc' from Paludis, which is:
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

#include <honei/util/condition_variable.hh>
#include <honei/util/instantiation_policy-impl.hh>
#include <honei/util/log.hh>
#include <honei/util/lock.hh>
#include <honei/util/exception.hh>
#include <honei/util/stringify.hh>
#include <honei/util/thread.hh>
#include <honei/util/assertion.hh>
#include <honei/util/configuration.hh>

#include <iostream>
#include <fstream>
#include <list>
#include <memory>
#include <set>
#include <string>

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

            /// Our log category.
            LogCategory category;

            /// Our message.
            std::string message;

            /// Constructor.
            LogData(const std::string & c, LogCategory l, const std::string & m) :
                finish_message(false),
                context(c),
                category(l),
                message(m)
            {
            }

            /// Constructor.
            LogData() :
                finish_message(true),
                context(""),
                category(lc_none),
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

        /// Our category selections
        std::set<LogCategory> selections;

        /// Our output stream
        std::ostream output;

        /// Our output file stream, if present
        std::ofstream output_file;

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

                if (selections.size() == 0)
                {
                    if (data->context != "" && previous_context == data->context)
                        output << "(same context) " << data->message << std::endl;
                    else if(data->context == "")
                        output << data->message << std::endl;
                    else
                        output << data->context << std::endl << data->message << std::endl;
                    previous_context = data->context;
                }
                else
                {
                    if (selections.count(data->category) == 1)
                    {
                        if (data->context != "" && previous_context == data->context)
                            output << "(same context) " << data->message << std::endl;
                        else if(data->context == "")
                            output << data->message << std::endl;
                        else
                            output << data->context << std::endl << data->message << std::endl;
                        previous_context = data->context;
                    }
                }

                delete data;
            }
        }

        LogQueue() :
            mutex(new Mutex),
            work_pending(new ConditionVariable),
            previous_context("(none)"),
            output(std::cout.rdbuf())
        {
            std::string output_string(Configuration::instance()->get_value("log::output", "cerr"));
            if (output_string == "cerr")
            {
                std::streambuf * buf;
                buf = std::cerr.rdbuf();
                output.rdbuf(buf);
            }
            else if (output_string == "cout")
            {
                std::streambuf * buf;
                buf = std::cout.rdbuf();
                output.rdbuf(buf);
            }
            else
            {
                output_file.open(output_string.c_str(), std::ios_base::out | std::ios_base::app);
                if(!output_file.is_open())
                    throw InternalError("Log: Cannot open log output file " + output_string);
                std::streambuf * buf;
                buf = output_file.rdbuf();
                output.rdbuf(buf);
            }

            std::string config_string(Configuration::instance()->get_value("log::categories", "none"));
            if (config_string.find("transfer", 0) != std::string::npos)
            {
                selections.insert(lc_transfer);
            }
            if (config_string.find("backend", 0) != std::string::npos)
            {
                selections.insert(lc_backend);
            }
            if (config_string.find("application", 0) != std::string::npos)
            {
                selections.insert(lc_application);
            }
            if (config_string.find("none", 0) != std::string::npos)
            {
                selections.insert(lc_none);
            }

            thread = new Thread(std::bind(std::mem_fn(&LogQueue::log_function), this));
        }

        ~LogQueue()
        {
            enqueue(new intern::LogData);

            if(output_file.is_open())
                    output_file.close();

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

    LogMessage::LogMessage(const LogCategory category, const std::string & message)
    {
#ifdef DEBUG
        intern::log_queue.enqueue(new intern::LogData(
                    "In thread ID '" + stringify(syscall(SYS_gettid)) + "':\n ... " + Context::backtrace("\n ... "),
                    category, message));
#else
        intern::log_queue.enqueue(new intern::LogData(
                    "",
                    category, message));
#endif
    }

    LogMessage::~LogMessage()
    {
    }
}
