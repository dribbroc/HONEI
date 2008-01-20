/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#include <honei/libutil/assertion.hh>
#include <honei/libutil/condition_variable.hh>
#include <honei/libutil/configuration.hh>
#include <honei/libutil/exception.hh>
#include <honei/libutil/lock.hh>
#include <honei/libutil/profiler.hh>
#include <honei/libutil/stringify.hh>
#include <honei/libutil/thread.hh>
#include <honei/libutil/time_stamp.hh>

#include <fstream>
#include <list>
#include <map>
#include <string>
#include <tr1/functional>

#include <syscall.h>

namespace honei
{
    namespace intern
    {
        struct ProfileId
        {
            /// Our function.
            const std::string function;

            /// Our tag.
            const std::string tag;

            /// Constructor.
            ProfileId(const std::string & f, const std::string & t) :
                function(f),
                tag(t)
            {
            }

            /// Comparison operator for less-than.
            bool operator< (const ProfileId & other) const
            {
                if (function < other.function)
                    return true;

                if (function == other.function)
                {
                    if (tag < other.tag)
                        return true;
                }

                return false;
            }
        };

        struct ProfilerData
        {
            /// Our finish-the-thread flag.
            bool finish;

            /// Our function.
            std::string function;

            /// Our tag.
            std::string tag;

            /// Our type.
            ProfilerMessageType type;

            /// Our thread id.
            int thread;

            /// Our time stamp.
            TimeStamp stamp;

            /// Our time.
            unsigned time;

            /// Constructor.
            ProfilerData(const std::string & f, const std::string & ta, ProfilerMessageType ty, int tid, TimeStamp s) :
                finish(false),
                function(f),
                tag(ta),
                type(ty),
                thread(tid),
                stamp(s)
            {
            }

            /// Constructor.
            ProfilerData() :
                finish(true),
                function(""),
                tag(""),
                type(pmt_stop),
                thread(0)
            {
            }

            /// Constructor.
            ProfilerData(const std::string & f, const std::string & ta, unsigned ti) :
                finish(false),
                function(f),
                tag(ta),
                type(pmt_direct),
                thread(0),
                time(ti)
            {
            }
        };

        struct ProfilerSession
        {
            /// Our function.
            const std::string function;

            /// Our tag.
            const std::string tag;

            /// Our thread id.
            const int thread;

            /// Our time stamp.
            TimeStamp stamp;

            /// Constructor.
            ProfilerSession(const std::string & f, const std::string & t, int tid, TimeStamp s) :
                function(f),
                tag(t),
                thread(tid),
                stamp(s)
            {
            }
        };
    }

    /**
     * Profiler evaluates ProfilerData in a separate thread.
     */
    struct Profiler
    {
        typedef std::list<intern::ProfilerData *> DataQueue;

        typedef std::list<intern::ProfilerSession> SessionList;

        typedef std::list<float> TimingList;

        typedef std::map<intern::ProfileId, TimingList> ResultMap;

        typedef std::tr1::function<void (const std::string &, const std::string &, unsigned, float, float, float)> EvaluationFunction;

        typedef std::list<EvaluationFunction> EvaluationFunctions;

        /// Our queue's mutex.
        Mutex * const queue_mutex;

        /// Our evaluation's mutex.
        Mutex * const mutex;

        /// Our work-pending condition variable.
        ConditionVariable * const work_pending;

        /// Our profiling thread.
        Thread * thread;

        /// Our list of profiler data.
        DataQueue queue;

        /// Our list of open sessions.
        SessionList sessions;

        /// Our profiling result.
        ResultMap results;

        /// Our output stream.
        std::fstream output;

        /// Our evaluation functions.
        EvaluationFunctions evaluation_functions;

        /// Process our profiler data.
        void profiler_function()
        {
            bool finish(false);

            while (true)
            {
                intern::ProfilerData * data(0);

                {
                    Lock l(*queue_mutex);

                    if (queue.empty())
                    {
                        if (finish)
                        {
                            break;
                        }

                        work_pending->wait(*queue_mutex);
                        continue;
                    }
                    else
                    {
                        ASSERT(! queue.empty(), "queue should not be empty!");

                        data = queue.front();
                        queue.pop_front();

                        if (data->finish)
                        {
                            finish = true;

                            delete data;
                            data = 0;
                            continue;
                        }
                    }
                }

                if (data->type == pmt_start)
                {
                    Lock l(*mutex);

                    sessions.push_back(intern::ProfilerSession(data->function, data->tag, data->thread, data->stamp));
                }
                else if (data->type == pmt_stop)
                {
                    Lock l(*mutex);

                    SessionList::const_iterator i(sessions.begin()), i_end(sessions.end());
                    for ( ; i != i_end ; ++i)
                    {
                        // Compare thread first, as integer comparison is fast.
                        if (i->thread != data->thread)
                            continue;

                        if (i->function != data->function)
                            continue;

                        if (i->tag != data->tag)
                            continue;

                        break;
                    }

                    if (i == i_end)
                    {
                        output << "Received profiler stop message for '" << data->function << "@" << stringify(data->thread)
                                << ":" << data->tag << "' that has no been started." << std::endl;
                    }
                    else
                    {
                        intern::ProfileId id(data->function, data->tag);
                        float timing((data->stamp.usec() - i->stamp.usec()));

                        ResultMap::iterator i(results.find(id)), i_end(results.end());
                        if (i_end == i)
                        {
                            TimingList list;
                            list.push_back(timing);
                            results.insert(ResultMap::value_type(id, list));
                        }
                        else
                        {
                            i->second.push_back(timing);
                        }
                    }
                }
                else if (data->type == pmt_direct)
                {
                    Lock l(*mutex);

                    intern::ProfileId id(data->function, data->tag);
                    float timing(data->time);

                    ResultMap::iterator i(results.find(id)), i_end(results.end());
                    if (i_end == i)
                    {
                        TimingList list;
                        list.push_back(timing);
                        results.insert(ResultMap::value_type(id, list));
                    }
                    else
                    {
                        i->second.push_back(timing);
                    }
                }
                else if (data->type == pmt_evaluate)
                {
                    Lock l(*mutex);

                    output << "--------------------------------------------" << std::endl;
                    output << "function[tag] : count highest average lowest" << std::endl;
                    output << "--------------------------------------------" << std::endl;

                    for (ResultMap::const_iterator r(results.begin()), r_end(results.end()) ; r != r_end ; ++r)
                    {
                        unsigned count(0);
                        float highest(0.0f), average(0.0f), lowest(std::numeric_limits<float>::max());
                        for (TimingList::const_iterator t(r->second.begin()), t_end(r->second.end()) ; t != t_end ; ++t)
                        {
                            ++count;
                            highest = std::max(highest, *t);
                            average += *t;
                            lowest = std::min(lowest, *t);
                        }
                        if (count > 0)
                            average /= count;

                        for (EvaluationFunctions::iterator f(evaluation_functions.begin()), f_end(evaluation_functions.end()) ; f != f_end ; ++f)
                        {
                            EvaluationFunction ef(*f);
                            ef(r->first.function, r->first.tag, count, highest, average, lowest);
                        }
                    }
                }

                delete data;
            }
        }

        void evaluation_printer(const std::string & function, const std::string & tag, unsigned count, float highest, float average, float lowest)
        {
            output << function << "[" << tag << "] : " << count << " " << highest << " " << average << " " << lowest << std::endl;
        }

        Profiler() :
            mutex(new Mutex),
            queue_mutex(new Mutex),
            work_pending(new ConditionVariable),
            output(Configuration::instance()->get_value("profiler::output", "/dev/null").c_str(), std::ios_base::out | std::ios_base::app)
        {
            output << "== New profiler session started ==" << std::endl;

            evaluation_functions.push_back(EvaluationFunction(std::tr1::bind(
                        std::tr1::mem_fn(&Profiler::evaluation_printer),
                        this, std::tr1::placeholders::_1, std::tr1::placeholders::_2,
                        std::tr1::placeholders::_3, std::tr1::placeholders::_4,
                        std::tr1::placeholders::_5, std::tr1::placeholders::_5)));
            thread = new Thread(std::tr1::bind(std::tr1::mem_fn(&Profiler::profiler_function), this));
        }

        ~Profiler()
        {
            enqueue(new intern::ProfilerData("", "", pmt_evaluate, 0, TimeStamp()));
            enqueue(new intern::ProfilerData);

            delete thread;
            delete work_pending;
            delete queue_mutex;
            delete mutex;
        }

        void enqueue(intern::ProfilerData * data)
        {
            Lock l(*queue_mutex);

            queue.push_back(data);
            work_pending->signal();
        }

        void enqueue(const EvaluationFunction & eval)
        {
            Lock l(*mutex);

            evaluation_functions.push_back(eval);
        }
    };

    namespace intern
    {
        Profiler profiler;
    }

    ProfilerMessage::ProfilerMessage(const char * function, const std::string & tag, ProfilerMessageType type)
    {
        intern::profiler.enqueue(new intern::ProfilerData(std::string(function), tag, type, syscall(SYS_gettid), TimeStamp().take()));
    }

    ProfilerMessage::ProfilerMessage(const EvaluationFunction & eval)
    {
        intern::profiler.enqueue(eval);
    }

    ProfilerMessage::ProfilerMessage(const std::string & function, const std::string & tag, unsigned time)
    {
        intern::profiler.enqueue(new intern::ProfilerData(function, tag, time));
    }
}
