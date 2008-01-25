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

#ifndef LIBUTIL_GUARD_PROFILER_HH
#define LIBUTIL_GUARD_PROFILER_HH 1

#include <string>
#include <tr1/functional>

namespace honei
{
    class ProfilerQueue;

    /**
     * ProfilerMessageType distinguishes between the various types of ProfilerMessage.
     *
     * \ingroup grpprofiler
     */
    enum ProfilerMessageType
    {
        pmt_start, //< A message that starts a profiling session.
        pmt_stop, //< A message that stops a profiling seesion and records it.
        pmt_direct, //< A message that directly inject a result, e.g. from the SPU.
        pmt_evaluate //< A message that triggers evaluation of all recorded sessions.
    };

    /**
     * ProfilerMessage enqueues a profiling message with the ProfilerQueue.
     *
     * \ingroup grpprofiler
     */
    class ProfilerMessage
    {
        private:
            /// \name Unwanted operations
            /// \{

            /// Unwanted default constructor: Do not implement. See EffC++, Item 27.
            ProfilerMessage();

            /// Unwanted assignment operator: Do not implement. See EffC++, Item 27.
            ProfilerMessage & operator= (ProfilerMessage & other);

            /// \}

        public:
            friend struct LogQueue;

            typedef std::tr1::function<void (const std::string &, const std::string, unsigned, float, float, float)> EvaluationFunction;

            /**
             * Constructor.
             *
             * \param function Name of the function that shall be profiled.
             * \param tag Named tag for a particular profiling purpose.
             * \param type Whether profiling was started or stopped.
             */
            ProfilerMessage(const char * function, const std::string & tag, ProfilerMessageType type);

            /**
             * Constructor.
             *
             * \param eval Evaluation-function object.
             */
            ProfilerMessage(const EvaluationFunction & eval);

            /**
             * Constructor.
             *
             * \param function Name of the function that was profiled.
             * \param tag Named tag for a particular profiling purpose.
             * \param time Time that the named session took to execute.
             */
            ProfilerMessage(const std::string & function, const std::string & tag, unsigned time);
    };

    /// Profile a function.
    static inline void profile(const char * function, const std::string & tag, ProfilerMessageType type)
    {
        ProfilerMessage(function, tag, type);
    }

/**
 * \def PROFILER_EVALUATE
 *
 * \brief Convenience definition that provides a way to evalute the profilers
 * data.
 *
 * The created ProfilerMessage will be automatically enqueue with the Profiler.
 *
 * \warning Will only be compiled in when profiler support is enabled.
 *
 * \ingroup grpprofiler
 */
#if defined (HONEI_PROFILER)
#define PROFILER_EVALUATE profile("", "", honei::pmt_evaluate)
#else
#define PROFILER_EVALUATE
#endif

/**
 * \def PROFILER_START
 *
 * \brief Convenience definition that provides a way to declare uniquely-named
 * instances of class ProfilerMessage that starts profiling of a given function.
 *
 * The created ProfilerMessage will be automatically enqueued with the Profiler.
 *
 * \param t The tag.
 *
 * \warning Will only be compiled in when profiler support is enabled.
 *
 * \ingroup grpprofiler
 */
#if defined (HONEI_PROFILER)
#define PROFILER_START(t) profile(__FUNCTION__, t, honei::pmt_start)
#else
#define PROFILER_START(t)
#endif

/**
 * \def PROFILER_STOP
 *
 * \brief Convenience definition that provides a way to declare uniquely-named
 * instances of class ProfilerMessage that stops profiling of a given function.
 *
 * The created ProfilerMessage will be automatically enqueued with the Profiler.
 *
 * \param t The tag.
 *
 * \warning Will only be compiled in when profiler support is enabled.
 *
 * \ingroup grpprofiler
 */
#if defined (HONEI_PROFILER)
#define PROFILER_STOP(t) profile(__FUNCTION__, t, honei::pmt_stop)
#else
#define PROFILER_STOP(t)
#endif
}

#endif
