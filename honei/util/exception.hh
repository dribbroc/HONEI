/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * Based upon 'exception.hh' from Paludis, which is:
 *     Copyright (c) 2005, 2006, 2007 Ciaran McCreesh
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

#ifndef LIBUTIL_GUARD_EXCEPTION_HH
#define LIBUTIL_GUARD_EXCEPTION_HH 1

#include <honei/util/instantiation_policy.hh>

#include <string>

namespace honei
{
    /**
     * Backtrace class context.
     *
     * \ingroup grpexceptions
     * \nosubgrouping
     */
    class Context :
        public InstantiationPolicy<Context, NonCopyable>
    {
        public:
            /**
             * Constructor.
             *
             * \param file Name of the source file that contains the context.
             * \param line Line number of the context.
             * \param context A sentence that describes this context.
             */
            Context(const char * const file, const long line, const std::string & context);

            /// Desctructor.
            ~Context();

            /**
             * Current context (forwards to libebt).
             */
            static std::string backtrace(const std::string & delimiter);
    };

/**
 * \def CONTEXT
 *
 * \brief Convenience definition that provides a way to declare uniquely-named
 * instances of class Context.
 *
 * The created Context will be automatically provided with the correct filename and
 * line number.
 *
 * \param s Context message that can be display by an exception-triggered backtrace.
 *
 * \warning Will only be compiled in when debug support is enabled.
 *
 * \ingroup grpdebug
 */
#if defined (DEBUG)
// C preprocessor abomination following...
#define CONTEXT_NAME_(x) ctx_##x
#define CONTEXT_NAME(x) CONTEXT_NAME_(x)
#define CONTEXT(s) \
    Context CONTEXT_NAME(__LINE__)(__FILE__, __LINE__, (s))
#else
#define CONTEXT(s)
#endif

    /**
     * Base exception class.
     *
     * \ingroup grpexceptions
     * \nosubgrouping
     */
    class Exception :
        public std::exception
    {
        private:
            struct ContextData;

            /// Our (local) context data.
            ContextData * const _context_data;

            /// Our message.
            const std::string _message;

            /// Our what string (for std::exception).
            mutable std::string _what_str;

        protected:
            /**
             * Constructor.
             *
             * \param message The exception's message.
             */
            Exception(const std::string & message) throw ();

            /// Copy-constructor.
            Exception(const Exception & e);

        public:
            /// Destructor.
            virtual ~Exception() throw ();

            /// Return a descriptive error message.
            const std::string & message() const throw ();

            /// Return a backtrace.
            std::string backtrace(const std::string & delimiter) const;

            /// Return true if the backtrace is empty.
            bool empty() const;

            /// Return a descriptive exception name.
            const char * what() const throw ();
    };

    /**
     * InternalError is an Exception that is thrown if something that is
     * never supposed to happen happens.
     *
     * \ingroup grpexceptions
     * \nosubgrouping
     */
    class InternalError :
        public Exception
    {
        public:
            /**
             * Constructor.
             *
             * \param message A short error message.
             */
            InternalError(const std::string & message) throw ();
    };

    /**
     * ExternalError is an Exception that is thrown if an external library
     * errors out.
     */
    class ExternalError :
        public Exception
    {
        public:
            /**
             * Constructor.
             *
             * \param library The name of the external library.
             * \param message A short error message.
             */
            ExternalError(const std::string & library, const std::string & message) throw ();
    };

    /**
     * PThreadError is an Exception that is thrown if one of libpthread's functions
     * fails.
     */
    class PThreadError :
        public ExternalError
    {
        public:
            /**
             * Constructor.
             *
             * \param function The name of the function that failed.
             * \param errno The error number or return value of the failed functions.
             */
            PThreadError(const std::string & function, int errno) throw ();
    };

    /**
     * SPEError is thrown by SPEManager and related classes whenever an error
     * occurs in interfacing Libspe2.
     *
     * \ingroup grpexceptions
     * \ingroup grpcell
     */
    struct SPEError :
        public ExternalError
    {
        /**
         * Constructor.
         *
         * \param msg The error message.
         * \param reason The reason for the error message.
         */
        SPEError(const std::string & function, int errno) throw ();
    };
}

#endif
