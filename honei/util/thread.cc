/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * Based upon 'thread.cc' from Paludis, which is:
 *     Copyright (c) 2007 Ciaran McCreesh
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

#include <honei/util/exception.hh>
#include <honei/util/instantiation_policy-impl.hh>
#include <honei/util/lock.hh>
#include <honei/util/mutex.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <honei/util/stringify.hh>
#include <honei/util/thread.hh>
#include <honei/util/log.hh>

#include <cstring>
#include <cerrno>

#include <pthread.h>

namespace honei
{
    template <> struct Implementation<Thread>
    {
        /// Our thread function.
        const Thread::Function function;

        /// Our thread.
        pthread_t * const thread;

        /// Our mutex.
        Mutex * const mutex;

        /// Our completion state.
        bool completed;

        static void * thread_function(void * argument)
        {
            CONTEXT("When runing libutil-thread");
            Implementation * imp(static_cast<Implementation *>(argument));

            /// \todo Implement exception handling for the call to function.
            try
            {
                imp->function();
            }
            catch (Exception & e)
            {
                throw InternalError("Exception in Thread: " + stringify(e.what()));
            }
            catch (...)
            {
                LOGMESSAGE(ll_minimal, "Unexpected std::exception or similar in Thread!");
            }

            Lock l(*imp->mutex);
            imp->completed = true;

            pthread_exit(0);
        }

        Implementation(const Thread::Function & f) :
            function(f),
            thread(new pthread_t),
            mutex(new Mutex),
            completed(false)
        {
            int retval;

            if (0 != (retval = pthread_create(thread, 0, &thread_function, this)))
                throw ExternalError("libpthread", "pthread_create failed, " + stringify(strerror(retval)));
        }

        ~Implementation()
        {
            pthread_join(*thread, 0);

            delete thread;
            delete mutex;
        }
    };
}

using namespace honei;

Thread::Thread(const Function & function) :
    PrivateImplementationPattern<Thread, Single>(new Implementation<Thread>(function))
{
}

Thread::~Thread()
{
}

bool Thread::completed() const
{
    Lock l(*_imp->mutex);

    return _imp->completed;
}
