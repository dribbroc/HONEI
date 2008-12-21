/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Sven Mallach <sven.mallach@cs.tu-dortmund.de>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef MULTICORE_GUARD_THREAD_HH
#define MULTICORE_GUARD_THREAD_HH 1

#include <honei/util/condition_variable.hh>
#include <honei/util/mutex.hh>
#include <honei/util/private_implementation_pattern.hh>
#include <honei/util/ticket.hh>

#include <list>
#include <tr1/functional>

namespace honei
{
    namespace mc
    {
        struct ThreadTask
        {
            typedef std::tr1::function<void () throw ()> WorkFunctor;

            WorkFunctor * functor;
            Ticket * ticket;
            unsigned * thread_id;

            template <typename WorkerTask> ThreadTask(WorkerTask & task, Ticket * tick, unsigned * tid) :
                functor(new WorkFunctor(task)),
                ticket(tick),
                thread_id(tid)
            {
            }

            ~ThreadTask()
            {
                delete functor;
            }
        };

        class Thread :
            public PrivateImplementationPattern<Thread, Shared>
        {
            private:
                /// Unwanted copy-constructor: Do not implement. See EffCpp, Item 27.
                Thread(const Thread &);

                /// Unwanted assignment operator: Do not implement. See EffCpp, Item 27.
                Thread & operator= (const Thread &);

            public:
                /// Constructor (cpu_id is the core to assign the thread to if affinity is enabled)
                Thread(Mutex * const pool_mutex, std::list<ThreadTask *> * list, ConditionVariable * cv);

                /// Destructor.
                ~Thread();

                /// Return the thread id given by the operating system
                const unsigned tid();
        };
    }
}
#endif
