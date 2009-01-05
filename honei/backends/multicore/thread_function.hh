/* vim: set sw=4 sts=4 et foldmethod=syntax : */
/*
 * Copyright (c) 2009 Sven Mallach <sven.mallach@cs.uni-dortmund.de>
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

#ifndef MULTICORE_GUARD_THREAD_OBJECT_HH
#define MULTICORE_GUARD_THREAD_OBJECT_HH 1

#include <honei/backends/multicore/ticket.hh>
#include <honei/util/condition_variable.hh>
#include <honei/util/private_implementation_pattern.hh>

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
            std::tr1::shared_ptr<Ticket<tags::CPU::MultiCore> > & ticket;

            template <typename WorkerTask> ThreadTask(WorkerTask & task, std::tr1::shared_ptr<Ticket<tags::CPU::MultiCore> > & tick) :
                functor(new WorkFunctor(task)),
                ticket(tick)
            {
            }

            ~ThreadTask()
            {
                delete functor;
            }
        };

        class ThreadFunction :
            public PrivateImplementationPattern<ThreadFunction, Shared>
        {
            private:

            public:

                ThreadFunction(Mutex * const mutex, ConditionVariable * const barrier, std::list<ThreadTask *> * const list);

                ~ThreadFunction();

                /// The threads' main function
                void operator() ();

                void stop();

                const unsigned tid() const;
        };
    }
}
#endif
