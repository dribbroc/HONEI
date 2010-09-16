/* vim: set sw=4 sts=4 et foldmethod=syntax : */
/*
 * Copyright (c) 2009, 2010 Sven Mallach <mallach@honei.org>
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

#ifndef MULTICORE_GUARD_THREAD_TASK_HH
#define MULTICORE_GUARD_THREAD_TASK_HH 1

#include <honei/backends/multicore/ticket.hh>

#include <tr1/functional>

namespace honei
{
    namespace mc
    {
        struct ThreadTask
        {
            typedef std::tr1::function<void () throw ()> WorkFunctor;

            WorkFunctor * const functor;
            Ticket<tags::CPU::MultiCore> * const ticket;

            template <typename WorkerTask> ThreadTask(WorkerTask & task, Ticket<tags::CPU::MultiCore> * const tick) :
                functor(new WorkFunctor(task)),
                ticket(tick)
            {
            }

            ~ThreadTask()
            {
                delete functor;
            }
        };

        struct TaskComp
        {
            const unsigned sched_lpu;

            TaskComp(unsigned slpu) :
                sched_lpu(slpu)
            {
            }

            bool operator () (ThreadTask * const t) const
            {
                const unsigned sched_min = t->ticket->sid_min();
                const unsigned sched_max = t->ticket->sid_max();

                if (sched_min == 0xFFFF || (sched_min <= sched_lpu && sched_lpu <= sched_max))
                    return true;

                return false;
            }
        };
    }
}
#endif
