/* vim: set sw=4 sts=4 et foldmethod=syntax : */
/*
 * Copyright (c) 2009, 2010, 2011 Sven Mallach <mallach@honei.org>
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
#include <honei/util/tr1_boost.hh>

namespace honei
{
    namespace mc
    {
        struct ThreadTask
        {
            const function<void ()> * functor;
            Ticket<tags::CPU::MultiCore> ticket;

            ThreadTask(const function<void ()> & task, Ticket<tags::CPU::MultiCore> & tick) :
                functor(new function<void ()>(task)),
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
                const unsigned sched_min = t->ticket.sid_min();
                const unsigned sched_max = t->ticket.sid_max();

                if (sched_min == 0xFFFF || (sched_min <= sched_lpu && sched_lpu <= sched_max))
                    return true;

                return false;
            }
        };

        struct PoolSyncData
        {
            Mutex * const mutex;
            ConditionVariable * const barrier;
            Mutex * const steal_mutex; // Currently only used with work stealing

            PoolSyncData() :
                mutex(new Mutex),
                barrier(new ConditionVariable),
                steal_mutex(new Mutex)
            {
            }

            ~PoolSyncData()
            {
                delete steal_mutex;
                delete barrier;
                delete mutex;
            }
        };

        struct ThreadData
        {
            volatile bool terminate;

            // Mutex for making the local task list secure
            // currently only used with work stealing with std::deque
            // as a work around
            Mutex * const local_mutex;

            ThreadData() :
                terminate(false),
                local_mutex(new Mutex)
            {
            }

            ~ThreadData()
            {
                delete local_mutex;
            }
        };
    }
}
#endif
