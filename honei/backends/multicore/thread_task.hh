/* vim: set sw=4 sts=4 et foldmethod=syntax : */
/*
 * Copyright (c) 2009 - 2012 Sven Mallach <mallach@honei.org>
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

#include <honei/backends/multicore/lpu.hh>
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
            LPU * const lpu;

            TaskComp(LPU * const pu) :
                lpu(pu)
            {
            }

            bool operator () (ThreadTask * const t) const
            {
                const unsigned req_socket = t->ticket.req_socket();
                const unsigned req_sched = t->ticket.req_sched();

                if (req_sched == 0xFFFFu)
                {
                    if (req_socket == 0xFFFFu || req_socket == (unsigned) lpu->socket_id)
                        return true;
                    else
                        return false;
                }
                else if (req_sched == (unsigned) lpu->sched_id)
                {
                    return true;
                }
                else
                {
                    return false;
                }
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
