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

#include <honei/backends/multicore/thread.hh>
#include <honei/util/exception.hh>
#include <honei/util/lock.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <honei/util/stringify.hh>

#include <errno.h>
#include <pthread.h>
#include <sys/syscall.h>

namespace honei
{
    template <> struct Implementation<mc::Thread>
    {
        /// Our underlying POSIX thread.
        pthread_t * thread;

        /// Our Thread ID (given by the operating system)
        unsigned tid;

        /// The thread pool's mutex (for grabbing work).
        Mutex * const pool_mutex;

        /// The task list administrated by the thread pool.
        std::list<mc::ThreadTask *> * const tasklist;

        /// ConditionVariable for work-status.
        ConditionVariable * const global_barrier;

        /// Flag if the Thread shall stop.
        bool terminate;

        /// The current ThreadTask to execute
        mc::ThreadTask * task;

        /// Helper function to pick work out of the pool's task list
        inline void pick_work() __attribute__((always_inline))
        {
            pthread_mutex_lock(pool_mutex->mutex());

            task = 0;

            for (std::list<mc::ThreadTask *>::iterator i(tasklist->begin()) , i_end(tasklist->end()) ; i != i_end ; ++i)
            {
                unsigned & thread_id = (*i)->ticket->tid();

                if (thread_id == 0 || thread_id == tid)
                {
                    task = *i;
                    thread_id = tid;
                    tasklist->remove(*i);
                    break;
                }
            }

            pthread_mutex_unlock(pool_mutex->mutex());
        }

        /// The threads' main function
        static void * thread_fn(void * argument)
        {
            Implementation * imp(static_cast<Implementation *>(argument));

            imp->tid = syscall(__NR_gettid);

            do
            {
                // Check work-status and condition again under same mutex protection to avoid race conditions
                // between pick / if and wait
                {
                    Lock l(*imp->pool_mutex);
                    imp->pick_work();
                    if (imp->task == 0 && ! imp->terminate)
                    {
                        imp->global_barrier->wait(*imp->pool_mutex);
                        imp->pick_work();
                    }
                }

                while (imp->task != 0)
                {
                    (*imp->task->functor)();
                    imp->task->ticket->mark();
                    delete imp->task;
                    imp->pick_work();
                }
            }
            while (! imp->terminate);

            pthread_exit(0);
        }

        Implementation(Mutex * const pool_m, std::list<mc::ThreadTask *> * l, ConditionVariable * g_b) :
            thread(new pthread_t),
            global_barrier(g_b),
            pool_mutex(pool_m),
            tasklist(l),
            tid(0),
            terminate(false),
            task((mc::ThreadTask *) 1)
        {
            int retval;

            if (0 != (retval = pthread_create(thread, 0, &thread_fn, this)))
                throw ExternalError("libpthread", "pthread_create failed with return value: " + stringify(retval));

            // Wait busy until thread is really set up.
            while (task != 0)
            {
                sleep(1);
            }
        }

        ~Implementation()
        {
            pthread_mutex_lock(pool_mutex->mutex());
            terminate = true;
            global_barrier->broadcast();
            pthread_mutex_unlock(pool_mutex->mutex());

            pthread_join(*thread, 0);

            delete thread;
        }
    } __attribute__((aligned(128)));
}

using namespace honei::mc;

Thread::Thread(Mutex * const pool_mutex, std::list<ThreadTask *> * list, ConditionVariable * cv) :
    PrivateImplementationPattern<Thread, Shared>(new Implementation<Thread>(pool_mutex, list, cv))
{
}

Thread::~Thread()
{
}

const unsigned Thread::tid()
{
    return _imp->tid;
}
