/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Volker Jung <volker.jung@uni-dortmund.de>
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

#include <libutil/assertion.hh>
#include <libutil/condition_variable.hh>
#include <libutil/exception.hh>
#include <libutil/lock.hh>
#include <libutil/mutex.hh>
#include <libutil/stringify.hh>
#include <libutil/pool_task.hh>
#include <libutil/pool_thread.hh>
#include <libutil/thread_pool.hh>

#include <list>

#include <cerrno>

using namespace honei;

struct PoolThread::Implementation
{
    private:
        /// Our thread.
        pthread_t * const _thread;

        /// The PoolThread we belong to.
        PoolThread * const _pool_thread;

        /// Our mutex for read/write access to a PoolThread instance.
        Mutex * const _mutex;

        /// Our condition variable.
        ConditionVariable * const _wake_up;

        /// Our exit condition.
        bool _exit;

        /// Our actual task.
        PoolTask * _task;

        bool _task_has_run;

        /// The pool we belong to.
        ThreadPool * const _pool;

    public:
        /// Set task and wake up this thread for computation.
        inline void run(PoolTask * task)
        {
            CONTEXT("In PoolThread, when performing run(PoolTask):");

            Lock l(*_mutex);

            _task = task;
            _task_has_run = false;

            _wake_up->signal(); //Notify our thread that new work is waiting.
        }

        /// Our thread's main function.
        static void * thread_function(void * argument)
        {
            CONTEXT("In PoolThread::thread_function :");

            Implementation * imp(static_cast<Implementation *>(argument));
            bool exit(false);
            do
            {
                Lock l(*imp->_mutex);

                // If there is no task to be run, wait
                if (imp->_task_has_run)
                {
                    imp->_wake_up->wait(*imp->_mutex);
                }

                // Update our exit condition
                exit = imp->_exit;

                // Run our task
                if (! imp->_task_has_run)
                {
                    (*imp->_task)();
                    imp->_task_has_run = true;
                    imp->_pool->notify(imp->_pool_thread);
                }
            }
            while (! exit);

            pthread_exit(0);
        }

        Implementation(ThreadPool * pool, PoolThread * pool_thread) :
            _thread(new pthread_t),
            _pool_thread(pool_thread),
            _mutex(new Mutex),
            _exit(false),
            _pool(pool),
            _task_has_run(true),
            _wake_up(new ConditionVariable),
            _task(0)
        {
            CONTEXT("When creating PoolThread::Implementation :");

            int retval;
            if (0 != (retval = pthread_create(_thread, 0, &thread_function, this)))
                throw ExternalError("libpthread", "pthread_create failed, " + stringify(strerror(retval)));
        }

        ~Implementation()
        {
            CONTEXT("When destroying PoolThrad::Implementation :");

            // Flag for exit.
            {
                Lock l(*_mutex);
                _exit = true;
                _wake_up->signal();
            }

            pthread_join(*_thread, 0);

            delete _mutex;
            delete _thread;
        }
};

PoolThread::PoolThread(ThreadPool * pool) :
    _imp(new Implementation(pool, this))
{
    CONTEXT("When creating PoolThread :");
}

PoolThread::~PoolThread()
{
    CONTEXT("When destroying PoolThread :");
    delete _imp;
}

void
PoolThread::run(PoolTask * task)
{
    CONTEXT("In PoolThread::run(PoolTask) :");
    _imp->run(task);
}
