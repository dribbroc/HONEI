/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Thorsten Deinert <thorsten.deinert@uni-dortmund.de>
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

#ifndef LIBUTIL_GUARD_THREAD_POOL_HH
#define LIBUTIL_GUARD_THREAD_POOL_HH 1

#include <libutil/assertion.hh>
#include <libutil/condition_variable.hh>
#include <libutil/lock.hh>
#include <libutil/mutex.hh>
#include <libutil/pool_task.hh>
#include <libutil/pool_thread.hh>
#include <libutil/stringify.hh>

#include <list>
#include <tr1/functional>

#define DEFAULT_NUM_THREADS 4

namespace honei
{
    typedef std::tr1::function<void () throw ()> WorkerTask;

    class ThreadPool
    {
        private:
            ///Our threads.
            std::list<PoolThread *> _pool;

            ///Our list of idle threads.
            std::list<PoolThread *> _idle_threads;

            ///Our list of untreated tasks.
            std::list<PoolTask *> _task_list;

            ///Our mutex.
            Mutex * const _mutex;

            ///Our condition variable.
            ConditionVariable * const _available;

            /**
             * Constructor.
             *
             * \param num_threads Number of threads the threadpool consists of.
             */
            ThreadPool(unsigned long num_threads) :
                _mutex(new Mutex),
                _available(new ConditionVariable)
            {
                CONTEXT("When creating ThreadPool:");

                for (unsigned long i(0); i < num_threads; ++i)
                {
                    PoolThread * poolthread(new PoolThread(this));
                    _pool.push_back(poolthread);
                    _idle_threads.push_back(poolthread);
                }
            }


        public:

            /// Returns the singleton instance of the pool.
            //If it doesn't exist yet, it is created with numThreads threads. numThreads is ignored otherwise.
            static ThreadPool * get_instance(unsigned long num_threads)
            {
                CONTEXT("When getting ThreadPool:");

                static ThreadPool result(num_threads);
                return &result;
            }

            /// Returns the singleton instance of the pool.
            //If it doesn't exist yet, it is created with a default number of threads. Useful for convenience.
            static ThreadPool * get_instance()
            {
                CONTEXT("When getting ThreadPool:");

                static ThreadPool result(DEFAULT_NUM_THREADS);
                return &result;
            }

           ///Destructor.
            ~ThreadPool()
            {
                CONTEXT("When destroying ThreadPool:");

                delete _mutex;
                delete _available;
            }

            ///Processes incoming tasks.
            PoolTask * dispatch(WorkerTask task)
            {
                CONTEXT("When dispatching in ThreadPool:");

                Lock l(*_mutex);

                PoolTask * pt(new PoolTask(task));
                if(! _idle_threads.empty())
                {
                    //Fetch idle thread and dispatch task.
                    PoolThread * result = _idle_threads.front();
                    _idle_threads.pop_front();
                    result->run(pt);
                }
                else
                {
                    //Since there is no idle thread, enqueue task.
                    _task_list.push_back(pt);
                }
                return pt;
            }

            void notify(PoolThread * returning_thread)
            {
                CONTEXT("When notifying ThreadPool:");

                Lock l(*_mutex);

                if(! _task_list.empty())
                {
                    //Fetch enqueued task and dispatch to returning thread.
                    PoolTask * next_task = _task_list.front();
                    _task_list.pop_front();
                    returning_thread->run(next_task);
                }
                else
                {
                    //Enqueue returning thread in our idle list.
                    _idle_threads.push_back(returning_thread);
                }
            }
    };
}
#endif
