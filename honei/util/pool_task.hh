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

#ifndef LIBUTIL_GUARD_POOL_TASK_HH
#define LIBUTIL_GUARD_POOL_TASK_HH 1

#include <honei/util/assertion.hh>
#include <honei/util/condition_variable.hh>
#include <honei/util/instantiation_policy-impl.hh>
#include <honei/util/lock.hh>
#include <honei/util/mutex.hh>
#include <honei/util/pool_task.hh>

#include <tr1/functional>

namespace honei
{
    typedef std::tr1::function<void () throw ()> WorkerTask;

    class PoolTask :
        public InstantiationPolicy<PoolTask, NonCopyable>
    {
        private:
            /// Our task.
            WorkerTask const _task;

            /// Our mutex.
            Mutex * const _mutex;

            /// Our condition variable.
            ConditionVariable * _task_finished;

            /// Indicates whether our task is done yet.
            bool _task_finished_flag;

        public:
            PoolTask() :
                _mutex(new Mutex)
            {
                CONTEXT("When creating PoolTask:");
            }

            PoolTask(WorkerTask & task) :
                _task(task),
                _mutex(new Mutex),
                _task_finished(new ConditionVariable),
                _task_finished_flag(false)
            {
                CONTEXT("When creating PoolTask:");
            }

            ///Destructor. (todo)
            ~PoolTask()
            {
                CONTEXT("When destroying PoolTask:");

                delete _mutex;
                delete _task_finished;
            }

            void operator() ()
            {
                CONTEXT("In PoolTask, when executing operator()");

                Lock l(*_mutex);

                _task();
                _task_finished_flag = true;
                _task_finished->broadcast();
            }

            void wait_on()
            {
                CONTEXT("In PoolTask, when performing wait_on():");

                Lock l(*_mutex);
                if (! _task_finished_flag)
                    _task_finished->wait(*_mutex);
            }
    };
}

#endif
