/* vim: set sw=4 sts=4 et foldmethod=syntax : */
/*
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@tu-dortmund.de>
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

#include <honei/backends/cuda/gpu_function.hh>
#include <honei/util/attributes.hh>
#include <honei/util/exception.hh>
#include <honei/util/lock.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <honei/util/stringify.hh>
#include <honei/backends/cuda/multi_gpu.hh>

#include <errno.h>
#include <utility>
#include <sys/syscall.h>

namespace honei
{
    template <> struct Implementation<cuda::GPUFunction>
    {
        /// Our GPU id
        int gpu_id;

        /// The gpu pool's mutex (for grabbing work).
        Mutex * const pool_mutex;

        /// The task queue administrated by the gpu pool.
        std::queue<cuda::GPUTask *> * const taskqueue;

        /// ConditionVariable for work-status.
        ConditionVariable * const global_barrier;

        /// Flag if the Thread shall stop.
        bool terminate;

        /// The current GPUTask to execute
        cuda::GPUTask * task;

        /// Helper function to pick work out of the pool's task queue
        inline void pick_work() HONEI_INLINE
        {
            pthread_mutex_lock(pool_mutex->mutex());

            task = 0;

            if(taskqueue->size() > 0)
            {
                task = taskqueue->front();
                taskqueue->pop();
            }

            pthread_mutex_unlock(pool_mutex->mutex());
        }

        Implementation(int id, Mutex * const mutex, ConditionVariable * const barrier, std::queue<cuda::GPUTask *> * const queue) :
            gpu_id(id),
            pool_mutex(mutex),
            taskqueue(queue),
            global_barrier(barrier),
            terminate(false),
            task((cuda::GPUTask *) 1)
            {
            }

            void stop()
            {
                pthread_mutex_lock(pool_mutex->mutex());
                terminate = true;
                global_barrier->broadcast();
                pthread_mutex_unlock(pool_mutex->mutex());
            }

            ~Implementation()
            {
            }
    };
}

using namespace honei::cuda;

GPUFunction::GPUFunction(int id, Mutex * const mutex, ConditionVariable * const barrier, std::queue<GPUTask *> * const queue) :
    PrivateImplementationPattern<GPUFunction, Shared>(new Implementation<GPUFunction>(id, mutex, barrier, queue))
{
}

GPUFunction::~GPUFunction()
{
}

void GPUFunction::stop()
{
    _imp->stop();
}

void GPUFunction::operator() ()
{
    cuda_set_device(_imp->gpu_id);
    do
    {
        // Check work-status and condition again under same mutex protection to avoid race conditions
        // between pick / if and wait
        {
            Lock l(*_imp->pool_mutex);
            _imp->pick_work();
            if (_imp->task == 0 && ! _imp->terminate)
            {
                _imp->global_barrier->wait(*_imp->pool_mutex);
                _imp->pick_work();
            }
        }

        while (_imp->task != 0)
        {
            (*_imp->task->functor)();
            _imp->task->ticket->mark();
            delete _imp->task;
            _imp->pick_work();
        }
    }
    while (! _imp->terminate);
}
