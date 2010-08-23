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

#pragma once
#ifndef CUDA_GUARD_GPU_FUNCTION_HH
#define CUDA_GUARD_GPU_FUNCTION_HH 1

#include <honei/util/condition_variable.hh>
#include <honei/util/private_implementation_pattern.hh>
#include <honei/backends/cuda/ticket.hh>

#include <queue>
#include <tr1/functional>

namespace honei
{
    namespace cuda
    {
        struct GPUTask
        {
            typedef std::tr1::function<void () throw ()> WorkFunctor;

            WorkFunctor * functor;
            Ticket<tags::GPU::MultiCore> * ticket;

            template <typename WorkerTask> GPUTask(WorkerTask & task, Ticket<tags::GPU::MultiCore> * tick) :
                functor(new WorkFunctor(task)),
                ticket(tick)
            {
            }

            ~GPUTask()
            {
                delete functor;
            }
        };

        class GPUFunction :
            public PrivateImplementationPattern<GPUFunction, Shared>
        {
            private:

            public:

                GPUFunction(int id, Mutex * const mutex, ConditionVariable * const barrier, std::queue<GPUTask *> * const queue);

                ~GPUFunction();

                /// The threads' main function
                void operator() ();

                void stop();

                bool idle();
        };
    }
}
#endif
