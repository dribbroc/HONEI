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

#ifndef CUDA_GUARD_GPU_POOL_HH
#define CUDA_GUARD_GPU_POOL_HH 1

#include <honei/backends/cuda/gpu_function.hh>
#include <honei/util/attributes.hh>
#include <honei/util/instantiation_policy.hh>
#include <honei/util/thread.hh>

#include <vector>
#include <queue>

namespace honei
{
    namespace cuda
    {
        class GPUPool :
            public InstantiationPolicy<GPUPool, Singleton>
        {
            private:

                // Number of GPUs in use
                unsigned num_gpus;

                // List of user POSIX threads
                std::vector<std::pair<Thread *, GPUFunction *> > threads HONEI_ALIGNED(128);

                // Waiting queue of worker tasks to be executed
                std::queue<GPUTask *> tasks_gpu0 HONEI_ALIGNED(128);
                std::queue<GPUTask *> tasks_gpu1 HONEI_ALIGNED(128);

                // Our Mutex
                Mutex * const mutex;

                // Condition Variable used to synchronize all threads
                ConditionVariable * const global_barrier;


            public:
                GPUPool();

                ~GPUPool();

                unsigned get_num_gpus() const;

                Ticket<tags::GPU::MultiCore> * enqueue(const std::tr1::function<void ()> & task, int device);

                bool idle();
        };
    }
}
#endif
