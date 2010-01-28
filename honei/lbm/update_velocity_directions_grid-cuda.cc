/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c)  2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/lbm/update_velocity_directions_grid.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/backends/cuda/gpu_pool.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/util/configuration.hh>


using namespace honei;

namespace
{
    class cudaUpVelDirGridfloat
    {
        private:
            PackedGridInfo<D2Q9> & info;
            PackedGridData<D2Q9, float> & data;
            unsigned long blocksize;
        public:
            cudaUpVelDirGridfloat(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data, unsigned long blocksize) :
                info(info),
                data(data),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                info.limits->lock(lm_read_only);

                void * cuda_types_gpu(info.cuda_types->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * f_temp_1_gpu(data.f_temp_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_2_gpu(data.f_temp_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_3_gpu(data.f_temp_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_4_gpu(data.f_temp_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_5_gpu(data.f_temp_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_6_gpu(data.f_temp_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_7_gpu(data.f_temp_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_8_gpu(data.f_temp_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

                unsigned long start((*info.limits)[0]);
                unsigned long end((*info.limits)[info.limits->size() - 1]);
                cuda_up_vel_dir_grid_float(start, end, cuda_types_gpu,
                        f_temp_1_gpu, f_temp_2_gpu,
                        f_temp_3_gpu, f_temp_4_gpu, f_temp_5_gpu,
                        f_temp_6_gpu, f_temp_7_gpu, f_temp_8_gpu,
                        blocksize);

                info.cuda_types->unlock(lm_read_only);

                data.f_temp_1->unlock(lm_read_and_write);
                data.f_temp_2->unlock(lm_read_and_write);
                data.f_temp_3->unlock(lm_read_and_write);
                data.f_temp_4->unlock(lm_read_and_write);
                data.f_temp_5->unlock(lm_read_and_write);
                data.f_temp_6->unlock(lm_read_and_write);
                data.f_temp_7->unlock(lm_read_and_write);
                data.f_temp_8->unlock(lm_read_and_write);
            }
    };
}

void UpdateVelocityDirectionsGrid<tags::GPU::CUDA, lbm_boundary_types::NOSLIP>::value(
        PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data)
{
    CONTEXT("When updating velocity directions (CUDA):");


    unsigned long blocksize(Configuration::instance()->get_value("cuda::up_vel_dir_grid_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaUpVelDirGridfloat task(info, data, blocksize);
        task();
    }
    else
    {
        cudaUpVelDirGridfloat task(info, data, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0)->wait();
    }
}

