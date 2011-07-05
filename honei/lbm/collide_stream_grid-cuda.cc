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

#include <honei/lbm/collide_stream_grid.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/backends/cuda/gpu_pool.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/util/configuration.hh>


using namespace honei;

namespace
{
    class cudaColStreamGridfloat
    {
        private:
            PackedGridInfo<D2Q9> & info;
            PackedGridData<D2Q9, float> & data;
            float tau;
            unsigned long blocksize;
        public:
            cudaColStreamGridfloat(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data, float tau, unsigned long blocksize) :
                info(info),
                data(data),
                tau(tau),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                info.limits->lock(lm_read_only);

                void * cuda_dir_1_gpu(info.cuda_dir_1->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * cuda_dir_2_gpu(info.cuda_dir_2->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * cuda_dir_3_gpu(info.cuda_dir_3->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * cuda_dir_4_gpu(info.cuda_dir_4->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * cuda_dir_5_gpu(info.cuda_dir_5->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * cuda_dir_6_gpu(info.cuda_dir_6->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * cuda_dir_7_gpu(info.cuda_dir_7->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * cuda_dir_8_gpu(info.cuda_dir_8->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * f_eq_0_gpu(data.f_eq_0->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_eq_1_gpu(data.f_eq_1->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_eq_2_gpu(data.f_eq_2->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_eq_3_gpu(data.f_eq_3->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_eq_4_gpu(data.f_eq_4->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_eq_5_gpu(data.f_eq_5->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_eq_6_gpu(data.f_eq_6->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_eq_7_gpu(data.f_eq_7->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_eq_8_gpu(data.f_eq_8->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * f_0_gpu(data.f_0->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_1_gpu(data.f_1->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_2_gpu(data.f_2->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_3_gpu(data.f_3->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_4_gpu(data.f_4->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_5_gpu(data.f_5->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_6_gpu(data.f_6->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_7_gpu(data.f_7->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_8_gpu(data.f_8->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * f_temp_0_gpu(data.f_temp_0->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_temp_1_gpu(data.f_temp_1->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_temp_2_gpu(data.f_temp_2->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_temp_3_gpu(data.f_temp_3->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_temp_4_gpu(data.f_temp_4->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_temp_5_gpu(data.f_temp_5->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_temp_6_gpu(data.f_temp_6->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_temp_7_gpu(data.f_temp_7->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_temp_8_gpu(data.f_temp_8->lock(lm_write_only, tags::GPU::CUDA::memory_value));

                unsigned long start((*info.limits)[0]);
                unsigned long end((*info.limits)[info.limits->size() - 1]);

                cuda_collide_stream_grid_float(start, end,
                        cuda_dir_1_gpu, cuda_dir_2_gpu, cuda_dir_3_gpu, cuda_dir_4_gpu,
                        cuda_dir_5_gpu, cuda_dir_6_gpu, cuda_dir_7_gpu, cuda_dir_8_gpu,
                        f_eq_0_gpu, f_eq_1_gpu, f_eq_2_gpu,
                        f_eq_3_gpu, f_eq_4_gpu, f_eq_5_gpu,
                        f_eq_6_gpu, f_eq_7_gpu, f_eq_8_gpu,
                        f_0_gpu, f_1_gpu, f_2_gpu,
                        f_3_gpu, f_4_gpu, f_5_gpu,
                        f_6_gpu, f_7_gpu, f_8_gpu,
                        f_temp_0_gpu, f_temp_1_gpu, f_temp_2_gpu,
                        f_temp_3_gpu, f_temp_4_gpu, f_temp_5_gpu,
                        f_temp_6_gpu, f_temp_7_gpu, f_temp_8_gpu,
                        tau, data.h->size(),
                        blocksize);

                info.limits->unlock(lm_read_only);

                info.cuda_dir_1->unlock(lm_read_only);
                info.cuda_dir_2->unlock(lm_read_only);
                info.cuda_dir_3->unlock(lm_read_only);
                info.cuda_dir_4->unlock(lm_read_only);
                info.cuda_dir_5->unlock(lm_read_only);
                info.cuda_dir_6->unlock(lm_read_only);
                info.cuda_dir_7->unlock(lm_read_only);
                info.cuda_dir_8->unlock(lm_read_only);

                data.f_eq_0->unlock(lm_read_only);
                data.f_eq_1->unlock(lm_read_only);
                data.f_eq_2->unlock(lm_read_only);
                data.f_eq_3->unlock(lm_read_only);
                data.f_eq_4->unlock(lm_read_only);
                data.f_eq_5->unlock(lm_read_only);
                data.f_eq_6->unlock(lm_read_only);
                data.f_eq_7->unlock(lm_read_only);
                data.f_eq_8->unlock(lm_read_only);

                data.f_0->unlock(lm_read_only);
                data.f_1->unlock(lm_read_only);
                data.f_2->unlock(lm_read_only);
                data.f_3->unlock(lm_read_only);
                data.f_4->unlock(lm_read_only);
                data.f_5->unlock(lm_read_only);
                data.f_6->unlock(lm_read_only);
                data.f_7->unlock(lm_read_only);
                data.f_8->unlock(lm_read_only);

                data.f_temp_0->unlock(lm_write_only);
                data.f_temp_1->unlock(lm_write_only);
                data.f_temp_2->unlock(lm_write_only);
                data.f_temp_3->unlock(lm_write_only);
                data.f_temp_4->unlock(lm_write_only);
                data.f_temp_5->unlock(lm_write_only);
                data.f_temp_6->unlock(lm_write_only);
                data.f_temp_7->unlock(lm_write_only);
                data.f_temp_8->unlock(lm_write_only);
            }
    };
}

void CollideStreamGrid<tags::GPU::CUDA, lbm_boundary_types::NOSLIP,
     lbm_lattice_types::D2Q9>::value(
             PackedGridInfo<lbm_lattice_types::D2Q9> & info,
             PackedGridData<lbm_lattice_types::D2Q9, float> & data,
             float tau)
{
    CONTEXT("When performing collision and streaming (CUDA):");


    unsigned long blocksize(Configuration::instance()->get_value("cuda::collide_stream_grid_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaColStreamGridfloat task(info, data, tau, blocksize);
        task();
    }
    else
    {
        cudaColStreamGridfloat task(info, data, tau, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }
}

