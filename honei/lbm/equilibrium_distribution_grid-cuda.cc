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

#include <honei/lbm/equilibrium_distribution_grid.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/backends/cuda/gpu_pool.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/util/configuration.hh>


using namespace honei;

namespace
{
    class cudaEqDistGridfloat
    {
        private:
            float g;
            float e;
            PackedGridInfo<D2Q9> & info;
            PackedGridData<D2Q9, float> & data;
            unsigned long blocksize;
        public:
            cudaEqDistGridfloat(float g, float e, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data, unsigned long blocksize) :
                g(g),
                e(e),
                info(info),
                data(data),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                info.limits->lock(lm_read_only);

                void * u_gpu(data.u->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * v_gpu(data.v->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * h_gpu(data.h->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * distribution_x_gpu(data.distribution_x->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * distribution_y_gpu(data.distribution_y->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * f_eq_0_gpu(data.f_eq_0->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_eq_1_gpu(data.f_eq_1->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_eq_2_gpu(data.f_eq_2->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_eq_3_gpu(data.f_eq_3->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_eq_4_gpu(data.f_eq_4->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_eq_5_gpu(data.f_eq_5->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_eq_6_gpu(data.f_eq_6->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_eq_7_gpu(data.f_eq_7->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_eq_8_gpu(data.f_eq_8->lock(lm_write_only, tags::GPU::CUDA::memory_value));

                unsigned long start((*info.limits)[0]);
                unsigned long end((*info.limits)[info.limits->size() - 1]);
                cuda_eq_dist_grid_float(start, end, u_gpu, v_gpu, h_gpu,
                        distribution_x_gpu, distribution_y_gpu,
                        f_eq_0_gpu, f_eq_1_gpu, f_eq_2_gpu,
                        f_eq_3_gpu, f_eq_4_gpu, f_eq_5_gpu,
                        f_eq_6_gpu, f_eq_7_gpu, f_eq_8_gpu,
                        g, e,
                        blocksize);

                info.limits->unlock(lm_read_only);

                data.u->unlock(lm_read_only);
                data.v->unlock(lm_read_only);
                data.h->unlock(lm_read_only);

                data.distribution_x->unlock(lm_read_only);
                data.distribution_y->unlock(lm_read_only);

                data.f_eq_0->unlock(lm_write_only);
                data.f_eq_1->unlock(lm_write_only);
                data.f_eq_2->unlock(lm_write_only);
                data.f_eq_3->unlock(lm_write_only);
                data.f_eq_4->unlock(lm_write_only);
                data.f_eq_5->unlock(lm_write_only);
                data.f_eq_6->unlock(lm_write_only);
                data.f_eq_7->unlock(lm_write_only);
                data.f_eq_8->unlock(lm_write_only);
            }
    };
}

void EquilibriumDistributionGrid<tags::GPU::CUDA, lbm_applications::LABSWE>::value(float g, float e,
        PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data)
{
    CONTEXT("When computing LABSWE local equilibrium distribution function (CUDA):");


    unsigned long blocksize(Configuration::instance()->get_value("cuda::eq_dist_grid_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaEqDistGridfloat task(g, e, info, data, blocksize);
        task();
    }
    else
    {
        cudaEqDistGridfloat task(g, e, info, data, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }
}

