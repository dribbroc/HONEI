/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c)  2012 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/lbm/force_grid_pollutant.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/backends/cuda/gpu_pool.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/util/configuration.hh>


using namespace honei;

namespace
{
    class cudaForceGridPollutantfloat
    {
        private:
            PackedGridInfo<D2Q9> & info;
            PackedGridData<D2Q9, float> & data_flow;
            PackedGridData<D2Q9, float> & data_poll;
            float d_t;
            float k;
            float s_0;
            unsigned long blocksize;
        public:
            cudaForceGridPollutantfloat(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data_flow, PackedGridData<D2Q9, float> & data_poll, float d_t, float k, float s_0, unsigned long blocksize) :
                info(info),
                data_flow(data_flow),
                data_poll(data_poll),
                d_t(d_t),
                k(k),
                s_0(s_0),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * h_gpu(data_flow.h->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * c_gpu(data_poll.h->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * f_temp_1_gpu(data_poll.f_temp_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_2_gpu(data_poll.f_temp_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_3_gpu(data_poll.f_temp_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_4_gpu(data_poll.f_temp_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_5_gpu(data_poll.f_temp_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_6_gpu(data_poll.f_temp_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_7_gpu(data_poll.f_temp_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_8_gpu(data_poll.f_temp_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

                unsigned long start((*info.limits)[0]);
                unsigned long end((*info.limits)[info.limits->size() - 1]);
                cuda_force_grid_pollutant_float(start, end,
                        h_gpu, c_gpu,
                        f_temp_1_gpu, f_temp_2_gpu,
                        f_temp_3_gpu, f_temp_4_gpu, f_temp_5_gpu,
                        f_temp_6_gpu, f_temp_7_gpu, f_temp_8_gpu,
                        d_t, k, s_0,
                        blocksize);


                data_flow.h->unlock(lm_read_only);
                data_poll.h->unlock(lm_read_only);

                data_poll.f_temp_1->unlock(lm_read_and_write);
                data_poll.f_temp_2->unlock(lm_read_and_write);
                data_poll.f_temp_3->unlock(lm_read_and_write);
                data_poll.f_temp_4->unlock(lm_read_and_write);
                data_poll.f_temp_5->unlock(lm_read_and_write);
                data_poll.f_temp_6->unlock(lm_read_and_write);
                data_poll.f_temp_7->unlock(lm_read_and_write);
                data_poll.f_temp_8->unlock(lm_read_and_write);
            }
    };

#ifdef HONEI_CUDA_DOUBLE
    class cudaForceGridPollutantdouble
    {
        private:
            PackedGridInfo<D2Q9> & info;
            PackedGridData<D2Q9, double> & data_flow;
            PackedGridData<D2Q9, double> & data_poll;
            double d_t;
            double k;
            double s_0;
            unsigned long blocksize;
        public:
            cudaForceGridPollutantdouble(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, double> & data_flow, PackedGridData<D2Q9, double> & data_poll, double d_t, double k, double s_0, unsigned long blocksize) :
                info(info),
                data_flow(data_flow),
                data_poll(data_poll),
                d_t(d_t),
                k(k),
                s_0(s_0),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * h_gpu(data_flow.h->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * c_gpu(data_poll.h->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * f_temp_1_gpu(data_poll.f_temp_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_2_gpu(data_poll.f_temp_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_3_gpu(data_poll.f_temp_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_4_gpu(data_poll.f_temp_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_5_gpu(data_poll.f_temp_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_6_gpu(data_poll.f_temp_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_7_gpu(data_poll.f_temp_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_8_gpu(data_poll.f_temp_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

                unsigned long start((*info.limits)[0]);
                unsigned long end((*info.limits)[info.limits->size() - 1]);
                cuda_force_grid_pollutant_double(start, end,
                        h_gpu, c_gpu,
                        f_temp_1_gpu, f_temp_2_gpu,
                        f_temp_3_gpu, f_temp_4_gpu, f_temp_5_gpu,
                        f_temp_6_gpu, f_temp_7_gpu, f_temp_8_gpu,
                        d_t, k, s_0,
                        blocksize);


                data_flow.h->unlock(lm_read_only);
                data_poll.h->unlock(lm_read_only);

                data_poll.f_temp_1->unlock(lm_read_and_write);
                data_poll.f_temp_2->unlock(lm_read_and_write);
                data_poll.f_temp_3->unlock(lm_read_and_write);
                data_poll.f_temp_4->unlock(lm_read_and_write);
                data_poll.f_temp_5->unlock(lm_read_and_write);
                data_poll.f_temp_6->unlock(lm_read_and_write);
                data_poll.f_temp_7->unlock(lm_read_and_write);
                data_poll.f_temp_8->unlock(lm_read_and_write);
            }
    };
#endif

}

void ForceGridPollutant<tags::GPU::CUDA, lbm_applications::LABSWE>::value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data_flow, PackedGridData<D2Q9, float> & data_poll, float dt, float k, float s_0)
{
    CONTEXT("When computing LABSWE source term (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::force_grid_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaForceGridPollutantfloat task(info, data_flow, data_poll, dt, k, s_0, blocksize);
        task();
    }
    else
    {
        cudaForceGridPollutantfloat task(info, data_flow, data_poll, dt, k, s_0, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }
}

#ifdef HONEI_CUDA_DOUBLE
void ForceGridPollutant<tags::GPU::CUDA, lbm_applications::LABSWE>::value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, double> & data_flow, PackedGridData<D2Q9, double> & data_poll, double dt, double k, double s_0)
{
    CONTEXT("When computing LABSWE source term (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::force_grid_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaForceGridPollutantdouble task(info, data_flow, data_poll, dt, k, s_0, blocksize);
        task();
    }
    else
    {
        cudaForceGridPollutantdouble task(info, data_flow, data_poll, dt, k, s_0, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }
}
#endif
