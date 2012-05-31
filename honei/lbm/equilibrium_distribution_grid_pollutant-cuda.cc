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

#include <honei/lbm/equilibrium_distribution_grid_pollutant.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/backends/cuda/gpu_pool.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/util/configuration.hh>


using namespace honei;

namespace
{
    class cudaEqDistGridPollutantfloat
    {
        private:
            float e;
            PackedGridInfo<D2Q9> & info;
            PackedGridData<D2Q9, float> & data_flow;
            PackedGridData<D2Q9, float> & data_poll;
            unsigned long blocksize;
        public:
            cudaEqDistGridPollutantfloat(float e, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data_flow, PackedGridData<D2Q9, float> & data_poll, unsigned long blocksize) :
                e(e),
                info(info),
                data_flow(data_flow),
                data_poll(data_poll),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                info.limits->lock(lm_read_only);

                void * u_gpu(data_flow.u->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * v_gpu(data_flow.v->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * h_gpu(data_flow.h->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * c_gpu(data_poll.h->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * distribution_x_gpu(data_flow.distribution_x->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * distribution_y_gpu(data_flow.distribution_y->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * f_eq_0_gpu(data_poll.f_eq_0->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_eq_1_gpu(data_poll.f_eq_1->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_eq_2_gpu(data_poll.f_eq_2->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_eq_3_gpu(data_poll.f_eq_3->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_eq_4_gpu(data_poll.f_eq_4->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_eq_5_gpu(data_poll.f_eq_5->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_eq_6_gpu(data_poll.f_eq_6->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_eq_7_gpu(data_poll.f_eq_7->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_eq_8_gpu(data_poll.f_eq_8->lock(lm_write_only, tags::GPU::CUDA::memory_value));

                unsigned long start((*info.limits)[0]);
                unsigned long end((*info.limits)[info.limits->size() - 1]);
                cuda_eq_dist_grid_pollutant_float(start, end, u_gpu, v_gpu, h_gpu, c_gpu,
                        distribution_x_gpu, distribution_y_gpu,
                        f_eq_0_gpu, f_eq_1_gpu, f_eq_2_gpu,
                        f_eq_3_gpu, f_eq_4_gpu, f_eq_5_gpu,
                        f_eq_6_gpu, f_eq_7_gpu, f_eq_8_gpu,
                        e,
                        blocksize);

                info.limits->unlock(lm_read_only);

                data_flow.u->unlock(lm_read_only);
                data_flow.v->unlock(lm_read_only);
                data_flow.h->unlock(lm_read_only);
                data_poll.h->unlock(lm_read_only);

                data_flow.distribution_x->unlock(lm_read_only);
                data_flow.distribution_y->unlock(lm_read_only);

                data_poll.f_eq_0->unlock(lm_write_only);
                data_poll.f_eq_1->unlock(lm_write_only);
                data_poll.f_eq_2->unlock(lm_write_only);
                data_poll.f_eq_3->unlock(lm_write_only);
                data_poll.f_eq_4->unlock(lm_write_only);
                data_poll.f_eq_5->unlock(lm_write_only);
                data_poll.f_eq_6->unlock(lm_write_only);
                data_poll.f_eq_7->unlock(lm_write_only);
                data_poll.f_eq_8->unlock(lm_write_only);
            }
    };

#ifdef HONEI_CUDA_DOUBLE
    class cudaEqDistGridPollutantdouble
    {
        private:
            double e;
            PackedGridInfo<D2Q9> & info;
            PackedGridData<D2Q9, double> & data_flow;
            PackedGridData<D2Q9, double> & data_poll;
            unsigned long blocksize;
        public:
            cudaEqDistGridPollutantdouble(double e, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, double> & data_flow, PackedGridData<D2Q9, double> & data_poll, unsigned long blocksize) :
                e(e),
                info(info),
                data_flow(data_flow),
                data_poll(data_poll),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                info.limits->lock(lm_read_only);

                void * u_gpu(data_flow.u->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * v_gpu(data_flow.v->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * h_gpu(data_flow.h->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * c_gpu(data_poll.h->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * distribution_x_gpu(data_flow.distribution_x->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * distribution_y_gpu(data_flow.distribution_y->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * f_eq_0_gpu(data_poll.f_eq_0->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_eq_1_gpu(data_poll.f_eq_1->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_eq_2_gpu(data_poll.f_eq_2->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_eq_3_gpu(data_poll.f_eq_3->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_eq_4_gpu(data_poll.f_eq_4->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_eq_5_gpu(data_poll.f_eq_5->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_eq_6_gpu(data_poll.f_eq_6->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_eq_7_gpu(data_poll.f_eq_7->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_eq_8_gpu(data_poll.f_eq_8->lock(lm_write_only, tags::GPU::CUDA::memory_value));

                unsigned long start((*info.limits)[0]);
                unsigned long end((*info.limits)[info.limits->size() - 1]);
                cuda_eq_dist_grid_pollutant_double(start, end, u_gpu, v_gpu, h_gpu, c_gpu,
                        distribution_x_gpu, distribution_y_gpu,
                        f_eq_0_gpu, f_eq_1_gpu, f_eq_2_gpu,
                        f_eq_3_gpu, f_eq_4_gpu, f_eq_5_gpu,
                        f_eq_6_gpu, f_eq_7_gpu, f_eq_8_gpu,
                        e,
                        blocksize);

                info.limits->unlock(lm_read_only);

                data_flow.u->unlock(lm_read_only);
                data_flow.v->unlock(lm_read_only);
                data_flow.h->unlock(lm_read_only);
                data_poll.h->unlock(lm_read_only);

                data_flow.distribution_x->unlock(lm_read_only);
                data_flow.distribution_y->unlock(lm_read_only);

                data_poll.f_eq_0->unlock(lm_write_only);
                data_poll.f_eq_1->unlock(lm_write_only);
                data_poll.f_eq_2->unlock(lm_write_only);
                data_poll.f_eq_3->unlock(lm_write_only);
                data_poll.f_eq_4->unlock(lm_write_only);
                data_poll.f_eq_5->unlock(lm_write_only);
                data_poll.f_eq_6->unlock(lm_write_only);
                data_poll.f_eq_7->unlock(lm_write_only);
                data_poll.f_eq_8->unlock(lm_write_only);
            }
    };
#endif
}

void EquilibriumDistributionGridPollutant<tags::GPU::CUDA, lbm_applications::LABSWE>::value(float e,
        PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data_flow, PackedGridData<D2Q9, float> & data_poll)
{
    CONTEXT("When computing LABSWE local equilibrium distribution function (CUDA):");


    unsigned long blocksize(Configuration::instance()->get_value("cuda::eq_dist_grid_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaEqDistGridPollutantfloat task(e, info, data_flow, data_poll, blocksize);
        task();
    }
    else
    {
        cudaEqDistGridPollutantfloat task(e, info, data_flow, data_poll, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }
}

#ifdef HONEI_CUDA_DOUBLE
void EquilibriumDistributionGridPollutant<tags::GPU::CUDA, lbm_applications::LABSWE>::value(double e,
        PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, double> & data_flow, PackedGridData<D2Q9, double> & data_poll)
{
    CONTEXT("When computing LABSWE local equilibrium distribution function (CUDA):");


    unsigned long blocksize(Configuration::instance()->get_value("cuda::eq_dist_grid_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaEqDistGridPollutantdouble task(e, info, data_flow, data_poll, blocksize);
        task();
    }
    else
    {
        cudaEqDistGridPollutantdouble task(e, info, data_flow, data_poll, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }
}
#endif
