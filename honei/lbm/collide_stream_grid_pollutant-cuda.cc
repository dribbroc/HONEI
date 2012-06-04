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

#include <honei/lbm/collide_stream_grid_pollutant.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/backends/cuda/gpu_pool.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/util/configuration.hh>


using namespace honei;

namespace
{
    class cudaColStreamGridPollutantfloat
    {
        private:
            PackedGridInfo<D2Q9> & info;
            PackedGridData<D2Q9, float> & data_flow;
            PackedGridData<D2Q9, float> & data_poll;
            float tau;
            unsigned long blocksize;
        public:
            cudaColStreamGridPollutantfloat(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data_flow, PackedGridData<D2Q9, float> & data_poll, float tau, unsigned long blocksize) :
                info(info),
                data_flow(data_flow),
                data_poll(data_poll),
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

                void * f_eq_0_gpu(data_poll.f_eq_0->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_eq_1_gpu(data_poll.f_eq_1->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_eq_2_gpu(data_poll.f_eq_2->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_eq_3_gpu(data_poll.f_eq_3->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_eq_4_gpu(data_poll.f_eq_4->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_eq_5_gpu(data_poll.f_eq_5->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_eq_6_gpu(data_poll.f_eq_6->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_eq_7_gpu(data_poll.f_eq_7->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_eq_8_gpu(data_poll.f_eq_8->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * f_0_gpu(data_poll.f_0->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_1_gpu(data_poll.f_1->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_2_gpu(data_poll.f_2->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_3_gpu(data_poll.f_3->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_4_gpu(data_poll.f_4->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_5_gpu(data_poll.f_5->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_6_gpu(data_poll.f_6->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_7_gpu(data_poll.f_7->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_8_gpu(data_poll.f_8->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * f_temp_0_gpu(data_poll.f_temp_0->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_temp_1_gpu(data_poll.f_temp_1->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_temp_2_gpu(data_poll.f_temp_2->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_temp_3_gpu(data_poll.f_temp_3->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_temp_4_gpu(data_poll.f_temp_4->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_temp_5_gpu(data_poll.f_temp_5->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_temp_6_gpu(data_poll.f_temp_6->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_temp_7_gpu(data_poll.f_temp_7->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_temp_8_gpu(data_poll.f_temp_8->lock(lm_write_only, tags::GPU::CUDA::memory_value));

                void * h_gpu(data_flow.h->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                unsigned long start((*info.limits)[0]);
                unsigned long end((*info.limits)[info.limits->size() - 1]);

                cuda_collide_stream_grid_pollutant_float(start, end, h_gpu,
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
                        tau, data_flow.h->size(),
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

                data_poll.f_eq_0->unlock(lm_read_only);
                data_poll.f_eq_1->unlock(lm_read_only);
                data_poll.f_eq_2->unlock(lm_read_only);
                data_poll.f_eq_3->unlock(lm_read_only);
                data_poll.f_eq_4->unlock(lm_read_only);
                data_poll.f_eq_5->unlock(lm_read_only);
                data_poll.f_eq_6->unlock(lm_read_only);
                data_poll.f_eq_7->unlock(lm_read_only);
                data_poll.f_eq_8->unlock(lm_read_only);

                data_poll.f_0->unlock(lm_read_only);
                data_poll.f_1->unlock(lm_read_only);
                data_poll.f_2->unlock(lm_read_only);
                data_poll.f_3->unlock(lm_read_only);
                data_poll.f_4->unlock(lm_read_only);
                data_poll.f_5->unlock(lm_read_only);
                data_poll.f_6->unlock(lm_read_only);
                data_poll.f_7->unlock(lm_read_only);
                data_poll.f_8->unlock(lm_read_only);

                data_poll.f_temp_0->unlock(lm_write_only);
                data_poll.f_temp_1->unlock(lm_write_only);
                data_poll.f_temp_2->unlock(lm_write_only);
                data_poll.f_temp_3->unlock(lm_write_only);
                data_poll.f_temp_4->unlock(lm_write_only);
                data_poll.f_temp_5->unlock(lm_write_only);
                data_poll.f_temp_6->unlock(lm_write_only);
                data_poll.f_temp_7->unlock(lm_write_only);
                data_poll.f_temp_8->unlock(lm_write_only);

                data_flow.h->unlock(lm_read_only);
            }
    };

#ifdef HONEI_CUDA_DOUBLE
    class cudaColStreamGridPollutantdouble
    {
        private:
            PackedGridInfo<D2Q9> & info;
            PackedGridData<D2Q9, double> & data_flow;
            PackedGridData<D2Q9, double> & data_poll;
            double tau;
            unsigned long blocksize;
        public:
            cudaColStreamGridPollutantdouble(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, double> & data_flow, PackedGridData<D2Q9, double> & data_poll, double tau, unsigned long blocksize) :
                info(info),
                data_flow(data_flow),
                data_poll(data_poll),
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

                void * f_eq_0_gpu(data_poll.f_eq_0->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_eq_1_gpu(data_poll.f_eq_1->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_eq_2_gpu(data_poll.f_eq_2->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_eq_3_gpu(data_poll.f_eq_3->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_eq_4_gpu(data_poll.f_eq_4->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_eq_5_gpu(data_poll.f_eq_5->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_eq_6_gpu(data_poll.f_eq_6->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_eq_7_gpu(data_poll.f_eq_7->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_eq_8_gpu(data_poll.f_eq_8->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * f_0_gpu(data_poll.f_0->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_1_gpu(data_poll.f_1->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_2_gpu(data_poll.f_2->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_3_gpu(data_poll.f_3->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_4_gpu(data_poll.f_4->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_5_gpu(data_poll.f_5->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_6_gpu(data_poll.f_6->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_7_gpu(data_poll.f_7->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * f_8_gpu(data_poll.f_8->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * f_temp_0_gpu(data_poll.f_temp_0->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_temp_1_gpu(data_poll.f_temp_1->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_temp_2_gpu(data_poll.f_temp_2->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_temp_3_gpu(data_poll.f_temp_3->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_temp_4_gpu(data_poll.f_temp_4->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_temp_5_gpu(data_poll.f_temp_5->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_temp_6_gpu(data_poll.f_temp_6->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_temp_7_gpu(data_poll.f_temp_7->lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * f_temp_8_gpu(data_poll.f_temp_8->lock(lm_write_only, tags::GPU::CUDA::memory_value));

                void * h_gpu(data_flow.h->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                unsigned long start((*info.limits)[0]);
                unsigned long end((*info.limits)[info.limits->size() - 1]);

                cuda_collide_stream_grid_pollutant_double(start, end, h_gpu,
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
                        tau, data_flow.h->size(),
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

                data_poll.f_eq_0->unlock(lm_read_only);
                data_poll.f_eq_1->unlock(lm_read_only);
                data_poll.f_eq_2->unlock(lm_read_only);
                data_poll.f_eq_3->unlock(lm_read_only);
                data_poll.f_eq_4->unlock(lm_read_only);
                data_poll.f_eq_5->unlock(lm_read_only);
                data_poll.f_eq_6->unlock(lm_read_only);
                data_poll.f_eq_7->unlock(lm_read_only);
                data_poll.f_eq_8->unlock(lm_read_only);

                data_poll.f_0->unlock(lm_read_only);
                data_poll.f_1->unlock(lm_read_only);
                data_poll.f_2->unlock(lm_read_only);
                data_poll.f_3->unlock(lm_read_only);
                data_poll.f_4->unlock(lm_read_only);
                data_poll.f_5->unlock(lm_read_only);
                data_poll.f_6->unlock(lm_read_only);
                data_poll.f_7->unlock(lm_read_only);
                data_poll.f_8->unlock(lm_read_only);

                data_poll.f_temp_0->unlock(lm_write_only);
                data_poll.f_temp_1->unlock(lm_write_only);
                data_poll.f_temp_2->unlock(lm_write_only);
                data_poll.f_temp_3->unlock(lm_write_only);
                data_poll.f_temp_4->unlock(lm_write_only);
                data_poll.f_temp_5->unlock(lm_write_only);
                data_poll.f_temp_6->unlock(lm_write_only);
                data_poll.f_temp_7->unlock(lm_write_only);
                data_poll.f_temp_8->unlock(lm_write_only);

                data_flow.h->unlock(lm_read_only);
            }
    };
#endif
}

void CollideStreamGridPollutant<tags::GPU::CUDA, lbm_boundary_types::NOSLIP,
     lbm_lattice_types::D2Q9>::value(
             PackedGridInfo<lbm_lattice_types::D2Q9> & info,
             PackedGridData<lbm_lattice_types::D2Q9, float> & data_flow,
             PackedGridData<lbm_lattice_types::D2Q9, float> & data_poll,
             float tau)
{
    CONTEXT("When performing collision and streaming (CUDA):");


    unsigned long blocksize(Configuration::instance()->get_value("cuda::collide_stream_grid_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaColStreamGridPollutantfloat task(info, data_flow, data_poll, tau, blocksize);
        task();
    }
    else
    {
        cudaColStreamGridPollutantfloat task(info, data_flow, data_poll, tau, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }
}

#ifdef HONEI_CUDA_DOUBLE
void CollideStreamGridPollutant<tags::GPU::CUDA, lbm_boundary_types::NOSLIP,
     lbm_lattice_types::D2Q9>::value(
             PackedGridInfo<lbm_lattice_types::D2Q9> & info,
             PackedGridData<lbm_lattice_types::D2Q9, double> & data_flow,
             PackedGridData<lbm_lattice_types::D2Q9, double> & data_poll,
             double tau)
{
    CONTEXT("When performing collision and streaming (CUDA):");


    unsigned long blocksize(Configuration::instance()->get_value("cuda::collide_stream_grid_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaColStreamGridPollutantdouble task(info, data_flow, data_poll, tau, blocksize);
        task();
    }
    else
    {
        cudaColStreamGridPollutantdouble task(info, data_flow, data_poll, tau, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }
}
#endif
