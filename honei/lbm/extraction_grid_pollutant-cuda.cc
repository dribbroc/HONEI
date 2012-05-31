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

#include <honei/lbm/extraction_grid_pollutant.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/backends/cuda/gpu_pool.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/util/configuration.hh>


using namespace honei;

namespace
{
    class cudaExtractionGridPollutantDryfloat
    {
        private:
            PackedGridInfo<D2Q9> & info;
            PackedGridData<D2Q9, float> & data_flow;
            PackedGridData<D2Q9, float> & data_poll;
            float epsilon;
            unsigned long blocksize;
        public:
            cudaExtractionGridPollutantDryfloat(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data_flow, PackedGridData<D2Q9, float> & data_poll ,float epsilon, unsigned long blocksize) :
                info(info),
                data_flow(data_flow),
                data_poll(data_poll),
                epsilon(epsilon),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                info.limits->lock(lm_read_only);

                void * f_0_gpu(data_poll.f_0->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_1_gpu(data_poll.f_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_2_gpu(data_poll.f_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_3_gpu(data_poll.f_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_4_gpu(data_poll.f_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_5_gpu(data_poll.f_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_6_gpu(data_poll.f_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_7_gpu(data_poll.f_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_8_gpu(data_poll.f_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

                void * h_flow(data_flow.h->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * h_poll(data_poll.h->lock(lm_write_only, tags::GPU::CUDA::memory_value));

                unsigned long start((*info.limits)[0]);
                unsigned long end((*info.limits)[info.limits->size() - 1]);

                cuda_extraction_grid_pollutant_dry_float(start, end,
                        f_0_gpu, f_1_gpu, f_2_gpu,
                        f_3_gpu, f_4_gpu, f_5_gpu,
                        f_6_gpu, f_7_gpu, f_8_gpu,
                        h_flow, h_poll, epsilon,
                        blocksize);

                info.limits->unlock(lm_read_only);

                data_poll.f_0->unlock(lm_read_and_write);
                data_poll.f_1->unlock(lm_read_and_write);
                data_poll.f_2->unlock(lm_read_and_write);
                data_poll.f_3->unlock(lm_read_and_write);
                data_poll.f_4->unlock(lm_read_and_write);
                data_poll.f_5->unlock(lm_read_and_write);
                data_poll.f_6->unlock(lm_read_and_write);
                data_poll.f_7->unlock(lm_read_and_write);
                data_poll.f_8->unlock(lm_read_and_write);

                data_flow.h->unlock(lm_read_only);
                data_poll.h->unlock(lm_write_only);
            }
    };

#ifdef HONEI_CUDA_DOUBLE
    class cudaExtractionGridPollutantDrydouble
    {
        private:
            PackedGridInfo<D2Q9> & info;
            PackedGridData<D2Q9, double> & data_flow;
            PackedGridData<D2Q9, double> & data_poll;
            double epsilon;
            unsigned long blocksize;
        public:
            cudaExtractionGridPollutantDrydouble(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, double> & data_flow, PackedGridData<D2Q9, double> & data_poll ,double epsilon, unsigned long blocksize) :
                info(info),
                data_flow(data_flow),
                data_poll(data_poll),
                epsilon(epsilon),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                info.limits->lock(lm_read_only);

                void * f_0_gpu(data_poll.f_0->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_1_gpu(data_poll.f_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_2_gpu(data_poll.f_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_3_gpu(data_poll.f_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_4_gpu(data_poll.f_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_5_gpu(data_poll.f_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_6_gpu(data_poll.f_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_7_gpu(data_poll.f_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_8_gpu(data_poll.f_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

                void * h_flow(data_flow.h->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * h_poll(data_poll.h->lock(lm_write_only, tags::GPU::CUDA::memory_value));

                unsigned long start((*info.limits)[0]);
                unsigned long end((*info.limits)[info.limits->size() - 1]);

                cuda_extraction_grid_pollutant_dry_double(start, end,
                        f_0_gpu, f_1_gpu, f_2_gpu,
                        f_3_gpu, f_4_gpu, f_5_gpu,
                        f_6_gpu, f_7_gpu, f_8_gpu,
                        h_flow, h_poll, epsilon,
                        blocksize);

                info.limits->unlock(lm_read_only);

                data_poll.f_0->unlock(lm_read_and_write);
                data_poll.f_1->unlock(lm_read_and_write);
                data_poll.f_2->unlock(lm_read_and_write);
                data_poll.f_3->unlock(lm_read_and_write);
                data_poll.f_4->unlock(lm_read_and_write);
                data_poll.f_5->unlock(lm_read_and_write);
                data_poll.f_6->unlock(lm_read_and_write);
                data_poll.f_7->unlock(lm_read_and_write);
                data_poll.f_8->unlock(lm_read_and_write);

                data_flow.h->unlock(lm_read_only);
                data_poll.h->unlock(lm_write_only);
            }
    };
#endif

}

void ExtractionGridPollutant<tags::GPU::CUDA, lbm_modes::DRY>::value(
        PackedGridInfo<lbm_lattice_types::D2Q9> & info,
        PackedGridData<lbm_lattice_types::D2Q9, float> & data_flow,
        PackedGridData<lbm_lattice_types::D2Q9, float> & data_poll, float epsilon)
{
    CONTEXT("When extracting h, u and v (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::extraction_grid_float", 128ul));

    //set f to t_temp
    DenseVector<float> * swap;
    swap = data_poll.f_0;
    data_poll.f_0 = data_poll.f_temp_0;
    data_poll.f_temp_0 = swap;
    swap = data_poll.f_1;
    data_poll.f_1 = data_poll.f_temp_1;
    data_poll.f_temp_1 = swap;
    swap = data_poll.f_2;
    data_poll.f_2 = data_poll.f_temp_2;
    data_poll.f_temp_2 = swap;
    swap = data_poll.f_3;
    data_poll.f_3 = data_poll.f_temp_3;
    data_poll.f_temp_3 = swap;
    swap = data_poll.f_4;
    data_poll.f_4 = data_poll.f_temp_4;
    data_poll.f_temp_4 = swap;
    swap = data_poll.f_5;
    data_poll.f_5 = data_poll.f_temp_5;
    data_poll.f_temp_5 = swap;
    swap = data_poll.f_6;
    data_poll.f_6 = data_poll.f_temp_6;
    data_poll.f_temp_6 = swap;
    swap = data_poll.f_7;
    data_poll.f_7 = data_poll.f_temp_7;
    data_poll.f_temp_7 = swap;
    swap = data_poll.f_8;
    data_poll.f_8 = data_poll.f_temp_8;
    data_poll.f_temp_8 = swap;

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaExtractionGridPollutantDryfloat task(info, data_flow, data_poll, epsilon, blocksize);
        task();
    }
    else
    {
        cudaExtractionGridPollutantDryfloat task(info, data_flow, data_poll, epsilon, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }
}

#ifdef HONEI_CUDA_DOUBLE
void ExtractionGridPollutant<tags::GPU::CUDA, lbm_modes::DRY>::value(
        PackedGridInfo<lbm_lattice_types::D2Q9> & info,
        PackedGridData<lbm_lattice_types::D2Q9, double> & data_flow,
        PackedGridData<lbm_lattice_types::D2Q9, double> & data_poll, double epsilon)
{
    CONTEXT("When extracting h, u and v (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::extraction_grid_double", 128ul));

    //set f to t_temp
    DenseVector<double> * swap;
    swap = data_poll.f_0;
    data_poll.f_0 = data_poll.f_temp_0;
    data_poll.f_temp_0 = swap;
    swap = data_poll.f_1;
    data_poll.f_1 = data_poll.f_temp_1;
    data_poll.f_temp_1 = swap;
    swap = data_poll.f_2;
    data_poll.f_2 = data_poll.f_temp_2;
    data_poll.f_temp_2 = swap;
    swap = data_poll.f_3;
    data_poll.f_3 = data_poll.f_temp_3;
    data_poll.f_temp_3 = swap;
    swap = data_poll.f_4;
    data_poll.f_4 = data_poll.f_temp_4;
    data_poll.f_temp_4 = swap;
    swap = data_poll.f_5;
    data_poll.f_5 = data_poll.f_temp_5;
    data_poll.f_temp_5 = swap;
    swap = data_poll.f_6;
    data_poll.f_6 = data_poll.f_temp_6;
    data_poll.f_temp_6 = swap;
    swap = data_poll.f_7;
    data_poll.f_7 = data_poll.f_temp_7;
    data_poll.f_temp_7 = swap;
    swap = data_poll.f_8;
    data_poll.f_8 = data_poll.f_temp_8;
    data_poll.f_temp_8 = swap;

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaExtractionGridPollutantDrydouble task(info, data_flow, data_poll, epsilon, blocksize);
        task();
    }
    else
    {
        cudaExtractionGridPollutantDrydouble task(info, data_flow, data_poll, epsilon, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }
}
#endif
