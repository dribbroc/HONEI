/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c)  2008, 2010 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/math/prolongation.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/backends/cuda/gpu_pool.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/util/configuration.hh>


using namespace honei;

namespace
{
    class cudaProlongationDVfloat
    {
        private:
            DenseVector<float> & fine;
            const DenseVector<float> & coarse;
            const DenseVector<unsigned long> & mask;
            unsigned long blocksize;
        public:
            cudaProlongationDVfloat(DenseVector<float> & fine, const DenseVector<float> & coarse, const DenseVector<unsigned long> & mask, unsigned long blocksize) :
                fine(fine),
                coarse(coarse),
                mask(mask),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * fine_gpu (fine.lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * coarse_gpu (coarse.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                cuda_prolongation_float(fine_gpu, fine.size(), coarse_gpu, coarse.size(), mask.elements(), blocksize);
                coarse.unlock(lm_read_only);
                fine.unlock(lm_write_only);
            }
    };

    class cudaProlongationDVdouble
    {
        private:
            DenseVector<double> & fine;
            const DenseVector<double> & coarse;
            const DenseVector<unsigned long> & mask;
            unsigned long blocksize;
        public:
            cudaProlongationDVdouble(DenseVector<double> & fine, const DenseVector<double> & coarse, const DenseVector<unsigned long> & mask, unsigned long blocksize) :
                fine(fine),
                coarse(coarse),
                mask(mask),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * fine_gpu (fine.lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * coarse_gpu (coarse.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                cuda_prolongation_double(fine_gpu, fine.size(), coarse_gpu, coarse.size(), mask.elements(), blocksize);
                coarse.unlock(lm_read_only);
                fine.unlock(lm_write_only);
            }
    };
}

DenseVector<float> & Prolongation<tags::GPU::CUDA, NONE>::value(DenseVector<float> & fine,
        const DenseVector<float> & coarse, const DenseVector<unsigned long> & mask, BandedMatrixQ1<float> & prolmat)
{
    CONTEXT("When prolongating from coarse to fine (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::prolongation_float", 64ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaProlongationDVfloat task(fine, coarse, mask, blocksize);
        task();
    }
    else
    {
        cudaProlongationDVfloat task(fine, coarse, mask, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0)->wait();
    }

    return fine;
}

DenseVector<float> & Prolongation<tags::GPU::CUDA, NONE>::value(DenseVector<float> & fine,
        const DenseVector<float> & coarse, const DenseVector<unsigned long> & mask, SparseMatrixELL<float> & prolmat)
{
    CONTEXT("When prolongating from coarse to fine (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::prolongation_float", 64ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaProlongationDVfloat task(fine, coarse, mask, blocksize);
        task();
    }
    else
    {
        cudaProlongationDVfloat task(fine, coarse, mask, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0)->wait();
    }

    return fine;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVector<double> & Prolongation<tags::GPU::CUDA, NONE>::value(DenseVector<double> & fine,
        const DenseVector<double> & coarse, const DenseVector<unsigned long> & mask, BandedMatrixQ1<double> & prolmat)
{
    CONTEXT("When prolongating from coarse to fine (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::prolongation_double", 64ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaProlongationDVdouble task(fine, coarse, mask, blocksize);
        task();
    }
    else
    {
        cudaProlongationDVdouble task(fine, coarse, mask, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0)->wait();
    }

    return fine;
}
DenseVector<double> & Prolongation<tags::GPU::CUDA, NONE>::value(DenseVector<double> & fine,
        const DenseVector<double> & coarse, const DenseVector<unsigned long> & mask, SparseMatrixELL<double> & prolmat)
{
    CONTEXT("When prolongating from coarse to fine (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::prolongation_double", 64ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaProlongationDVdouble task(fine, coarse, mask, blocksize);
        task();
    }
    else
    {
        cudaProlongationDVdouble task(fine, coarse, mask, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0)->wait();
    }

    return fine;
}
#endif
