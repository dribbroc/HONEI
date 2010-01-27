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

#include <honei/la/element_inverse.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/backends/cuda/gpu_pool.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/util/configuration.hh>


using namespace honei;

namespace
{
    class cudaElementInverseDVfloat
    {
        private:
            DenseVectorContinuousBase<float> & a;
            unsigned long blocksize;
        public:
            cudaElementInverseDVfloat(DenseVectorContinuousBase<float> & a, unsigned long blocksize) :
                a(a),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * a_gpu (a.lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                cuda_element_inverse_one_float(a_gpu, a.size(), blocksize);
                a.unlock(lm_read_and_write);
            }
    };

    class cudaElementInverseDMfloat
    {
        private:
            DenseMatrix<float> & a;
            unsigned long blocksize;
        public:
            cudaElementInverseDMfloat(DenseMatrix<float> & a, unsigned long blocksize) :
                a(a),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * a_gpu (a.lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                cuda_element_inverse_one_float(a_gpu, a.size(), blocksize);
                a.unlock(lm_read_and_write);
            }
    };
}

DenseVectorContinuousBase<float> & ElementInverse<tags::GPU::CUDA>::value(DenseVectorContinuousBase<float> & x)
{
    CONTEXT("When inverting DenseVectorContinuousBase<float> (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::element_inverse_one_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaElementInverseDVfloat task(x, blocksize);
        task();
    }
    else
    {
        cudaElementInverseDVfloat task(x, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0)->wait();
    }

    return x;
}

DenseMatrix<float> & ElementInverse<tags::GPU::CUDA>::value(DenseMatrix<float> & x)
{
    CONTEXT("When inverting DenseMatrix<float> (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::element_inverse_one_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaElementInverseDMfloat task(x, blocksize);
        task();
    }
    else
    {
        cudaElementInverseDMfloat task(x, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0)->wait();
    }

    return x;
}
