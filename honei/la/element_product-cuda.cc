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

#include <honei/la/element_product.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/backends/cuda/gpu_pool.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/util/configuration.hh>


using namespace honei;

namespace
{
    class cudaElementProductDVfloat
    {
        private:
            DenseVectorContinuousBase<float> & a;
            const DenseVectorContinuousBase<float> & b;
            unsigned long blocksize;
        public:
            cudaElementProductDVfloat(DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b, unsigned long blocksize) :
                a(a),
                b(b),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * a_gpu(a.lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * b_gpu(b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                cuda_element_product_three_float(a_gpu, a_gpu, b_gpu, a.size(), blocksize);
                b.unlock(lm_read_only);
                a.unlock(lm_read_and_write);
            }
    };

    class cudaElementProductDVdouble
    {
        private:
            DenseVectorContinuousBase<double> & a;
            const DenseVectorContinuousBase<double> & b;
            unsigned long blocksize;
        public:
            cudaElementProductDVdouble(DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b, unsigned long blocksize) :
                a(a),
                b(b),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * a_gpu(a.lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * b_gpu(b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                cuda_element_product_three_double(a_gpu, a_gpu, b_gpu, a.size(), blocksize);
                b.unlock(lm_read_only);
                a.unlock(lm_read_and_write);
            }
    };

    class cudaElementProduct3DVfloat
    {
        private:
            DenseVectorContinuousBase<float> & result;
            const DenseVectorContinuousBase<float> & a;
            const DenseVectorContinuousBase<float> & b;
            unsigned long blocksize;
        public:
            cudaElementProduct3DVfloat(DenseVectorContinuousBase<float> & result, const DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b, unsigned long blocksize) :
                result(result),
                a(a),
                b(b),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * result_gpu(result.lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * a_gpu(a.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * b_gpu(b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                cuda_element_product_three_float(result_gpu, a_gpu, b_gpu, a.size(), blocksize);
                result.unlock(lm_write_only);
                b.unlock(lm_read_only);
                a.unlock(lm_read_only);
            }
    };

    class cudaElementProduct3DVdouble
    {
        private:
            DenseVectorContinuousBase<double> & result;
            const DenseVectorContinuousBase<double> & a;
            const DenseVectorContinuousBase<double> & b;
            unsigned long blocksize;
        public:
            cudaElementProduct3DVdouble(DenseVectorContinuousBase<double> & result, const DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b, unsigned long blocksize) :
                result(result),
                a(a),
                b(b),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * result_gpu(result.lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * a_gpu(a.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * b_gpu(b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                cuda_element_product_three_double(result_gpu, a_gpu, b_gpu, a.size(), blocksize);
                result.unlock(lm_write_only);
                b.unlock(lm_read_only);
                a.unlock(lm_read_only);
            }
    };

    class cudaElementProductDMfloat
    {
        private:
            DenseMatrix<float> & a;
            const DenseMatrix<float> & b;
            unsigned long blocksize;
        public:
            cudaElementProductDMfloat(DenseMatrix<float> & a, const DenseMatrix<float> & b, unsigned long blocksize) :
                a(a),
                b(b),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * a_gpu(a.lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * b_gpu(b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                cuda_element_product_three_float(a_gpu, a_gpu, b_gpu, a.size(), blocksize);
                b.unlock(lm_read_only);
                a.unlock(lm_read_and_write);
            }
    };
}

DenseVectorContinuousBase<float> & ElementProduct<tags::GPU::CUDA>::value(DenseVectorContinuousBase<float> & a,
        const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When multiplying DenseVectorContinuousBase<float> and DenseVectorContinuousBase<float> elementwise "
            "(CUDA):");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::element_inverse_one_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaElementProductDVfloat task(a, b, blocksize);
        task();
    }
    else
    {
        cudaElementProductDVfloat task(a, b, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0)->wait();
    }

    return a;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVectorContinuousBase<double> & ElementProduct<tags::GPU::CUDA>::value(DenseVectorContinuousBase<double> & a,
        const DenseVectorContinuousBase<double> & b)
{
    CONTEXT("When multiplying DenseVectorContinuousBase<double> and DenseVectorContinuousBase<double> elementwise "
            "(CUDA):");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::element_inverse_one_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaElementProductDVdouble task(a, b, blocksize);
        task();
    }
    else
    {
        cudaElementProductDVdouble task(a, b, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0)->wait();
    }

    return a;
}
#endif

DenseVectorContinuousBase<float> & ElementProduct<tags::GPU::MultiCore::CUDA>::value(DenseVectorContinuousBase<float> & a,
        const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When multiplying DenseVectorContinuousBase<float> and DenseVectorContinuousBase<float> elementwise "
            "(CUDA):");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::element_inverse_one_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        throw InternalError("You should not run this operation within any MC CUDA op!");
    }
    else
    {
        DenseVectorRange<float> a1(a.range(a.size()/2, 0));
        DenseVectorRange<float> b1(b.range(b.size()/2, 0));
        cudaElementProductDVfloat task1(a1, b1, blocksize);
        DenseVectorRange<float> a2(a.range(a.size()/2 + a.size()%2, a.size()/2));
        DenseVectorRange<float> b2(b.range(b.size()/2 + b.size()%2, b.size()/2));
        cudaElementProductDVfloat task2(a2, b2, blocksize);
        cuda::GPUPool::instance()->enqueue(task1, 0)->wait();
        cuda::GPUPool::instance()->enqueue(task2, 1)->wait();
    }

    return a;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVectorContinuousBase<double> & ElementProduct<tags::GPU::MultiCore::CUDA>::value(DenseVectorContinuousBase<double> & a,
        const DenseVectorContinuousBase<double> & b)
{
    CONTEXT("When multiplying DenseVectorContinuousBase<double> and DenseVectorContinuousBase<double> elementwise "
            "(CUDA):");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::element_inverse_one_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        throw InternalError("You should not run this operation within any MC CUDA op!");
    }
    else
    {
        DenseVectorRange<double> a1(a.range(a.size()/2, 0));
        DenseVectorRange<double> b1(b.range(b.size()/2, 0));
        cudaElementProductDVdouble task1(a1, b1, blocksize);
        DenseVectorRange<double> a2(a.range(a.size()/2 + a.size()%2, a.size()/2));
        DenseVectorRange<double> b2(b.range(b.size()/2 + b.size()%2, b.size()/2));
        cudaElementProductDVdouble task2(a2, b2, blocksize);
        cuda::GPUPool::instance()->enqueue(task1, 0)->wait();
        cuda::GPUPool::instance()->enqueue(task2, 1)->wait();
    }

    return a;
}
#endif

DenseVectorContinuousBase<float> & ElementProduct<tags::GPU::CUDA>::value(DenseVectorContinuousBase<float> & result, DenseVectorContinuousBase<float> & a,
        const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When multiplying DenseVectorContinuousBase<float> and DenseVectorContinuousBase<float> elementwise "
            "(CUDA):");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::element_inverse_one_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaElementProduct3DVfloat task(result, a, b, blocksize);
        task();
    }
    else
    {
        cudaElementProduct3DVfloat task(result, a, b, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0)->wait();
    }

    return result;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVectorContinuousBase<double> & ElementProduct<tags::GPU::CUDA>::value(DenseVectorContinuousBase<double> & result, DenseVectorContinuousBase<double> & a,
        const DenseVectorContinuousBase<double> & b)
{
    CONTEXT("When multiplying DenseVectorContinuousBase<double> and DenseVectorContinuousBase<double> elementwise "
            "(CUDA):");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::element_inverse_one_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaElementProduct3DVdouble task(result, a, b, blocksize);
        task();
    }
    else
    {
        cudaElementProduct3DVdouble task(result, a, b, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0)->wait();
    }

    return result;
}
#endif

DenseVectorContinuousBase<float> & ElementProduct<tags::GPU::MultiCore::CUDA>::value(DenseVectorContinuousBase<float> & result, DenseVectorContinuousBase<float> & a,
        const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When multiplying DenseVectorContinuousBase<float> and DenseVectorContinuousBase<float> elementwise "
            "(CUDA):");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::element_inverse_one_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        throw InternalError("You should not run this operation within any MC CUDA op!");
    }
    else
    {
        DenseVectorRange<float> r1(result.range(result.size()/2, 0));
        DenseVectorRange<float> a1(a.range(a.size()/2, 0));
        DenseVectorRange<float> b1(b.range(b.size()/2, 0));
        cudaElementProduct3DVfloat task1(r1, a1, b1, blocksize);
        DenseVectorRange<float> r2(result.range(result.size()/2 + result.size()%2, result.size()/2));
        DenseVectorRange<float> a2(a.range(a.size()/2 + a.size()%2, a.size()/2));
        DenseVectorRange<float> b2(b.range(b.size()/2 + b.size()%2, b.size()/2));
        cudaElementProduct3DVfloat task2(r2, a2, b2, blocksize);
        cuda::GPUPool::instance()->enqueue(task1, 0)->wait();
        cuda::GPUPool::instance()->enqueue(task2, 1)->wait();
    }

    return result;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVectorContinuousBase<double> & ElementProduct<tags::GPU::MultiCore::CUDA>::value(DenseVectorContinuousBase<double> & result, DenseVectorContinuousBase<double> & a,
        const DenseVectorContinuousBase<double> & b)
{
    CONTEXT("When multiplying DenseVectorContinuousBase<double> and DenseVectorContinuousBase<double> elementwise "
            "(CUDA):");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::element_inverse_one_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        throw InternalError("You should not run this operation within any MC CUDA op!");
    }
    else
    {
        DenseVectorRange<double> r1(result.range(result.size()/2, 0));
        DenseVectorRange<double> a1(a.range(a.size()/2, 0));
        DenseVectorRange<double> b1(b.range(b.size()/2, 0));
        cudaElementProduct3DVdouble task1(r1, a1, b1, blocksize);
        DenseVectorRange<double> r2(result.range(result.size()/2 + result.size()%2, result.size()/2));
        DenseVectorRange<double> a2(a.range(a.size()/2 + a.size()%2, a.size()/2));
        DenseVectorRange<double> b2(b.range(b.size()/2 + b.size()%2, b.size()/2));
        cudaElementProduct3DVdouble task2(r1, a2, b2, blocksize);
        cuda::GPUPool::instance()->enqueue(task1, 0)->wait();
        cuda::GPUPool::instance()->enqueue(task2, 1)->wait();
    }

    return result;
}
#endif

DenseMatrix<float> & ElementProduct<tags::GPU::CUDA>::value(DenseMatrix<float> & a, const DenseMatrix<float> & b)
{
    CONTEXT("When multiplying DenseMatrix<float> and DenseMatrix<float> elementwise (CUDA):");

    if (a.columns() != b.columns())
    {
        throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
    }

    if (a.rows() != b.rows())
    {
        throw MatrixRowsDoNotMatch(b.rows(), a.rows());
    }

    unsigned long blocksize(Configuration::instance()->get_value("cuda::element_inverse_one_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaElementProductDMfloat task(a, b, blocksize);
        task();
    }
    else
    {
        cudaElementProductDMfloat task(a, b, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0)->wait();
    }

    return a;
}

