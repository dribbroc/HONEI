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

#include <honei/la/difference.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/util/configuration.hh>
#include <honei/backends/cuda/gpu_pool.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/util/profiler.hh>


using namespace honei;

namespace
{
    class cudaDifferenceDVfloat
    {
        private:
            DenseVectorContinuousBase<float> & a;
            const DenseVectorContinuousBase<float> & b;
            unsigned long blocksize;
        public:
            cudaDifferenceDVfloat(DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b, unsigned long blocksize) :
                a(a),
                b(b),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * a_gpu (a.lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * b_gpu (b.lock(lm_read_only, tags::GPU::CUDA::memory_value));

                cuda_difference_two_float(a_gpu, b_gpu, a.size(), blocksize);

                b.unlock(lm_read_only);
                a.unlock(lm_read_and_write);
            }
    };

    class cudaDifferenceDVdouble
    {
        private:
            DenseVectorContinuousBase<double> & a;
            const DenseVectorContinuousBase<double> & b;
            unsigned long blocksize;
        public:
            cudaDifferenceDVdouble(DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b, unsigned long blocksize) :
                a(a),
                b(b),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * a_gpu (a.lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * b_gpu (b.lock(lm_read_only, tags::GPU::CUDA::memory_value));

                cuda_difference_two_double(a_gpu, b_gpu, a.size(), blocksize);

                b.unlock(lm_read_only);
                a.unlock(lm_read_and_write);
            }
    };

    class cudaDifference3DVfloat
    {
        private:
            DenseVectorContinuousBase<float> & r;
            const DenseVectorContinuousBase<float> & a;
            const DenseVectorContinuousBase<float> & b;
            unsigned long blocksize;
        public:
            cudaDifference3DVfloat(DenseVectorContinuousBase<float> & r, const DenseVectorContinuousBase<float> & a,
                    const DenseVectorContinuousBase<float> & b, unsigned long blocksize) :
                r(r),
                a(a),
                b(b),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * r_gpu (r.lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * a_gpu (a.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * b_gpu (b.lock(lm_read_only, tags::GPU::CUDA::memory_value));

                cuda_difference_three_float(r_gpu, a_gpu, b_gpu, a.size(), blocksize);

                r.unlock(lm_write_only);
                b.unlock(lm_read_only);
                a.unlock(lm_read_only);
            }
    };

#ifdef HONEI_CUDA_DOUBLE
    class cudaDifference3DVdouble
    {
        private:
            DenseVectorContinuousBase<double> & r;
            const DenseVectorContinuousBase<double> & a;
            const DenseVectorContinuousBase<double> & b;
            unsigned long blocksize;
        public:
            cudaDifference3DVdouble(DenseVectorContinuousBase<double> & r, const DenseVectorContinuousBase<double> & a,
                    const DenseVectorContinuousBase<double> & b, unsigned long blocksize) :
                r(r),
                a(a),
                b(b),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * r_gpu (r.lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * a_gpu (a.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * b_gpu (b.lock(lm_read_only, tags::GPU::CUDA::memory_value));

                cuda_difference_three_double(r_gpu, a_gpu, b_gpu, a.size(), blocksize);

                r.unlock(lm_write_only);
                b.unlock(lm_read_only);
                a.unlock(lm_read_only);
            }
    };
#endif

    class cudaDifferenceDMfloat
    {
        private:
            DenseMatrix<float> & a;
            const DenseMatrix<float> & b;
            unsigned long blocksize;
        public:
            cudaDifferenceDMfloat(DenseMatrix<float> & a, const DenseMatrix<float> & b, unsigned long blocksize) :
                a(a),
                b(b),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * a_gpu (a.lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * b_gpu (b.lock(lm_read_only, tags::GPU::CUDA::memory_value));

                cuda_difference_two_float(a_gpu, b_gpu, a.size(), blocksize);

                b.unlock(lm_read_only);
                a.unlock(lm_read_and_write);
            }
    };
}


DenseVectorContinuousBase<float> & Difference<tags::GPU::CUDA>::value(DenseVectorContinuousBase<float> & a,
        const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When subtracting DenseVectorContinuousBase<float> from DenseVectorContinuousBase<float> (CUDA):");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::difference_two_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaDifferenceDVfloat task(a, b, blocksize);
        task();
    }
    else
    {
        cudaDifferenceDVfloat task(a, b, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0)->wait();
    }

    return a;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVectorContinuousBase<double> & Difference<tags::GPU::CUDA>::value(DenseVectorContinuousBase<double> & a,
        const DenseVectorContinuousBase<double> & b)
{
    CONTEXT("When subtracting DenseVectorContinuousBase<double> from DenseVectorContinuousBase<double> (CUDA):");
    PROFILER_START("Difference DV double tags::GPU::CUDA");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::difference_two_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaDifferenceDVdouble task(a, b, blocksize);
        task();
    }
    else
    {
        cudaDifferenceDVdouble task(a, b, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0)->wait();
    }

    PROFILER_STOP("Difference DV double tags::GPU::CUDA");
    return a;
}
#endif

DenseVectorContinuousBase<float> & Difference<tags::GPU::MultiCore::CUDA>::value(DenseVectorContinuousBase<float> & a,
        const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When subtracting DenseVectorContinuousBase<float> from DenseVectorContinuousBase<float> (MC CUDA):");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::difference_two_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        throw InternalError("You should not run this operation within any MC CUDA op!");
    }
    else
    {
        DenseVectorRange<float> a1(a.range(a.size()/2, 0));
        DenseVectorRange<float> b1(b.range(b.size()/2, 0));
        cudaDifferenceDVfloat task1(a1, b1, blocksize);
        DenseVectorRange<float> a2(a.range(a.size()/2 + a.size()%2, a.size()/2));
        DenseVectorRange<float> b2(b.range(b.size()/2 + b.size()%2, b.size()/2));
        cudaDifferenceDVfloat task2(a2, b2, blocksize);
        cuda::GPUPool::instance()->enqueue(task1, 0)->wait();
        cuda::GPUPool::instance()->enqueue(task2, 1)->wait();
    }

    return a;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVectorContinuousBase<double> & Difference<tags::GPU::MultiCore::CUDA>::value(DenseVectorContinuousBase<double> & a,
        const DenseVectorContinuousBase<double> & b)
{
    CONTEXT("When subtracting DenseVectorContinuousBase<double> from DenseVectorContinuousBase<double> (MC CUDA):");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::difference_two_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        throw InternalError("You should not run this operation within any MC CUDA op!");
    }
    else
    {
        DenseVectorRange<double> a1(a.range(a.size()/2, 0));
        DenseVectorRange<double> b1(b.range(b.size()/2, 0));
        cudaDifferenceDVdouble task1(a1, b1, blocksize);
        DenseVectorRange<double> a2(a.range(a.size()/2 + a.size()%2, a.size()/2));
        DenseVectorRange<double> b2(b.range(b.size()/2 + b.size()%2, b.size()/2));
        cudaDifferenceDVdouble task2(a2, b2, blocksize);
        cuda::GPUPool::instance()->enqueue(task1, 0)->wait();
        cuda::GPUPool::instance()->enqueue(task2, 1)->wait();
    }

    return a;
}
#endif

DenseVectorContinuousBase<float> & Difference<tags::GPU::CUDA>::value(DenseVectorContinuousBase<float> & r,
        const DenseVectorContinuousBase<float> & a,
        const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When subtracting DenseVectorContinuousBase<float> from DenseVectorContinuousBase<float> (CUDA):");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());
    if (r.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), r.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::difference_two_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaDifference3DVfloat task(r, a, b, blocksize);
        task();
    }
    else
    {
        cudaDifference3DVfloat task(r, a, b, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0)->wait();
    }

    return r;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVectorContinuousBase<double> & Difference<tags::GPU::CUDA>::value(DenseVectorContinuousBase<double> & r,
        const DenseVectorContinuousBase<double> & a,
        const DenseVectorContinuousBase<double> & b)
{
    CONTEXT("When subtracting DenseVectorContinuousBase<double> from DenseVectorContinuousBase<double> (CUDA):");
    PROFILER_START("Difference DV double tags::GPU::CUDA");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());
    if (r.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), r.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::difference_two_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaDifference3DVdouble task(r, a, b, blocksize);
        task();
    }
    else
    {
        cudaDifference3DVdouble task(r, a, b, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0)->wait();
    }

    PROFILER_STOP("Difference DV double tags::GPU::CUDA");
    return r;
}
#endif

DenseVectorContinuousBase<float> & Difference<tags::GPU::MultiCore::CUDA>::value(DenseVectorContinuousBase<float> & r,
        const DenseVectorContinuousBase<float> & a,
        const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When subtracting DenseVectorContinuousBase<float> from DenseVectorContinuousBase<float> (MC CUDA):");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());
    if (r.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), r.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::difference_two_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        throw InternalError("You should not run this operation within any MC CUDA op!");
    }
    else
    {
        DenseVectorRange<float> r1(r.range(r.size()/2, 0));
        DenseVectorRange<float> a1(a.range(a.size()/2, 0));
        DenseVectorRange<float> b1(b.range(b.size()/2, 0));
        cudaDifference3DVfloat task1(r1, a1, b1, blocksize);
        DenseVectorRange<float> r2(r.range(r.size()/2 + r.size()%2, r.size()/2));
        DenseVectorRange<float> a2(a.range(a.size()/2 + a.size()%2, a.size()/2));
        DenseVectorRange<float> b2(b.range(b.size()/2 + b.size()%2, b.size()/2));
        cudaDifference3DVfloat task2(r2, a2, b2, blocksize);
        cuda::GPUPool::instance()->enqueue(task1, 0)->wait();
        cuda::GPUPool::instance()->enqueue(task2, 1)->wait();
    }

    return r;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVectorContinuousBase<double> & Difference<tags::GPU::MultiCore::CUDA>::value(DenseVectorContinuousBase<double> & r,
        const DenseVectorContinuousBase<double> & a,
        const DenseVectorContinuousBase<double> & b)
{
    CONTEXT("When subtracting DenseVectorContinuousBase<double> from DenseVectorContinuousBase<double> (MC CUDA):");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());
    if (r.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), r.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::difference_two_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        throw InternalError("You should not run this operation within any MC CUDA op!");
    }
    else
    {
        DenseVectorRange<double> r1(r.range(r.size()/2, 0));
        DenseVectorRange<double> a1(a.range(a.size()/2, 0));
        DenseVectorRange<double> b1(b.range(b.size()/2, 0));
        cudaDifference3DVdouble task1(r1, a1, b1, blocksize);
        DenseVectorRange<double> r2(r.range(r.size()/2 + r.size()%2, r.size()/2));
        DenseVectorRange<double> a2(a.range(a.size()/2 + a.size()%2, a.size()/2));
        DenseVectorRange<double> b2(b.range(b.size()/2 + b.size()%2, b.size()/2));
        cudaDifference3DVdouble task2(r2, a2, b2, blocksize);
        cuda::GPUPool::instance()->enqueue(task1, 0)->wait();
        cuda::GPUPool::instance()->enqueue(task2, 1)->wait();
    }

    return r;
}
#endif

DenseMatrix<float> & Difference<tags::GPU::CUDA>::value(DenseMatrix<float> & a, const DenseMatrix<float> & b)
{
    CONTEXT("When subtracting DenseMatrix<float> from DenseMatrix<float> (CUDA):");

    if (a.columns() != b.columns())
    {
        throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
    }

    if (a.rows() != b.rows())
    {
        throw MatrixRowsDoNotMatch(b.rows(), a.rows());
    }

    unsigned long blocksize(Configuration::instance()->get_value("cuda::difference_two_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaDifferenceDMfloat task(a, b, blocksize);
        task();
    }
    else
    {
        cudaDifferenceDMfloat task(a, b, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0)->wait();
    }

    return a;
}
