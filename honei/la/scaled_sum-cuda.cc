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

#include <honei/la/scaled_sum.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/backends/cuda/gpu_pool.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/util/configuration.hh>


using namespace honei;

namespace
{
    class cudaScaledSumDVfloat
    {
        private:
            DenseVectorContinuousBase<float> & x;
            const DenseVectorContinuousBase<float> & y;
            float b;
            unsigned long blocksize;
        public:
            cudaScaledSumDVfloat(DenseVectorContinuousBase<float> & x, const DenseVectorContinuousBase<float> & y, float b, unsigned long blocksize) :
                x(x),
                y(y),
                b(b),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * x_gpu (x.lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * y_gpu (y.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                cuda_scaled_sum_three_float_s(x_gpu, x_gpu, y_gpu, b, x.size(), blocksize);
                y.unlock(lm_read_only);
                x.unlock(lm_read_and_write);
            }
    };

    class cudaScaledSumDVdouble
    {
        private:
            DenseVectorContinuousBase<double> & x;
            const DenseVectorContinuousBase<double> & y;
            double b;
            unsigned long blocksize;
        public:
            cudaScaledSumDVdouble(DenseVectorContinuousBase<double> & x, const DenseVectorContinuousBase<double> & y, double b, unsigned long blocksize) :
                x(x),
                y(y),
                b(b),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * x_gpu (x.lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * y_gpu (y.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                cuda_scaled_sum_three_double_s(x_gpu, x_gpu, y_gpu, b, x.size(), blocksize);
                y.unlock(lm_read_only);
                x.unlock(lm_read_and_write);
            }
    };

    class cudaScaledSum3sDVfloat
    {
        private:
            DenseVectorContinuousBase<float> & result;
            const DenseVectorContinuousBase<float> & x;
            const DenseVectorContinuousBase<float> & y;
            float b;
            unsigned long blocksize;
        public:
            cudaScaledSum3sDVfloat(DenseVectorContinuousBase<float> & result, const DenseVectorContinuousBase<float> & x, const DenseVectorContinuousBase<float> & y, float b, unsigned long blocksize) :
                result(result),
                x(x),
                y(y),
                b(b),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * result_gpu (result.lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * x_gpu (x.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * y_gpu (y.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                cuda_scaled_sum_three_float_s(result_gpu, x_gpu, y_gpu, b, x.size(), blocksize);
                y.unlock(lm_read_only);
                result.unlock(lm_write_only);
                x.unlock(lm_read_only);
            }
    };

    class cudaScaledSum3sDVdouble
    {
        private:
            DenseVectorContinuousBase<double> & result;
            const DenseVectorContinuousBase<double> & x;
            const DenseVectorContinuousBase<double> & y;
            double b;
            unsigned long blocksize;
        public:
            cudaScaledSum3sDVdouble(DenseVectorContinuousBase<double> & result, const DenseVectorContinuousBase<double> & x, const DenseVectorContinuousBase<double> & y, double b, unsigned long blocksize) :
                result(result),
                x(x),
                y(y),
                b(b),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * result_gpu (result.lock(lm_write_only, tags::GPU::CUDA::memory_value));
                void * x_gpu (x.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * y_gpu (y.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                cuda_scaled_sum_three_double_s(result_gpu, x_gpu, y_gpu, b, x.size(), blocksize);
                y.unlock(lm_read_only);
                result.unlock(lm_write_only);
                x.unlock(lm_read_only);
            }
    };

    class cudaScaledSum3DVfloat
    {
        private:
            DenseVectorContinuousBase<float> & a;
            const DenseVectorContinuousBase<float> & b;
            const DenseVectorContinuousBase<float> & c;
            unsigned long blocksize;
        public:
            cudaScaledSum3DVfloat(DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b, const DenseVectorContinuousBase<float> & c, unsigned long blocksize) :
                a(a),
                b(b),
                c(c),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * a_gpu (a.lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * b_gpu (b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * c_gpu (c.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                cuda_scaled_sum_three_float(a_gpu, b_gpu, c_gpu, a.size(), blocksize);
                c.unlock(lm_read_only);
                b.unlock(lm_read_only);
                a.unlock(lm_read_and_write);
            }
    };

    class cudaScaledSum3DVdouble
    {
        private:
            DenseVectorContinuousBase<double> & a;
            const DenseVectorContinuousBase<double> & b;
            const DenseVectorContinuousBase<double> & c;
            unsigned long blocksize;
        public:
            cudaScaledSum3DVdouble(DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b, const DenseVectorContinuousBase<double> & c, unsigned long blocksize) :
                a(a),
                b(b),
                c(c),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * a_gpu (a.lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * b_gpu (b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * c_gpu (c.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                cuda_scaled_sum_three_double(a_gpu, b_gpu, c_gpu, a.size(), blocksize);
                c.unlock(lm_read_only);
                b.unlock(lm_read_only);
                a.unlock(lm_read_and_write);
            }
    };
}

DenseVectorContinuousBase<float> & ScaledSum<tags::GPU::CUDA>::value(DenseVectorContinuousBase<float> & x,
        const DenseVectorContinuousBase<float> & y, float b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<float> (CUDA):");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::scaled_sum_two_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaScaledSumDVfloat task(x, y, b, blocksize);
        task();
    }
    else
    {
        cudaScaledSumDVfloat task(x, y, b, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0)->wait();
    }

    return x;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVectorContinuousBase<double> & ScaledSum<tags::GPU::CUDA>::value(DenseVectorContinuousBase<double> & x,
        const DenseVectorContinuousBase<double> & y, double b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<double> (CUDA):");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::scaled_sum_two_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaScaledSumDVdouble task(x, y, b, blocksize);
        task();
    }
    else
    {
        cudaScaledSumDVdouble task(x, y, b, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0)->wait();
    }

    return x;
}
#endif

DenseVectorContinuousBase<float> & ScaledSum<tags::GPU::MultiCore::CUDA>::value(DenseVectorContinuousBase<float> & x,
        const DenseVectorContinuousBase<float> & y, float b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<float> (MC CUDA):");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::scaled_sum_two_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        throw InternalError("You should not run this operation within any MC CUDA op!");
    }
    else
    {
        DenseVectorRange<float> x1(x.range(x.size()/2, 0));
        DenseVectorRange<float> y1(y.range(y.size()/2, 0));
        cudaScaledSumDVfloat task1(x1, y1, b, blocksize);
        DenseVectorRange<float> x2(x.range(x.size()/2 + x.size()%2, x.size()/2));
        DenseVectorRange<float> y2(y.range(y.size()/2 + y.size()%2, y.size()/2));
        cudaScaledSumDVfloat task2(x2, y2, b, blocksize);
        cuda::GPUPool::instance()->enqueue(task1, 0)->wait();
        cuda::GPUPool::instance()->enqueue(task2, 1)->wait();
    }

    return x;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVectorContinuousBase<double> & ScaledSum<tags::GPU::MultiCore::CUDA>::value(DenseVectorContinuousBase<double> & x,
        const DenseVectorContinuousBase<double> & y, double b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<double> (MC CUDA):");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::scaled_sum_two_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        throw InternalError("You should not run this operation within any MC CUDA op!");
    }
    else
    {
        DenseVectorRange<double> x1(x.range(x.size()/2, 0));
        DenseVectorRange<double> y1(y.range(y.size()/2, 0));
        cudaScaledSumDVdouble task1(x1, y1, b, blocksize);
        DenseVectorRange<double> x2(x.range(x.size()/2 + x.size()%2, x.size()/2));
        DenseVectorRange<double> y2(y.range(y.size()/2 + y.size()%2, y.size()/2));
        cudaScaledSumDVdouble task2(x2, y2, b, blocksize);
        cuda::GPUPool::instance()->enqueue(task1, 0)->wait();
        cuda::GPUPool::instance()->enqueue(task2, 1)->wait();
    }

    return x;
}
#endif

DenseVectorContinuousBase<float> & ScaledSum<tags::GPU::CUDA>::value(DenseVectorContinuousBase<float> & result, const DenseVectorContinuousBase<float> & x,
        const DenseVectorContinuousBase<float> & y, float b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<float> (CUDA):");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::scaled_sum_two_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaScaledSum3sDVfloat task(result, x, y, b, blocksize);
        task();
    }
    else
    {
        cudaScaledSum3sDVfloat task(result, x, y, b, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0)->wait();
    }

    return result;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVectorContinuousBase<double> & ScaledSum<tags::GPU::CUDA>::value(DenseVectorContinuousBase<double> & result, const DenseVectorContinuousBase<double> & x,
        const DenseVectorContinuousBase<double> & y, double b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<double> (CUDA):");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::scaled_sum_two_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaScaledSum3sDVdouble task(result, x, y, b, blocksize);
        task();
    }
    else
    {
        cudaScaledSum3sDVdouble task(result, x, y, b, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0)->wait();
    }

    return result;
}
#endif

DenseVectorContinuousBase<float> & ScaledSum<tags::GPU::MultiCore::CUDA>::value(DenseVectorContinuousBase<float> & result, const DenseVectorContinuousBase<float> & x,
        const DenseVectorContinuousBase<float> & y, float b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<float> (MC CUDA):");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::scaled_sum_two_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        throw InternalError("You should not run this operation within any MC CUDA op!");
    }
    else
    {
        DenseVectorRange<float> x1(x.range(x.size()/2, 0));
        DenseVectorRange<float> y1(y.range(y.size()/2, 0));
        DenseVectorRange<float> r1(result.range(result.size()/2, 0));
        cudaScaledSum3sDVfloat task1(r1, x1, y1, b, blocksize);
        DenseVectorRange<float> x2(x.range(x.size()/2 + x.size()%2, x.size()/2));
        DenseVectorRange<float> y2(y.range(y.size()/2 + y.size()%2, y.size()/2));
        DenseVectorRange<float> r2(result.range(result.size()/2 + result.size()%2, result.size()/2));
        cudaScaledSum3sDVfloat task2(r2, x2, y2, b, blocksize);
        cuda::GPUPool::instance()->enqueue(task1, 0)->wait();
        cuda::GPUPool::instance()->enqueue(task2, 1)->wait();
    }

    return result;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVectorContinuousBase<double> & ScaledSum<tags::GPU::MultiCore::CUDA>::value(DenseVectorContinuousBase<double> & result, const DenseVectorContinuousBase<double> & x,
        const DenseVectorContinuousBase<double> & y, double b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<double> (MC CUDA):");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::scaled_sum_two_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        throw InternalError("You should not run this operation within any MC CUDA op!");
    }
    else
    {
        DenseVectorRange<double> x1(x.range(x.size()/2, 0));
        DenseVectorRange<double> y1(y.range(y.size()/2, 0));
        DenseVectorRange<double> r1(result.range(result.size()/2, 0));
        cudaScaledSum3sDVdouble task1(r1, x1, y1, b, blocksize);
        DenseVectorRange<double> x2(x.range(x.size()/2 + x.size()%2, x.size()/2));
        DenseVectorRange<double> y2(y.range(y.size()/2 + y.size()%2, y.size()/2));
        DenseVectorRange<double> r2(result.range(result.size()/2 + result.size()%2, result.size()/2));
        cudaScaledSum3sDVdouble task2(r2, x2, y2, b, blocksize);
        cuda::GPUPool::instance()->enqueue(task1, 0)->wait();
        cuda::GPUPool::instance()->enqueue(task2, 1)->wait();
    }

    return result;
}
#endif

DenseVectorContinuousBase<float> & ScaledSum<tags::GPU::CUDA>::value(DenseVectorContinuousBase<float> & a,
        const DenseVectorContinuousBase<float> & b, const DenseVectorContinuousBase<float> & c)
{
    CONTEXT("When calculating ScaledSum (DenseVectorContinuousBase<float>, DenseVectorContinuousBase<float>, "
            "DenseVectorContinuousBase<float>) (CUDA):");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    if (a.size() != c.size())
        throw VectorSizeDoesNotMatch(c.size(), a.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::scaled_sum_three_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaScaledSum3DVfloat task(a, b, c, blocksize);
        task();
    }
    else
    {
        cudaScaledSum3DVfloat task(a, b, c, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0)->wait();
    }

    return a;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVectorContinuousBase<double> & ScaledSum<tags::GPU::CUDA>::value(DenseVectorContinuousBase<double> & a,
        const DenseVectorContinuousBase<double> & b, const DenseVectorContinuousBase<double> & c)
{
    CONTEXT("When calculating ScaledSum (DenseVectorContinuousBase<double>, DenseVectorContinuousBase<double>, "
            "DenseVectorContinuousBase<double>) (CUDA):");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    if (a.size() != c.size())
        throw VectorSizeDoesNotMatch(c.size(), a.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::scaled_sum_three_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaScaledSum3DVdouble task(a, b, c, blocksize);
        task();
    }
    else
    {
        cudaScaledSum3DVdouble task(a, b, c, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0)->wait();
    }

    return a;
}
#endif

DenseVectorContinuousBase<float> & ScaledSum<tags::GPU::MultiCore::CUDA>::value(DenseVectorContinuousBase<float> & x,
        const DenseVectorContinuousBase<float> & y, const DenseVectorContinuousBase<float> & z)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<float> (MC CUDA):");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());
    if (x.size() != z.size())
        throw VectorSizeDoesNotMatch(x.size(), z.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::scaled_sum_three_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        throw InternalError("You should not run this operation within any MC CUDA op!");
    }
    else
    {
        DenseVectorRange<float> x1(x.range(x.size()/2, 0));
        DenseVectorRange<float> y1(y.range(y.size()/2, 0));
        DenseVectorRange<float> z1(z.range(z.size()/2, 0));
        cudaScaledSum3DVfloat task1(x1, y1, z1, blocksize);
        DenseVectorRange<float> x2(x.range(x.size()/2 + x.size()%2, x.size()/2));
        DenseVectorRange<float> y2(y.range(y.size()/2 + y.size()%2, y.size()/2));
        DenseVectorRange<float> z2(z.range(z.size()/2 + z.size()%2, z.size()/2));
        cudaScaledSum3DVfloat task2(x2, y2, z2, blocksize);
        cuda::GPUPool::instance()->enqueue(task1, 0)->wait();
        cuda::GPUPool::instance()->enqueue(task2, 1)->wait();
    }

    return x;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVectorContinuousBase<double> & ScaledSum<tags::GPU::MultiCore::CUDA>::value(DenseVectorContinuousBase<double> & x,
        const DenseVectorContinuousBase<double> & y, const DenseVectorContinuousBase<double> & z)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<double> (MC CUDA):");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());
    if (x.size() != z.size())
        throw VectorSizeDoesNotMatch(x.size(), z.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::scaled_sum_three_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        throw InternalError("You should not run this operation within any MC CUDA op!");
    }
    else
    {
        DenseVectorRange<double> x1(x.range(x.size()/2, 0));
        DenseVectorRange<double> y1(y.range(y.size()/2, 0));
        DenseVectorRange<double> z1(z.range(z.size()/2, 0));
        cudaScaledSum3DVdouble task1(x1, y1, z1, blocksize);
        DenseVectorRange<double> x2(x.range(x.size()/2 + x.size()%2, x.size()/2));
        DenseVectorRange<double> y2(y.range(y.size()/2 + y.size()%2, y.size()/2));
        DenseVectorRange<double> z2(z.range(z.size()/2 + z.size()%2, z.size()/2));
        cudaScaledSum3DVdouble task2(x2, y2, z2, blocksize);
        cuda::GPUPool::instance()->enqueue(task1, 0)->wait();
        cuda::GPUPool::instance()->enqueue(task2, 1)->wait();
    }

    return x;
}
#endif
