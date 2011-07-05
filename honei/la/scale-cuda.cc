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

#include <honei/la/scale.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/backends/cuda/gpu_pool.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/util/configuration.hh>
#include <honei/util/profiler.hh>


using namespace honei;

namespace
{
    class cudaScaleDVfloat
    {
        private:
            DenseVectorContinuousBase<float> & a;
            float scal;
            unsigned long blocksize;
        public:
            cudaScaleDVfloat(DenseVectorContinuousBase<float> & a, float scal, unsigned long blocksize) :
                a(a),
                scal(scal),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * a_gpu (a.lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                cuda_scale_one_float(a_gpu, scal, a.size(), blocksize);
                a.unlock(lm_read_and_write);
            }
    };

    class cudaScaleDVdouble
    {
        private:
            DenseVectorContinuousBase<double> & a;
            double scal;
            unsigned long blocksize;
        public:
            cudaScaleDVdouble(DenseVectorContinuousBase<double> & a, double scal, unsigned long blocksize) :
                a(a),
                scal(scal),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * a_gpu (a.lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                cuda_scale_one_double(a_gpu, scal, a.size(), blocksize);
                a.unlock(lm_read_and_write);
            }
    };

    class cudaScaleDMfloat
    {
        private:
            DenseMatrix<float> & a;
            float scal;
            unsigned long blocksize;
        public:
            cudaScaleDMfloat(DenseMatrix<float> & a, float scal, unsigned long blocksize) :
                a(a),
                scal(scal),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * a_gpu (a.lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                cuda_scale_one_float(a_gpu, scal, a.size(), blocksize);
                a.unlock(lm_read_and_write);
            }
    };
}


DenseVectorContinuousBase<float> & Scale<tags::GPU::CUDA>::value(DenseVectorContinuousBase<float> & x, const float a)
{
    CONTEXT("When scaling DenseVectorContinuousBase<float> by float (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::scale_one_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaScaleDVfloat task(x, a, blocksize);
        task();
    }
    else
    {
        cudaScaleDVfloat task(x, a, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }

    return x;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVectorContinuousBase<double> & Scale<tags::GPU::CUDA>::value(DenseVectorContinuousBase<double> & x, const double a)
{
    CONTEXT("When scaling DenseVectorContinuousBase<double> by double (CUDA):");
    PROFILER_START("Scale DV double tags::GPU::CUDA");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::scale_one_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaScaleDVdouble task(x, a, blocksize);
        task();
    }
    else
    {
        cudaScaleDVdouble task(x, a, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }

    PROFILER_STOP("Scale DV double tags::GPU::CUDA");
    return x;
}
#endif


DenseVectorContinuousBase<float> & Scale<tags::GPU::MultiCore::CUDA>::value(DenseVectorContinuousBase<float> & x, const float a)
{
    CONTEXT("When scaling DenseVectorContinuousBase<float> by float (MC CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::scale_one_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        throw InternalError("You should not run this operation within any MC CUDA op!");
    }
    else
    {
        DenseVectorRange<float> x1(x.range(x.size()/2, 0));
        cudaScaleDVfloat task1(x1, a, blocksize);
        DenseVectorRange<float> x2(x.range(x.size()/2 + x.size()%2, x.size()/2));
        cudaScaleDVfloat task2(x2, a, blocksize);
        cuda::GPUPool::instance()->enqueue(task1, 0).wait();
        cuda::GPUPool::instance()->enqueue(task2, 1).wait();
    }

    return x;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVectorContinuousBase<double> & Scale<tags::GPU::MultiCore::CUDA>::value(DenseVectorContinuousBase<double> & x, const double a)
{
    CONTEXT("When scaling DenseVectorContinuousBase<double> by double (MC CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::scale_one_double", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        throw InternalError("You should not run this operation within any MC CUDA op!");
    }
    else
    {
        DenseVectorRange<double> x1(x.range(x.size()/2, 0));
        cudaScaleDVdouble task1(x1, a, blocksize);
        DenseVectorRange<double> x2(x.range(x.size()/2 + x.size()%2, x.size()/2));
        cudaScaleDVdouble task2(x2, a, blocksize);
        cuda::GPUPool::instance()->enqueue(task1, 0).wait();
        cuda::GPUPool::instance()->enqueue(task2, 1).wait();
    }

    return x;
}
#endif

DenseMatrix<float> & Scale<tags::GPU::CUDA>::value(DenseMatrix<float> & x, const float a)
{
    CONTEXT("When scaling DenseMatrix<float> by float (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::scale_one_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaScaleDMfloat task(x, a, blocksize);
        task();
    }
    else
    {
        cudaScaleDMfloat task(x, a, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }

    return x;
}
