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

#include <honei/la/norm.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/backends/cuda/gpu_pool.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/util/configuration.hh>


using namespace honei;

namespace
{
    class cudaNormL2oneDVfloat
    {
        private:
            const DenseVectorContinuousBase<float> & a;
            float * result;
            unsigned long blocksize;
            unsigned long gridsize;
        public:
            cudaNormL2oneDVfloat(const DenseVectorContinuousBase<float> & a, float * result, unsigned long blocksize, unsigned long gridsize) :
                a(a),
                result(result),
                blocksize(blocksize),
                gridsize(gridsize)
            {
            }

            void operator() ()
            {
                void * a_gpu(a.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                *result = cuda_norm_l2_one_float(a_gpu, a.size(), blocksize, gridsize);
                a.unlock(lm_read_only);
            }
    };

    class cudaNormL2oneDVdouble
    {
        private:
            const DenseVectorContinuousBase<double> & a;
            double * result;
            unsigned long blocksize;
            unsigned long gridsize;
        public:
            cudaNormL2oneDVdouble(const DenseVectorContinuousBase<double> & a, double * result, unsigned long blocksize, unsigned long gridsize) :
                a(a),
                result(result),
                blocksize(blocksize),
                gridsize(gridsize)
            {
            }

            void operator() ()
            {
                void * a_gpu(a.lock(lm_read_only, tags::GPU::CUDA::memory_value));
                *result = cuda_norm_l2_one_double(a_gpu, a.size(), blocksize, gridsize);
                a.unlock(lm_read_only);
            }
    };
}

float Norm<vnt_l_two, false, tags::GPU::CUDA>::value(const DenseVectorContinuousBase<float> & a)
{
    CONTEXT("When calculating L2 norm (false) of DenseVectorContinuousBase<float> (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::norm_l2_one_float", 128ul));
    unsigned long gridsize(Configuration::instance()->get_value("cuda::norm_l2_one_float_grid", 16ul));

    float result (0.);

    if (a.size() < gridsize * blocksize)
    {
        a.lock(lm_read_only);
        for (unsigned long i(0) ; i < a.size() ; ++i)
        {
            result += a[i] * a[i];
        }
        a.unlock(lm_read_only);
    }
    else
    {
        if (! cuda::GPUPool::instance()->idle())
        {
            cudaNormL2oneDVfloat task(a, &result, blocksize, gridsize);
            task();
        }
        else
        {
            cudaNormL2oneDVfloat task(a, &result, blocksize, gridsize);
            cuda::GPUPool::instance()->enqueue(task, 0)->wait();
        }
    }

    return result;
}

float Norm<vnt_l_two, true, tags::GPU::CUDA>::value(const DenseVectorContinuousBase<float> & a)
{
    CONTEXT("When calculating L2 norm (true) of DenseVectorContinuousBase<float> (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::norm_l2_one_float", 128ul));
    unsigned long gridsize(Configuration::instance()->get_value("cuda::norm_l2_one_float_grid", 16ul));

    float result (0.);

    if (a.size() < gridsize * blocksize)
    {
        a.lock(lm_read_only);
        for (unsigned long i(0) ; i < a.size() ; ++i)
        {
            result += a[i] * a[i];
        }
        a.unlock(lm_read_only);
    }
    else
    {
        if (! cuda::GPUPool::instance()->idle())
        {
            cudaNormL2oneDVfloat task(a, &result, blocksize, gridsize);
            task();
        }
        else
        {
            cudaNormL2oneDVfloat task(a, &result, blocksize, gridsize);
            cuda::GPUPool::instance()->enqueue(task, 0)->wait();
        }
    }

    return sqrt(result);
}

#ifdef HONEI_CUDA_DOUBLE
double Norm<vnt_l_two, false, tags::GPU::CUDA>::value(const DenseVectorContinuousBase<double> & a)
{
    CONTEXT("When calculating L2 norm (false) of DenseVectorContinuousBase<double> (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::norm_l2_one_double", 128ul));
    unsigned long gridsize(Configuration::instance()->get_value("cuda::norm_l2_one_double", 16ul));

    double result (0.);

    if (a.size() < gridsize * blocksize)
    {
        a.lock(lm_read_only);
        for (unsigned long i(0) ; i < a.size() ; ++i)
        {
            result += a[i] * a[i];
        }
        a.unlock(lm_read_only);
    }
    else
    {
        if (! cuda::GPUPool::instance()->idle())
        {
            cudaNormL2oneDVdouble task(a, &result, blocksize, gridsize);
            task();
        }
        else
        {
            cudaNormL2oneDVdouble task(a, &result, blocksize, gridsize);
            cuda::GPUPool::instance()->enqueue(task, 0)->wait();
        }
    }

    return result;
}

double Norm<vnt_l_two, true, tags::GPU::CUDA>::value(const DenseVectorContinuousBase<double> & a)
{
    CONTEXT("When calculating L2 norm (true) of DenseVectorContinuousBase<double> (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::norm_l2_one_double", 128ul));
    unsigned long gridsize(Configuration::instance()->get_value("cuda::norm_l2_one_double", 16ul));

    double result (0.);

    if (a.size() < gridsize * blocksize)
    {
        a.lock(lm_read_only);
        for (unsigned long i(0) ; i < a.size() ; ++i)
        {
            result += a[i] * a[i];
        }
        a.unlock(lm_read_only);
    }
    else
    {
        if (! cuda::GPUPool::instance()->idle())
        {
            cudaNormL2oneDVdouble task(a, &result, blocksize, gridsize);
            task();
        }
        else
        {
            cudaNormL2oneDVdouble task(a, &result, blocksize, gridsize);
            cuda::GPUPool::instance()->enqueue(task, 0)->wait();
        }
    }

    return sqrt(result);
}
#endif
