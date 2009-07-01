/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
#include <honei/util/configuration.hh>


using namespace honei;


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
        void * a_gpu(a.lock(lm_read_only, tags::GPU::CUDA::memory_value));
        result = cuda_norm_l2_one_float(a_gpu, a.size(), blocksize, gridsize);
        a.unlock(lm_read_only);
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
        void * a_gpu(a.lock(lm_read_only, tags::GPU::CUDA::memory_value));
        result = cuda_norm_l2_one_float(a_gpu, a.size(), blocksize, gridsize);
        a.unlock(lm_read_only);
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
        void * a_gpu(a.lock(lm_read_only, tags::GPU::CUDA::memory_value));
        result = cuda_norm_l2_one_double(a_gpu, a.size(), blocksize, gridsize);
        a.unlock(lm_read_only);
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
        void * a_gpu(a.lock(lm_read_only, tags::GPU::CUDA::memory_value));
        result = cuda_norm_l2_one_double(a_gpu, a.size(), blocksize, gridsize);
        a.unlock(lm_read_only);
    }

    return sqrt(result);
}
#endif
