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

#include <honei/la/scaled_sum.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/util/configuration.hh>


using namespace honei;

DenseVectorContinuousBase<float> & ScaledSum<tags::GPU::CUDA>::value(DenseVectorContinuousBase<float> & x,
        const DenseVectorContinuousBase<float> & y, float b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<float> (CUDA):");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::scaled_sum_two_float", 128ul));

    void * x_gpu (x.lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
    void * y_gpu (y.lock(lm_read_only, tags::GPU::CUDA::memory_value));
    cuda_scaled_sum_three_float_s(x_gpu, x_gpu, y_gpu, b, x.size(), blocksize);
    y.unlock(lm_read_only);
    x.unlock(lm_read_and_write);

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

    void * x_gpu (x.lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
    void * y_gpu (y.lock(lm_read_only, tags::GPU::CUDA::memory_value));
    cuda_scaled_sum_three_double_s(x_gpu, x_gpu, y_gpu, b, x.size(), blocksize);
    y.unlock(lm_read_only);
    x.unlock(lm_read_and_write);

    return x;
}
#endif

DenseVectorContinuousBase<float> & ScaledSum<tags::GPU::CUDA>::value(DenseVectorContinuousBase<float> & result, DenseVectorContinuousBase<float> & x,
        const DenseVectorContinuousBase<float> & y, float b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<float> (CUDA):");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::scaled_sum_two_float", 128ul));

    void * result_gpu (result.lock(lm_write_only, tags::GPU::CUDA::memory_value));
    void * x_gpu (x.lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * y_gpu (y.lock(lm_read_only, tags::GPU::CUDA::memory_value));
    cuda_scaled_sum_three_float_s(result_gpu, x_gpu, y_gpu, b, x.size(), blocksize);
    y.unlock(lm_read_only);
    result.unlock(lm_write_only);
    x.unlock(lm_read_only);

    return result;
}

#ifdef HONEI_CUDA_DOUBLE
DenseVectorContinuousBase<double> & ScaledSum<tags::GPU::CUDA>::value(DenseVectorContinuousBase<double> & result, DenseVectorContinuousBase<double> & x,
        const DenseVectorContinuousBase<double> & y, double b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<double> (CUDA):");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::scaled_sum_two_double", 128ul));

    void * result_gpu (result.lock(lm_write_only, tags::GPU::CUDA::memory_value));
    void * x_gpu (x.lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * y_gpu (y.lock(lm_read_only, tags::GPU::CUDA::memory_value));
    cuda_scaled_sum_three_double_s(result_gpu, x_gpu, y_gpu, b, x.size(), blocksize);
    y.unlock(lm_read_only);
    result.unlock(lm_write_only);
    x.unlock(lm_read_only);

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

    void * a_gpu (a.lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
    void * b_gpu (b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * c_gpu (c.lock(lm_read_only, tags::GPU::CUDA::memory_value));
    cuda_scaled_sum_three_float(a_gpu, b_gpu, c_gpu, a.size(), blocksize);
    c.unlock(lm_read_only);
    b.unlock(lm_read_only);
    a.unlock(lm_read_and_write);

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

    void * a_gpu (a.lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
    void * b_gpu (b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * c_gpu (c.lock(lm_read_only, tags::GPU::CUDA::memory_value));
    cuda_scaled_sum_three_double(a_gpu, b_gpu, c_gpu, a.size(), blocksize);
    c.unlock(lm_read_only);
    b.unlock(lm_read_only);
    a.unlock(lm_read_and_write);

    return a;
}
#endif
