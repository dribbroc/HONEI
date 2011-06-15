/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
 *
 * This file is part of the LA C++ library. LibLa is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibLa is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/la/scaled_sum.hh>
#include <honei/backends/opencl/operations.hh>
#include <honei/util/profiler.hh>


using namespace honei;

DenseVectorContinuousBase<float> & ScaledSum<tags::OpenCL::CPU>::value(DenseVectorContinuousBase<float> & x,
        const DenseVectorContinuousBase<float> & y, float b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<float> (OpenCL CPU):");
    PROFILER_START("ScaledSum tags::OpenCL::CPU");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());

    void * x_cl(x.lock(lm_read_and_write, tags::OpenCL::CPU::memory_value));
    void * y_cl(y.lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    opencl::scaled_sum<float>(x_cl, x_cl, y_cl, b, x.size(), CL_DEVICE_TYPE_CPU, "scaled_sum_three_float_s");
    x.unlock(lm_read_and_write);
    y.unlock(lm_read_only);

    PROFILER_STOP("ScaledSum tags::OpenCL::CPU");
    return x;
}

DenseVectorContinuousBase<double> & ScaledSum<tags::OpenCL::CPU>::value(DenseVectorContinuousBase<double> & x,
        const DenseVectorContinuousBase<double> & y, double b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<double> (OpenCL CPU):");
    PROFILER_START("ScaledSum tags::OpenCL::CPU");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());

    void * x_cl(x.lock(lm_read_and_write, tags::OpenCL::CPU::memory_value));
    void * y_cl(y.lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    opencl::scaled_sum<double>(x_cl, x_cl, y_cl, b, x.size(), CL_DEVICE_TYPE_CPU, "scaled_sum_three_double_s");
    x.unlock(lm_read_and_write);
    y.unlock(lm_read_only);

    PROFILER_STOP("ScaledSum tags::OpenCL::CPU");
    return x;
}

DenseVectorContinuousBase<float> & ScaledSum<tags::OpenCL::GPU>::value(DenseVectorContinuousBase<float> & x,
        const DenseVectorContinuousBase<float> & y, float b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<float> (OpenCL GPU):");
    PROFILER_START("ScaledSum tags::OpenCL::GPU");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());

    void * x_cl(x.lock(lm_read_and_write, tags::OpenCL::GPU::memory_value));
    void * y_cl(y.lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    opencl::scaled_sum<float>(x_cl, x_cl, y_cl, b, x.size(), CL_DEVICE_TYPE_GPU, "scaled_sum_three_float_s");
    x.unlock(lm_read_and_write);
    y.unlock(lm_read_only);

    PROFILER_STOP("ScaledSum tags::OpenCL::GPU");
    return x;
}

DenseVectorContinuousBase<double> & ScaledSum<tags::OpenCL::GPU>::value(DenseVectorContinuousBase<double> & x,
        const DenseVectorContinuousBase<double> & y, double b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<double> (OpenCL GPU):");
    PROFILER_START("ScaledSum tags::OpenCL::GPU");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());

    void * x_cl(x.lock(lm_read_and_write, tags::OpenCL::GPU::memory_value));
    void * y_cl(y.lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    opencl::scaled_sum<double>(x_cl, x_cl, y_cl, b, x.size(), CL_DEVICE_TYPE_GPU, "scaled_sum_three_double_s");
    x.unlock(lm_read_and_write);
    y.unlock(lm_read_only);

    PROFILER_STOP("ScaledSum tags::OpenCL::GPU");
    return x;
}

DenseVectorContinuousBase<float> & ScaledSum<tags::OpenCL::CPU>::value(DenseVectorContinuousBase<float> & r, const DenseVectorContinuousBase<float> & x,
        const DenseVectorContinuousBase<float> & y, float b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<float> (OpenCL CPU):");
    PROFILER_START("ScaledSum tags::OpenCL::CPU");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());
    if (r.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), r.size());

    void * r_cl(r.lock(lm_write_only, tags::OpenCL::CPU::memory_value));
    void * x_cl(x.lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    void * y_cl(y.lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    opencl::scaled_sum<float>(r_cl, x_cl, y_cl, b, x.size(), CL_DEVICE_TYPE_CPU, "scaled_sum_three_float_s");
    r.unlock(lm_write_only);
    x.unlock(lm_read_only);
    y.unlock(lm_read_only);

    PROFILER_STOP("ScaledSum tags::OpenCL::CPU");
    return r;
}

DenseVectorContinuousBase<double> & ScaledSum<tags::OpenCL::CPU>::value(DenseVectorContinuousBase<double> & r, const DenseVectorContinuousBase<double> & x,
        const DenseVectorContinuousBase<double> & y, double b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<double> (OpenCL CPU):");
    PROFILER_START("ScaledSum tags::OpenCL::CPU");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());
    if (r.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), r.size());

    void * r_cl(r.lock(lm_write_only, tags::OpenCL::CPU::memory_value));
    void * x_cl(x.lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    void * y_cl(y.lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    opencl::scaled_sum<double>(r_cl, x_cl, y_cl, b, x.size(), CL_DEVICE_TYPE_CPU, "scaled_sum_three_double_s");
    r.unlock(lm_write_only);
    x.unlock(lm_read_only);
    y.unlock(lm_read_only);

    PROFILER_STOP("ScaledSum tags::OpenCL::CPU");
    return r;
}

DenseVectorContinuousBase<float> & ScaledSum<tags::OpenCL::GPU>::value(DenseVectorContinuousBase<float> & r, const DenseVectorContinuousBase<float> & x,
        const DenseVectorContinuousBase<float> & y, float b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<float> (OpenCL GPU):");
    PROFILER_START("ScaledSum tags::OpenCL::GPU");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());
    if (r.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), r.size());

    void * r_cl(r.lock(lm_write_only, tags::OpenCL::GPU::memory_value));
    void * x_cl(x.lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    void * y_cl(y.lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    opencl::scaled_sum<float>(r_cl, x_cl, y_cl, b, x.size(), CL_DEVICE_TYPE_GPU, "scaled_sum_three_float_s");
    r.unlock(lm_write_only);
    x.unlock(lm_read_only);
    y.unlock(lm_read_only);

    PROFILER_STOP("ScaledSum tags::OpenCL::GPU");
    return r;
}

DenseVectorContinuousBase<double> & ScaledSum<tags::OpenCL::GPU>::value(DenseVectorContinuousBase<double> & r, const DenseVectorContinuousBase<double> & x,
        const DenseVectorContinuousBase<double> & y, double b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<double> (OpenCL GPU):");
    PROFILER_START("ScaledSum tags::OpenCL::GPU");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());
    if (r.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), r.size());

    void * r_cl(r.lock(lm_write_only, tags::OpenCL::GPU::memory_value));
    void * x_cl(x.lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    void * y_cl(y.lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    opencl::scaled_sum<double>(r_cl, x_cl, y_cl, b, x.size(), CL_DEVICE_TYPE_GPU, "scaled_sum_three_double_s");
    r.unlock(lm_write_only);
    x.unlock(lm_read_only);
    y.unlock(lm_read_only);

    PROFILER_STOP("ScaledSum tags::OpenCL::GPU");
    return r;
}

DenseVectorContinuousBase<float> & ScaledSum<tags::OpenCL::CPU>::value(DenseVectorContinuousBase<float> & r, const DenseVectorContinuousBase<float> & x,
        const DenseVectorContinuousBase<float> & y)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<float> (OpenCL CPU):");
    PROFILER_START("ScaledSum tags::OpenCL::CPU");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());
    if (r.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), r.size());

    void * r_cl(r.lock(lm_read_and_write, tags::OpenCL::CPU::memory_value));
    void * x_cl(x.lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    void * y_cl(y.lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    opencl::scaled_sum<float>(r_cl, x_cl, y_cl, x.size(), CL_DEVICE_TYPE_CPU, "scaled_sum_three_float");
    r.unlock(lm_read_and_write);
    x.unlock(lm_read_only);
    y.unlock(lm_read_only);

    PROFILER_STOP("ScaledSum tags::OpenCL::CPU");
    return r;
}

DenseVectorContinuousBase<double> & ScaledSum<tags::OpenCL::CPU>::value(DenseVectorContinuousBase<double> & r, const DenseVectorContinuousBase<double> & x,
        const DenseVectorContinuousBase<double> & y)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<double> (OpenCL CPU):");
    PROFILER_START("ScaledSum tags::OpenCL::CPU");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());
    if (r.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), r.size());

    void * r_cl(r.lock(lm_read_and_write, tags::OpenCL::CPU::memory_value));
    void * x_cl(x.lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    void * y_cl(y.lock(lm_read_only, tags::OpenCL::CPU::memory_value));
    opencl::scaled_sum<double>(r_cl, x_cl, y_cl, x.size(), CL_DEVICE_TYPE_CPU, "scaled_sum_three_double");
    r.unlock(lm_read_and_write);
    x.unlock(lm_read_only);
    y.unlock(lm_read_only);

    PROFILER_STOP("ScaledSum tags::OpenCL::CPU");
    return r;
}

DenseVectorContinuousBase<float> & ScaledSum<tags::OpenCL::GPU>::value(DenseVectorContinuousBase<float> & r, const DenseVectorContinuousBase<float> & x,
        const DenseVectorContinuousBase<float> & y)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<float> (OpenCL GPU):");
    PROFILER_START("ScaledSum tags::OpenCL::GPU");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());
    if (r.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), r.size());

    void * r_cl(r.lock(lm_read_and_write, tags::OpenCL::GPU::memory_value));
    void * x_cl(x.lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    void * y_cl(y.lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    opencl::scaled_sum<float>(r_cl, x_cl, y_cl, x.size(), CL_DEVICE_TYPE_GPU, "scaled_sum_three_float");
    r.unlock(lm_read_and_write);
    x.unlock(lm_read_only);
    y.unlock(lm_read_only);

    PROFILER_STOP("ScaledSum tags::OpenCL::GPU");
    return r;
}

DenseVectorContinuousBase<double> & ScaledSum<tags::OpenCL::GPU>::value(DenseVectorContinuousBase<double> & r, const DenseVectorContinuousBase<double> & x,
        const DenseVectorContinuousBase<double> & y)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<double> (OpenCL GPU):");
    PROFILER_START("ScaledSum tags::OpenCL::GPU");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());
    if (r.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), r.size());

    void * r_cl(r.lock(lm_read_and_write, tags::OpenCL::GPU::memory_value));
    void * x_cl(x.lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    void * y_cl(y.lock(lm_read_only, tags::OpenCL::GPU::memory_value));
    opencl::scaled_sum<double>(r_cl, x_cl, y_cl, x.size(), CL_DEVICE_TYPE_GPU, "scaled_sum_three_double");
    r.unlock(lm_read_and_write);
    x.unlock(lm_read_only);
    y.unlock(lm_read_only);

    PROFILER_STOP("ScaledSum tags::OpenCL::GPU");
    return r;
}
