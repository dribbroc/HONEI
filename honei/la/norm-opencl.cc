/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2011 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
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

#include <honei/la/norm.hh>
#include <honei/backends/opencl/operations.hh>
#include <honei/util/profiler.hh>


using namespace honei;

float Norm<vnt_l_two, false, tags::OpenCL::CPU>::value(const DenseVectorContinuousBase<float> & x)
{
    CONTEXT("When calculating L2 norm (false) of DenseVectorContinuousBase<float> (OpenCL CPU):");
    PROFILER_START("Norm l2 false tags::OpenCL::CPU");

    float result(0);
    if (x.size() < 32 * 128)
    {
        x.lock(lm_read_only);
        for (unsigned long i(0) ; i < x.size() ; ++i)
        {
            result += x[i] * x[i];
        }
        x.unlock(lm_read_only);
    }
    else
    {
        void * x_cl(x.lock(lm_read_only, tags::OpenCL::CPU::memory_value));
        result = opencl::norm_l2_false<float>(x_cl, x.size(), CL_DEVICE_TYPE_CPU, "norm_l2_false_float");
        x.unlock(lm_read_only);
    }
    PROFILER_STOP("Norm l2 false tags::OpenCL::CPU");
    return result;
}

double Norm<vnt_l_two, false, tags::OpenCL::CPU>::value(const DenseVectorContinuousBase<double> & x)
{
    CONTEXT("When calculating L2 norm (false) of DenseVectorContinuousBase<double> (OpenCL CPU):");
    PROFILER_START("Norm l2 false tags::OpenCL::CPU");

    double result(0);
    if (x.size() < 32 * 128)
    {
        x.lock(lm_read_only);
        for (unsigned long i(0) ; i < x.size() ; ++i)
        {
            result += x[i] * x[i];
        }
        x.unlock(lm_read_only);
    }
    else
    {
        void * x_cl(x.lock(lm_read_only, tags::OpenCL::CPU::memory_value));
        result = opencl::norm_l2_false<double>(x_cl, x.size(), CL_DEVICE_TYPE_CPU, "norm_l2_false_double");
        x.unlock(lm_read_only);
    }
    PROFILER_STOP("Norm l2 false tags::OpenCL::CPU");
    return result;
}

float Norm<vnt_l_two, true, tags::OpenCL::CPU>::value(const DenseVectorContinuousBase<float> & x)
{
    CONTEXT("When calculating L2 norm (false) of DenseVectorContinuousBase<float> (OpenCL CPU):");
    PROFILER_START("Norm l2 true tags::OpenCL::CPU");

    float result(0);
    if (x.size() < 32 * 128)
    {
        x.lock(lm_read_only);
        for (unsigned long i(0) ; i < x.size() ; ++i)
        {
            result += x[i] * x[i];
        }
        x.unlock(lm_read_only);
    }
    else
    {
        void * x_cl(x.lock(lm_read_only, tags::OpenCL::CPU::memory_value));
        result = opencl::norm_l2_false<float>(x_cl, x.size(), CL_DEVICE_TYPE_CPU, "norm_l2_false_float");
        x.unlock(lm_read_only);
    }
    PROFILER_STOP("Norm l2 true tags::OpenCL::CPU");
    return sqrt(result);
}

double Norm<vnt_l_two, true, tags::OpenCL::CPU>::value(const DenseVectorContinuousBase<double> & x)
{
    CONTEXT("When calculating L2 norm (false) of DenseVectorContinuousBase<double> (OpenCL CPU):");
    PROFILER_START("Norm l2 true tags::OpenCL::CPU");

    double result(0);
    if (x.size() < 32 * 128)
    {
        x.lock(lm_read_only);
        for (unsigned long i(0) ; i < x.size() ; ++i)
        {
            result += x[i] * x[i];
        }
        x.unlock(lm_read_only);
    }
    else
    {
        void * x_cl(x.lock(lm_read_only, tags::OpenCL::CPU::memory_value));
        result = opencl::norm_l2_false<double>(x_cl, x.size(), CL_DEVICE_TYPE_CPU, "norm_l2_false_double");
        x.unlock(lm_read_only);
    }
    PROFILER_STOP("Norm l2 true tags::OpenCL::CPU");
    return sqrt(result);
}

float Norm<vnt_l_two, false, tags::OpenCL::GPU>::value(const DenseVectorContinuousBase<float> & x)
{
    CONTEXT("When calculating L2 norm (false) of DenseVectorContinuousBase<float> (OpenCL GPU):");
    PROFILER_START("Norm l2 false tags::OpenCL::GPU");

    float result(0);
    if (x.size() < 32 * 128)
    {
        x.lock(lm_read_only);
        for (unsigned long i(0) ; i < x.size() ; ++i)
        {
            result += x[i] * x[i];
        }
        x.unlock(lm_read_only);
    }
    else
    {
        void * x_cl(x.lock(lm_read_only, tags::OpenCL::GPU::memory_value));
        result = opencl::norm_l2_false<float>(x_cl, x.size(), CL_DEVICE_TYPE_GPU, "norm_l2_false_float");
        x.unlock(lm_read_only);
    }
    PROFILER_STOP("Norm l2 false tags::OpenCL::GPU");
    return result;
}

double Norm<vnt_l_two, false, tags::OpenCL::GPU>::value(const DenseVectorContinuousBase<double> & x)
{
    CONTEXT("When calculating L2 norm (false) of DenseVectorContinuousBase<double> (OpenCL GPU):");
    PROFILER_START("Norm l2 false tags::OpenCL::GPU");

    double result(0);
    if (x.size() < 32 * 128)
    {
        x.lock(lm_read_only);
        for (unsigned long i(0) ; i < x.size() ; ++i)
        {
            result += x[i] * x[i];
        }
        x.unlock(lm_read_only);
    }
    else
    {
        void * x_cl(x.lock(lm_read_only, tags::OpenCL::GPU::memory_value));
        result = opencl::norm_l2_false<double>(x_cl, x.size(), CL_DEVICE_TYPE_GPU, "norm_l2_false_double");
        x.unlock(lm_read_only);
    }
    PROFILER_STOP("Norm l2 false tags::OpenCL::GPU");
    return result;
}

float Norm<vnt_l_two, true, tags::OpenCL::GPU>::value(const DenseVectorContinuousBase<float> & x)
{
    CONTEXT("When calculating L2 norm (false) of DenseVectorContinuousBase<float> (OpenCL GPU):");
    PROFILER_START("Norm l2 true tags::OpenCL::GPU");

    float result(0);
    if (x.size() < 32 * 128)
    {
        x.lock(lm_read_only);
        for (unsigned long i(0) ; i < x.size() ; ++i)
        {
            result += x[i] * x[i];
        }
        x.unlock(lm_read_only);
    }
    else
    {
        void * x_cl(x.lock(lm_read_only, tags::OpenCL::GPU::memory_value));
        result = opencl::norm_l2_false<float>(x_cl, x.size(), CL_DEVICE_TYPE_GPU, "norm_l2_false_float");
        x.unlock(lm_read_only);
    }
    PROFILER_STOP("Norm l2 true tags::OpenCL::GPU");
    return sqrt(result);
}

double Norm<vnt_l_two, true, tags::OpenCL::GPU>::value(const DenseVectorContinuousBase<double> & x)
{
    CONTEXT("When calculating L2 norm (false) of DenseVectorContinuousBase<double> (OpenCL GPU):");
    PROFILER_START("Norm l2 true tags::OpenCL::GPU");

    double result(0);
    if (x.size() < 32 * 128)
    {
        x.lock(lm_read_only);
        for (unsigned long i(0) ; i < x.size() ; ++i)
        {
            result += x[i] * x[i];
        }
        x.unlock(lm_read_only);
    }
    else
    {
        void * x_cl(x.lock(lm_read_only, tags::OpenCL::GPU::memory_value));
        result = opencl::norm_l2_false<double>(x_cl, x.size(), CL_DEVICE_TYPE_GPU, "norm_l2_false_double");
        x.unlock(lm_read_only);
    }
    PROFILER_STOP("Norm l2 true tags::OpenCL::GPU");
    return sqrt(result);
}
