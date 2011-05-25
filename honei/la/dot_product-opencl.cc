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

#include <honei/la/dot_product.hh>
#include <honei/backends/opencl/operations.hh>


using namespace honei;

float DotProduct<tags::OpenCL::CPU>::value(const DenseVectorContinuousBase<float> & x, const DenseVectorContinuousBase<float> & y)
{
    CONTEXT("When calculating dot product of DenseVectorContinuousBase<float>(OpenCL CPU):");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(y.size(), x.size());

    float result(0);
    if (x.size() < 16 * 128)
    {
        x.lock(lm_read_only);
        y.lock(lm_read_only);
        for (unsigned long i(0) ; i < x.size() ; ++i)
        {
            result += x[i] * y[i];
        }
        x.unlock(lm_read_only);
        y.unlock(lm_read_only);
    }
    else
    {
        void * x_cl(x.lock(lm_read_only, tags::OpenCL::CPU::memory_value));
        void * y_cl(y.lock(lm_read_only, tags::OpenCL::CPU::memory_value));
        result = opencl::dot_product_float(x_cl, y_cl, x.size(), CL_DEVICE_TYPE_CPU);
        x.unlock(lm_read_only);
        y.unlock(lm_read_only);
    }
    return result;
}

double DotProduct<tags::OpenCL::CPU>::value(const DenseVectorContinuousBase<double> & x, const DenseVectorContinuousBase<double> & y)
{
    CONTEXT("When calculating dot product of DenseVectorContinuousBase<double>(OpenCL CPU):");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(y.size(), x.size());

    double result(0);
    if (x.size() < 16 * 128)
    {
        x.lock(lm_read_only);
        y.lock(lm_read_only);
        for (unsigned long i(0) ; i < x.size() ; ++i)
        {
            result += x[i] * y[i];
        }
        x.unlock(lm_read_only);
        y.unlock(lm_read_only);
    }
    else
    {
        void * x_cl(x.lock(lm_read_only, tags::OpenCL::CPU::memory_value));
        void * y_cl(y.lock(lm_read_only, tags::OpenCL::CPU::memory_value));
        result = opencl::dot_product_double(x_cl, y_cl, x.size(), CL_DEVICE_TYPE_CPU);
        x.unlock(lm_read_only);
        y.unlock(lm_read_only);
    }
    return result;
}

float DotProduct<tags::OpenCL::GPU>::value(const DenseVectorContinuousBase<float> & x, const DenseVectorContinuousBase<float> & y)
{
    CONTEXT("When calculating dot product of DenseVectorContinuousBase<float>(OpenCL GPU):");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(y.size(), x.size());

    float result(0);
    if (x.size() < 16 * 128)
    {
        x.lock(lm_read_only);
        y.lock(lm_read_only);
        for (unsigned long i(0) ; i < x.size() ; ++i)
        {
            result += x[i] * y[i];
        }
        x.unlock(lm_read_only);
        y.unlock(lm_read_only);
    }
    else
    {
        void * x_cl(x.lock(lm_read_only, tags::OpenCL::GPU::memory_value));
        void * y_cl(y.lock(lm_read_only, tags::OpenCL::GPU::memory_value));
        result = opencl::dot_product_float(x_cl, y_cl, x.size(), CL_DEVICE_TYPE_GPU);
        x.unlock(lm_read_only);
        y.unlock(lm_read_only);
    }
    return result;
}

double DotProduct<tags::OpenCL::GPU>::value(const DenseVectorContinuousBase<double> & x, const DenseVectorContinuousBase<double> & y)
{
    CONTEXT("When calculating dot product of DenseVectorContinuousBase<double>(OpenCL GPU):");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(y.size(), x.size());

    double result(0);
    if (x.size() < 16 * 128)
    {
        x.lock(lm_read_only);
        y.lock(lm_read_only);
        for (unsigned long i(0) ; i < x.size() ; ++i)
        {
            result += x[i] * y[i];
        }
        x.unlock(lm_read_only);
        y.unlock(lm_read_only);
    }
    else
    {
        void * x_cl(x.lock(lm_read_only, tags::OpenCL::GPU::memory_value));
        void * y_cl(y.lock(lm_read_only, tags::OpenCL::GPU::memory_value));
        result = opencl::dot_product_double(x_cl, y_cl, x.size(), CL_DEVICE_TYPE_GPU);
        x.unlock(lm_read_only);
        y.unlock(lm_read_only);
    }
    return result;
}