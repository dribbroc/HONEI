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

#include <honei/la/scale.hh>
#include <honei/backends/opencl/operations.hh>


using namespace honei;

DenseVectorContinuousBase<float> & Scale<tags::OpenCL::CPU>::value(DenseVectorContinuousBase<float> & x, const float a)
{
    CONTEXT("When scaling DenseVectorContinuousBase<float> by scalar (OpenCL CPU):");

    void * x_cl(x.lock(lm_read_and_write, tags::OpenCL::CPU::memory_value));
    opencl::scale_float(x_cl, a, x.size(), CL_DEVICE_TYPE_CPU);
    x.unlock(lm_read_and_write);

    return x;
}

DenseVectorContinuousBase<double> & Scale<tags::OpenCL::CPU>::value(DenseVectorContinuousBase<double> & x, const double a)
{
    CONTEXT("When scaling DenseVectorContinuousBase<double> by scalar (OpenCL CPU):");

    void * x_cl(x.lock(lm_read_and_write, tags::OpenCL::CPU::memory_value));
    opencl::scale_double(x_cl, a, x.size(), CL_DEVICE_TYPE_CPU);
    x.unlock(lm_read_and_write);

    return x;
}

DenseVectorContinuousBase<float> & Scale<tags::OpenCL::GPU>::value(DenseVectorContinuousBase<float> & x, const float a)
{
    CONTEXT("When scaling DenseVectorContinuousBase<float> by scalar (OpenCL GPU):");

    void * x_cl(x.lock(lm_read_and_write, tags::OpenCL::GPU::memory_value));
    opencl::scale_float(x_cl, a, x.size(), CL_DEVICE_TYPE_GPU);
    x.unlock(lm_read_and_write);

    return x;
}

DenseVectorContinuousBase<double> & Scale<tags::OpenCL::GPU>::value(DenseVectorContinuousBase<double> & x, const double a)
{
    CONTEXT("When scaling DenseVectorContinuousBase<double> by scalar (OpenCL GPU):");

    void * x_cl(x.lock(lm_read_and_write, tags::OpenCL::GPU::memory_value));
    opencl::scale_double(x_cl, a, x.size(), CL_DEVICE_TYPE_GPU);
    x.unlock(lm_read_and_write);

    return x;
}

