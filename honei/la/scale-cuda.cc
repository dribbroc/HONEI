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

#include <honei/la/scale.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/util/configuration.hh>


using namespace honei;

DenseVectorContinuousBase<float> & Scale<tags::GPU::CUDA>::value(DenseVectorContinuousBase<float> & x, const float a)
{
    CONTEXT("When scaling DenseVectorContinuousBase<float> by float (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::scale_one_float", 128ul));

    void * x_gpu (x.lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
    cuda_scale_one_float(x_gpu, a, x.size(), blocksize);
    x.unlock(lm_read_and_write);

    return x;
}

DenseMatrix<float> & Scale<tags::GPU::CUDA>::value(DenseMatrix<float> & x, const float a)
{
    CONTEXT("When scaling DenseMatrix<float> by float (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::scale_one_float", 128ul));

    void * x_gpu (x.lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
    cuda_scale_one_float(x_gpu, a, x.size(), blocksize);
    x.unlock(lm_read_and_write);

    return x;
}
