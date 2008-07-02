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

#include <honei/la/element_inverse.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/util/configuration.hh>


using namespace honei;

DenseVectorContinuousBase<float> & ElementInverse<tags::GPU::CUDA>::value(DenseVectorContinuousBase<float> & x)
{
    CONTEXT("When inverting DenseVectorContinuousBase<float> (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::element_inverse_one_float", 128ul));

    void * x_gpu(x.write(tags::GPU::CUDA::memory_value));
    cuda_element_inverse_one_float(x_gpu, x.size(), blocksize);
    x.release_write();

    return x;
}

DenseMatrix<float> & ElementInverse<tags::GPU::CUDA>::value(DenseMatrix<float> & x)
{
    CONTEXT("When inverting DenseMatrix<float> (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::element_inverse_one_float", 128ul));

    void * x_gpu(x.write(tags::GPU::CUDA::memory_value));
    cuda_element_inverse_one_float(x_gpu, x.size(), blocksize);
    x.release_write();

    return x;
}
