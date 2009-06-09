/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c)  2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/la/dot_product.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/util/configuration.hh>

using namespace honei;

float DotProduct<tags::GPU::CUDA>::value(const DenseVectorContinuousBase<float> & a,
        const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When calculating dot-product of DenseVectorContinuousBase<float> with DenseVectorContinuousBase<float> "
            "(CUDA):");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::dot_product_two_float", 128ul));
    unsigned long gridsize(Configuration::instance()->get_value("cuda::dot_product_two_float_grid", 16ul));

    float result (0.);

    if (a.size() < gridsize * blocksize)
    {
        /// \todo run mini dot product in cuda
        a.lock(lm_read_only);
        b.lock(lm_read_only);
        for (unsigned long i(0) ; i < a.size() ; ++i)
        {
            result += a[i] * b[i];
        }
        a.unlock(lm_read_only);
        b.unlock(lm_read_only);
    }
    else
    {
        void * a_gpu(a.lock(lm_read_only, tags::GPU::CUDA::memory_value));
        void * b_gpu(b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
        result = cuda_dot_product_two_float(a_gpu, b_gpu, a.size(), blocksize, gridsize);
        b.unlock(lm_read_only);
        a.unlock(lm_read_only);
    }

    return result;
}

double DotProduct<tags::GPU::CUDA>::value(const DenseVectorContinuousBase<double> & a,
        const DenseVectorContinuousBase<double> & b)
{
    CONTEXT("When calculating dot-product of DenseVectorContinuousBase<double> with DenseVectorContinuousBase<double> "
            "(CUDA):");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::dot_product_two_double", 128ul));
    unsigned long gridsize(Configuration::instance()->get_value("cuda::dot_product_two_double", 16ul));

    double result (0.);

    if (a.size() < gridsize * blocksize)
    {
        /// \todo run mini dot product in cuda
        a.lock(lm_read_only);
        b.lock(lm_read_only);
        for (unsigned long i(0) ; i < a.size() ; ++i)
        {
            result += a[i] * b[i];
        }
        a.unlock(lm_read_only);
        b.unlock(lm_read_only);
    }
    else
    {
        void * a_gpu(a.lock(lm_read_only, tags::GPU::CUDA::memory_value));
        void * b_gpu(b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
        result = cuda_dot_product_two_double(a_gpu, b_gpu, a.size(), blocksize, gridsize);
        b.unlock(lm_read_only);
        a.unlock(lm_read_only);
    }

    return result;
}
