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

#include <honei/la/sum.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/backends/memory/memory_arbiter.hh>
#include <honei/util/configuration.hh>


using namespace honei;

DenseVectorContinuousBase<float> & Sum<tags::GPU::CUDA>::value(DenseVectorContinuousBase<float> & a,
        const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When adding DenseVectorContinuousBase<float> to DenseVectorContinuousBase<float> (CUDA):");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::sum_two_float", 128ul));

    void * a_gpu = MemoryArbiter::instance()->write<tags::GPU::CUDA>(a.memid(), a.address(), a.size() * sizeof(float));
    void * b_gpu = MemoryArbiter::instance()->read<tags::GPU::CUDA>(b.memid(), b.address(), b.size() * sizeof(float));
    cuda_sum_two_float(a_gpu, b_gpu, a.size(), blocksize);
    MemoryArbiter::instance()->release_read<tags::GPU::CUDA>(b.memid());
    MemoryArbiter::instance()->release_write<tags::GPU::CUDA>(a.memid());

    return a;
}

DenseMatrix<float> & Sum<tags::GPU::CUDA>::value(DenseMatrix<float> & a, const DenseMatrix<float> & b)
{
    CONTEXT("When adding DenseMatrix<float> to DenseMatrix<float> (CUDA):");

    if (a.columns() != b.columns())
    {
        throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
    }

    if (a.rows() != b.rows())
    {
        throw MatrixRowsDoNotMatch(b.rows(), a.rows());
    }

    unsigned long blocksize(Configuration::instance()->get_value("cuda::sum_two_float", 128ul));

    void * a_gpu = MemoryArbiter::instance()->write<tags::GPU::CUDA>(a.memid(), a.address(), a.size() * sizeof(float));
    void * b_gpu = MemoryArbiter::instance()->read<tags::GPU::CUDA>(b.memid(), b.address(), b.size() * sizeof(float));
    cuda_sum_two_float(a_gpu, b_gpu, a.size(), blocksize);
    MemoryArbiter::instance()->release_read<tags::GPU::CUDA>(b.memid());
    MemoryArbiter::instance()->release_write<tags::GPU::CUDA>(a.memid());

    return a;
}
