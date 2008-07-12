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

#include <honei/la/product.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/util/configuration.hh>

using namespace honei;

DenseVector<float> Product<tags::GPU::CUDA>::value(const BandedMatrix<float> & a, const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When multiplying BandedMatrix<float> with DenseVectorContinuousBase<float> (CUDA):");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }

    DenseVector<float> result(a.rows(), float(0));

    unsigned long blocksize(Configuration::instance()->get_value("cuda::scaled_sum_three_float", 128ul));

    unsigned long middle_index(a.rows() - 1);
    unsigned long op_offset;

    void * b_gpu(b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * result_gpu(result.lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

    for (BandedMatrix<float>::ConstBandIterator band(a.begin_non_zero_bands()), band_end(a.end_non_zero_bands()) ;
            band != band_end ; ++band)
    {
        // If we are above or on the diagonal band, we start at Element 0 and go on until Element band_size-band_index.
        if (band.index() >= middle_index)
        {
            op_offset = band.index() - middle_index;
            void * band_gpu(band->lock(lm_read_only, tags::GPU::CUDA::memory_value));
            void * temp_b_gpu = (void *)((float *)b_gpu + op_offset);
            cuda_scaled_sum_three_float(result_gpu, band_gpu, temp_b_gpu, a.size() - op_offset, blocksize);
            band->unlock(lm_read_only);
        }
        else // If we are below the diagonal band, we start at Element 'start' and go on until the last element.
        {
            op_offset = middle_index - band.index();
            void * band_gpu(band->lock(lm_read_only, tags::GPU::CUDA::memory_value));
            band_gpu = (void *)((float *)band_gpu + op_offset);
            void * temp_result_gpu = (void *)((float *)result_gpu + op_offset);
            cuda_scaled_sum_three_float(temp_result_gpu, band_gpu, b_gpu,
                    a.size() - op_offset, blocksize);
            band->unlock(lm_read_only);
        }
    }

    result.unlock(lm_read_and_write);
    b.unlock(lm_read_only);

    return result;
}


DenseVector<float> Product<tags::GPU::CUDA>::value(const BandedMatrixQ1<float> & a, const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When multiplying BandedMatrix<float> with DenseVectorContinuousBase<float> (CUDA):");

    if (b.size() != a.columns())
    {
        throw VectorSizeDoesNotMatch(b.size(), a.columns());
    }

    DenseVector<float> result(a.rows());

    unsigned long blocksize(Configuration::instance()->get_value("cuda::product_bmdv_q1_float", 128ul));
    unsigned long m(a.root());

    void * b_gpu(b.lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * result_gpu(result.lock(lm_write_only, tags::GPU::CUDA::memory_value));
    void * ll_gpu(a.band(LL).lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * ld_gpu(a.band(LD).lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * lu_gpu(a.band(LU).lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * dl_gpu(a.band(DL).lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * dd_gpu(a.band(DD).lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * du_gpu(a.band(DU).lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * ul_gpu(a.band(UL).lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * ud_gpu(a.band(UD).lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * uu_gpu(a.band(UU).lock(lm_read_only, tags::GPU::CUDA::memory_value));

    cuda_product_bmdv_q1_float(ll_gpu, ld_gpu, lu_gpu,
            dl_gpu, dd_gpu, du_gpu,
            ul_gpu, ud_gpu, uu_gpu,
            b_gpu, result_gpu, a.size(), blocksize, m);

    result.unlock(lm_write_only);
    b.unlock(lm_read_only);
    a.band(LL).unlock(lm_read_only);
    a.band(LD).unlock(lm_read_only);
    a.band(LU).unlock(lm_read_only);
    a.band(DL).unlock(lm_read_only);
    a.band(DD).unlock(lm_read_only);
    a.band(DU).unlock(lm_read_only);
    a.band(UL).unlock(lm_read_only);
    a.band(UD).unlock(lm_read_only);
    a.band(UU).unlock(lm_read_only);

    return result;
}

