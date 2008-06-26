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

    void * b_gpu(b.read(tags::GPU::CUDA::memory_value));
    void * result_gpu(result.write(tags::GPU::CUDA::memory_value));

    for (BandedMatrix<float>::ConstBandIterator band(a.begin_non_zero_bands()), band_end(a.end_non_zero_bands()) ;
            band != band_end ; ++band)
    {
        // If we are above or on the diagonal band, we start at Element 0 and go on until Element band_size-band_index.
        if (band.index() >= middle_index)
        {
            op_offset = band.index() - middle_index;
            void * band_gpu(band->read(tags::GPU::CUDA::memory_value));
            void * temp_b_gpu = (void *)((float *)b_gpu + op_offset);
            cuda_scaled_sum_three_float(result_gpu, band_gpu, temp_b_gpu, a.size() - op_offset, blocksize);
            band->release_read();
        }
        else // If we are below the diagonal band, we start at Element 'start' and go on until the last element.
        {
            op_offset = middle_index - band.index();
            void * band_gpu(band->read(tags::GPU::CUDA::memory_value));
            band_gpu = (void *)((float *)band_gpu + op_offset);
            void * temp_result_gpu = (void *)((float *)result_gpu + op_offset);
            cuda_scaled_sum_three_float(temp_result_gpu, band_gpu, b_gpu,
                    a.size() - op_offset, blocksize);
            band->release_read();
        }
    }

    result.release_write();
    b.release_read();

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

    void * b_gpu(b.read(tags::GPU::CUDA::memory_value));
    void * result_gpu(result.write_only(tags::GPU::CUDA::memory_value));
    void * ll_gpu(a.band(LL).read(tags::GPU::CUDA::memory_value));
    void * ld_gpu(a.band(LD).read(tags::GPU::CUDA::memory_value));
    void * lu_gpu(a.band(LU).read(tags::GPU::CUDA::memory_value));
    void * dl_gpu(a.band(DL).read(tags::GPU::CUDA::memory_value));
    void * dd_gpu(a.band(DD).read(tags::GPU::CUDA::memory_value));
    void * du_gpu(a.band(DU).read(tags::GPU::CUDA::memory_value));
    void * ul_gpu(a.band(UL).read(tags::GPU::CUDA::memory_value));
    void * ud_gpu(a.band(UD).read(tags::GPU::CUDA::memory_value));
    void * uu_gpu(a.band(UU).read(tags::GPU::CUDA::memory_value));

    cuda_product_bmdv_q1_float(ll_gpu, ld_gpu, lu_gpu,
            dl_gpu, dd_gpu, du_gpu,
            ul_gpu, ud_gpu, uu_gpu,
            b_gpu, result_gpu, a.size(), blocksize, m);

    result.release_write();
    b.release_read();
    a.band(LL).release_read();
    a.band(LD).release_read();
    a.band(LU).release_read();
    a.band(DL).release_read();
    a.band(DD).release_read();
    a.band(DU).release_read();
    a.band(UL).release_read();
    a.band(UD).release_read();
    a.band(UU).release_read();

    return result;
}

