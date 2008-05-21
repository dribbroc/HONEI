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

    for (BandedMatrix<float>::ConstVectorIterator band(a.begin_non_zero_bands()), band_end(a.end_non_zero_bands()) ;
            band != band_end ; ++band)
    {
        // If we are above or on the diagonal band, we start at Element 0 and go on until Element band_size-band_index.
        if (band.index() >= middle_index)
        {
            op_offset = band.index() - middle_index;
            cuda_scaled_sum_three_float(result.elements(), band->elements(), b.elements() + op_offset, a.size() - op_offset, blocksize);
        }
        else // If we are below the diagonal band, we start at Element 'start' and go on until the last element.
        {
            op_offset = middle_index - band.index();
            cuda_scaled_sum_three_float(result.elements() + op_offset, band->elements() + op_offset, b.elements(), a.size() - op_offset, blocksize);
        }
    }
    return result;
}
