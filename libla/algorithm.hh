/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_ALGORITHM_HH
#define LIBLA_GUARD_ALGORITHM_HH 1

#include <libla/dense_vector.hh>
#include <libla/dense_matrix.hh>
#include <libla/sparse_matrix.hh>
#include <libla/banded_matrix.hh>
#include <libutil/tags.hh>

#include <limits>

namespace honei
{
    /**
     * \{
     *
     * Copies data from a source container to a destination container.
     *
     * \ingroup grpalgorithm
     */
    template <typename DT_> void copy(const DenseVector<DT_> & source, DenseVector<DT_> & dest)
    {
        CONTEXT("When copying elements from DenseVector to DenseVector:");
        ASSERT(source.elements() != dest.elements(),
                "trying to copy data from a DenseVector to the very same DenseVector!");

        if (source.size() != dest.size())
            throw VectorSizeDoesNotMatch(dest.size(), source.size());

        TypeTraits<DT_>::copy(source.elements(), dest.elements(), source.size());
    }
    /// \}

    /**
     * \{
     *
     * Converts data from a source container to a destination container.
     *
     * \ingroup grpalgorithm
     */
    template <typename DataType_, typename OrigType_>
    void convert(DenseVectorContinuousBase<DataType_> & copy, const DenseVectorContinuousBase<OrigType_> & orig)
    {
        CONTEXT("When converting DenseVector to DenseVector:");
        ASSERT(orig.size() > 0, "size is zero!");

        if (copy.size() != orig.size())
            throw VectorSizeDoesNotMatch(orig.size(), copy.size());

        /// \todo use typetraits convert
        typename Vector<DataType_>::ElementIterator f(copy.begin_elements());
        for (typename Vector<OrigType_>::ConstElementIterator i(orig.begin_elements()),
                i_end(orig.end_elements()) ; i != i_end ; ++i, ++f)
        {
            *f = *i;
        }
    }

    template <typename DataType_, typename OrigType_>
    void convert(SparseVector<DataType_> & copy, const SparseVector<OrigType_> & orig)
    {
        CONTEXT("When converting SparseVector to SparseVector:");
        ASSERT(orig.size() > 0, "size is zero!");

        if (copy.size() != orig.size())
            throw VectorSizeDoesNotMatch(orig.size(), copy.size());

        /// \todo use typetraits convert
        for (typename Vector<OrigType_>::ConstElementIterator i(orig.begin_non_zero_elements()),
                i_end(orig.end_non_zero_elements()) ; i != i_end ; ++i)
        {
            copy[i.index()] = *i;
        }
    }

    template <typename DataType_, typename OrigType_>
    void convert(DenseMatrix<DataType_> & copy, const Matrix<OrigType_> & orig)
    {
        CONTEXT("When converting DenseMatrix to DenseMatrix:");
        ASSERT(orig.size() > 0, "size is zero!");

        if (copy.rows() != orig.rows())
            throw MatrixRowsDoNotMatch(orig.rows(), copy.rows());
        if (copy.columns() != orig.columns())
            throw MatrixColumnsDoNotMatch(orig.columns(), copy.columns());

        /// \todo use typetraits convert
        typename MutableMatrix<DataType_>::ElementIterator f(copy.begin_elements());
        for (typename Matrix<OrigType_>::ConstElementIterator i(orig.begin_elements()),
                i_end(orig.end_elements()) ; i != i_end ; ++i, ++f)
        {
            *f = *i;
        }
    }

    template <typename DataType_, typename OrigType_>
    void convert(SparseMatrix<DataType_> & copy, const SparseMatrix<OrigType_> & orig)
    {
        CONTEXT("When converting SparseMatrix to SparseMatrix:");
        ASSERT(orig.size() > 0, "size is zero!");

        if (copy.rows() != orig.rows())
            throw MatrixRowsDoNotMatch(orig.rows(), copy.rows());
        if (copy.columns() != orig.columns())
            throw MatrixColumnsDoNotMatch(orig.columns(), copy.columns());

        for (unsigned long i(0) ; i < orig.rows() ; ++i)
        {
            convert(copy[i], orig[i]);
        }
    }

    template <typename DataType_, typename OrigType_>
    void convert(BandedMatrix<DataType_> & copy, const BandedMatrix<OrigType_> & orig)
    {
        CONTEXT("When converting BandedMatrix to BandedMatrix:");
        ASSERT(orig.size() > 0, "size is zero!");

        if (copy.size() != orig.size())
            throw MatrixSizeDoesNotMatch(orig.size(), copy.size());

        for (BandedMatrix<float>::ConstVectorIterator band(orig.begin_non_zero_bands()), band_end(orig.end_non_zero_bands()) ;
                band != band_end ; ++band)
        {
            convert(copy.band_unsigned(band.index()), *band);
        }
    }

    /// \}
}

#endif
