/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Sven Mallach <sven.mallach@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_MATRIX_SUM_HH
#define LIBLA_GUARD_MATRIX_SUM_HH 1

#include <libutil/tags.hh>
#include <libla/dense_matrix.hh>
#include <libla/sparse_matrix.hh>
#include <libla/banded_matrix.hh>
#include <libla/matrix_error.hh>
#include <libla/vector_sum.hh>

#include <algorithm>

/**
 * \file
 *
 * Templatized definitions of operation MatrixSum.
 *
 * \ingroup grpoperations
 **/
namespace pg512
{
    /**
     * \brief Sum of two matrices.
     *
     * MatrixSum is the class template for the addition operation
     * \f[
     *     MatrixSum(a, b): \quad r \leftarrow a + b,
     * \f]
     * which yields r, the sum of matrices a and b.
     *
     * The return value is the summand a after modification.
     *
     * \ingroup grpoperations
     */

    template <typename Tag_ = tags::CPU> struct MatrixSum
    {
        /**
         * \brief Returns the the sum of two given matrices.
         *
         * \param a The matrix that is the left-hand summand of the operation.
         * \param b The matrix that is the right-hand summand of the operation.
         *
         * \retval Will modify the summand a and return it.
         */

        /// \{
        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
        {
            CONTEXT("When adding DenseMatrix to DenseMatrix:");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            typename Matrix<DT2_>::ConstElementIterator r(b.begin_elements());
            for (typename MutableMatrix<DT1_>::ElementIterator l(a.begin_elements()),
                    l_end(a.end_elements()) ; l != l_end ; ++l, ++r)
            {
                *l += *r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
        {
            CONTEXT("When adding SparseMatrix to DenseMatrix:");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            // DenseMatrix<>::operator[] is as fast as the iterator-based access.
            // Avoid using unnecessary branches by iterating over a's elements.
            for (typename Matrix<DT2_>::ConstElementIterator r(b.begin_non_zero_elements()),
                    r_end(b.end_non_zero_elements()) ; r != r_end ; ++r)
            {
                a[r.row()][r.column()] += *r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
        {
            CONTEXT("When adding SparseMatrix to SparseMatrix:");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            typename Matrix<DT2_>::ConstElementIterator r(b.begin_non_zero_elements());
            for (typename MutableMatrix<DT1_>::ElementIterator l(a.begin_non_zero_elements()),
                    l_end(a.end_elements()) ; l != l_end ; )
            {
                if (r.index() < l.index())
                {
                    a[r.row()][r.column()] = *r;
                    ++r;
                }
                else if (l.index() < r.index())
                {
                    ++l;
                }
                else
                {
                    *l += *r;
                    ++l; ++r;
                }
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static BandedMatrix<DT1_> & value(BandedMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
        {
            CONTEXT("When adding BandedMatrix to BandedMatrix:");

            if (a.size() != b.size())
            {
                throw MatrixSizeDoesNotMatch(b.size(), a.size());
            }

            typename BandedMatrix<DT1_>::VectorIterator l(a.begin_bands()), l_end(a.end_bands());
            typename BandedMatrix<DT2_>::ConstVectorIterator r(b.begin_bands()), r_end(b.end_bands());
            for ( ; ((l != l_end) && (r != r_end)) ; ++l, ++r)
            {
                if (! r.exists())
                    continue;

                if (l.exists())
                    VectorSum<>::value(*l, *r);
                else
                    a.band(r.index()) = *((*r).copy());
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
        {
            CONTEXT("When adding BandedMatrix to DenseMatrix:");

            if (a.columns() != a.rows())
            {
                throw MatrixIsNotSquare(a.rows(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixSizeDoesNotMatch(b.rows(), a.rows());
            }

            for (typename BandedMatrix<DT2_>::ConstVectorIterator r(b.begin_bands()), r_end(b.end_bands()) ;
                    r != r_end ; ++r)
            {
                if (! r.exists())
                    continue;

                unsigned long size(b.size());
                unsigned long row_index(-std::max(long(r.index() - size + 1), long(0)));
                unsigned long col_index(std::min(long(r.index() - size + 1), long(0)));

                for (typename Vector<DT2_>::ConstElementIterator c(r->begin_elements()), c_end(r->end_elements()) ;
                        c != c_end ; ++c)
                {
                    if (row_index + c.index() >= size)
                        break;

                    if (col_index + c.index() >= size)
                        break;

                    a[row_index + c.index()][col_index + c.index()] += *c;
                }
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
        {
            CONTEXT("When adding BandedMatrix to SparseMatrix:");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            for (typename MutableMatrix<DT1_>::ElementIterator l(a.begin_non_zero_elements()),
                    l_end(a.end_non_zero_elements()) ; l != l_end ; )
            {
                *l += b[l.row()][l.column()];
                ++l;
            }

            return a;
        }
        /// \}
    };
}
#endif
