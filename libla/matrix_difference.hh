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

#ifndef LIBLA_GUARD_MATRIX_DIFFERENCE_HH
#define LIBLA_GUARD_MATRIX_DIFFERENCE_HH 1

#include <libutil/tags.hh>
#include <libla/dense_matrix.hh>
#include <libla/banded_matrix.hh>
#include <libla/sparse_matrix.hh>
#include <libla/matrix_error.hh>
#include <libla/scalar_vector_product.hh>
#include <libla/vector_difference.hh>

/**
 * \file
 *
 * Templatized definitions of operation MatrixDifference.
 *
 * \ingroup grpoperations
 */
namespace pg512
{
    /**
     * \brief Difference of two matrices.
     *
     * MatrixDifference is the class template for the subtraction operation
     * \f[
     *     MatrixDifference(a, b): \quad r \leftarrow a - b,
     * \f]
     * which yields r, the difference of matrices a and b.
     *
     * Usually, the return value is the minuend a after modification. However,
     * there are signatures for which b is returned. For these cases a short
     * notice is added.
     *
     * \ingroup grpoperations
     */

    template <typename Tag_ = tags::CPU> struct MatrixDifference
    {
        /**
         * \brief Returns the the difference of two given matrices.
         *
         * \param a The matrix that is the minuend of the operation.
         * \param b The matrix that is the subtrahend of the operation.
         *
         * \retval Will modify the minuend a and return it.
         */

        /// \{
        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
        {
            CONTEXT("When subtracting DenseMatrix from DenseMatrix:");

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
                *l -= *r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
        {
            CONTEXT("When subtracting SparseMatrix from DenseMatrix:");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            typename Matrix<DT2_>::ConstElementIterator r(b.begin_non_zero_elements()), r_end(b.end_non_zero_elements());
            for ( ; r != r_end ; ++r; )
            {
                a[r.row()][r.column()] -= *r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
        {
            CONTEXT("When subtracting SparseMatrix from SparseMatrix:");

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
                    l_end(a.end_non_zero_elements()) ; l != l_end ; )
            {
                if (r.index() < l.index())
                {
                    a[r.row()][r.column()] = -(*r);
                    ++r;
                }
                else if (l.index() < r.index())
                {
                    ++l;
                }
                else
                {
                    *l -= *r;
                    ++l; ++r;
                }
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static BandedMatrix<DT1_> & value(BandedMatrix<DT1_> & a,
                const BandedMatrix<DT2_> & b)
        {
            CONTEXT("When subtracting BandedMatrix from BandedMatrix:");

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
                    VectorDifference<>::value(*l, *r);
                else
                    a.band(r.index()) = ScalarVectorProduct<>::value(DT1_(-1), *((*r).copy()));
            }

            return a;
        }
        /// \}

        /**
         * \brief Returns the the difference of two given matrices.
         *
         * \param a The matrix that is the minuend of the operation.
         * \param b The matrix that is the subtrahend of the operation.
         *
         * \retval Will modify the subtrahend b and return it.
         */

        /// \{
        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT2_> & value(const BandedMatrix<DT1_> & a, DenseMatrix<DT2_> & b)
        {
            CONTEXT("When subtracting DenseMatrix from BandedMatrix:");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            /// \todo Optimise algorithm. Iterate over a's bands.
            typename Matrix<DT1_>::ConstElementIterator l(a.begin_elements());
            for (typename MutableMatrix<DT2_>::ElementIterator r(b.begin_elements()),
                    r_end(b.end_elements()) ; r != r_end ; ++r, ++l)
            {
                *r = *l - *r;
            }

            return b;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> & value(const BandedMatrix<DT2_> & a, SparseMatrix<DT1_> & b)
        {
            CONTEXT("When subtracting SparseMatrix from BandedMatrix:");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            /// \todo Optimise algorithm. Iterate over a's bands.
            typename Matrix<DT1_>::ConstElementIterator l(a.begin_elements()), l_end(a.begin_elements());
            for (typename MutableMatrix<DT2_>::ElementIterator r(b.begin_non_zero_elements()),
                    r_end(b.end_non_zero_elements()) ; r != r_end ; ++r, ++l)
            {
                while (l.index() < r.index())
                {
                    ++l;
                }

                *r = *l - *r;
                ++r;
            }

            return b;
        }
        /// \}
    };
}

#endif
