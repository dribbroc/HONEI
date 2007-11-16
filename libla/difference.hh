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

#ifndef LIBLA_GUARD_DIFFERENCE_HH
#define LIBLA_GUARD_DIFFERENCE_HH 1

#include <libla/banded_matrix.hh>
#include <libla/dense_matrix.hh>
#include <libla/difference.hh>
#include <libla/matrix_error.hh>
#include <libla/product.hh>
#include <libla/scale.hh>
#include <libla/sparse_matrix.hh>
#include <libutil/tags.hh>

namespace honei
{
    template <typename Tag_ = tags::CPU> struct Difference;

    /**
     * \brief Difference of two entities.
     *
     * Difference is the class template for the subtraction operation
     * \f[
     *     \texttt{Difference}(a, b): \quad r \leftarrow a - b,
     * \f]
     * which yields r, the difference of entities a and b.
     *
     * Usually, the return value is the minuend a after modification. However,
     * there are signatures for which b is returned. For these cases a short
     * notice is added.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <> struct Difference<tags::CPU>
    {
        /**
         * \name Matrix differences that return minuend
         * \{
         *
         * \brief Returns the the difference of two given matrices
         *
         * \param a The entity that is the minuend of the operation.
         * \param b The entity that is the subtrahend of the operation.
         *
         * \retval r Will modify the minuend a and return it.
         *
         * \exception MatrixSizeDoesNotMatch is thrown if two banded matrices do not have the same size.
         * \exception MatrixRowsDoNotMatch is thrown if two matrices do not have the same number of rows.
         * \exception MatrixColumnsDoNotMatch is thrown if two matrices do not have the same number of columns.
         * \exception MatrixIsNotSquare is thrown if a RowAccessMatrix's number of rows does not equal its number of columns.
         */

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

            typename Matrix<DT2_>::ConstElementIterator r(b.begin_non_zero_elements()),
                     r_end(b.end_non_zero_elements());
            for ( ; r != r_end ; ++r )
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

            typename MutableMatrix<DT1_>::ElementIterator l(a.begin_non_zero_elements());
            for (typename Matrix<DT2_>::ConstElementIterator r(b.begin_non_zero_elements()),
                    r_end(b.end_non_zero_elements()) ; r != r_end ; )
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
        static BandedMatrix<DT1_> & value(BandedMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
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
                    Difference<>::value(*l, *r);
                else
                {
                    DenseVector<DT2_> band(r->copy());
                    Scale<>::value(DT1_(-1), band);
                    a.band(r.index()) = band;
                }
            }

            return a;
        }

        /// \}

        /**
         * \name Matrix differences that return subtrahend
         * \{
         *
         * \brief Returns the the difference of two given matrices.
         *
         * \param a The matrix that is the minuend of the operation.
         * \param b The matrix that is the subtrahend of the operation.
         *
         * \retval b Will modify the subtrahend b and return it.
         *
         * \exception MatrixSizeDoesNotMatch is thrown if two banded matrices do not have the same size.
         * \exception MatrixRowsDoNotMatch is thrown if two matrices do not have the same number of rows.
         * \exception MatrixColumnsDoNotMatch is thrown if two matrices do not have the same number of columns.
         * \exception MatrixIsNotSquare is thrown if a row access matrix's number of rows does not equal its number of columns.
         */

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT2_> & value(const BandedMatrix<DT1_> & a, DenseMatrix<DT2_> & b)
        {
            CONTEXT("When subtracting DenseMatrix from BandedMatrix:");

            if (b.columns() != b.rows())
            {
                throw MatrixIsNotSquare(b.rows(), b.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

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

            if (b.columns() != b.rows())
            {
                throw MatrixIsNotSquare(b.rows(), b.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            for (typename Matrix<DT1_>::ConstElementIterator l(a.begin_elements()),
                    l_end(a.end_elements()) ; l != l_end ; ++l)
            {
                    b[l.row()][l.column()] = (b[l.row()][l.column()] * (-1)) + *l;
            }

            return b;
        }

        /// \}

        /**
         * \name Vector differences
         * \{
         *
         * \brief Returns the the difference of two given vectors.
         *
         * \param a The vector that is the left-hand minuend of the operation.
         * \param b The vector that is the right-hand subtrahend of the operation.
         *
         * \retval a Will normally modify the minuend a and return it. Only for
         *           Difference(SparseVector, DenseVector) the subtrahend b is modified and returned.
         *
         * \exception VectorSizeDoesNotMatch is thrown if the two vectors don't have the same size.
         */

        template <typename DT1_, typename DT2_>
        static DenseVectorBase<DT1_> & value(DenseVectorBase<DT1_> & a, const DenseVectorBase<DT2_> & b)
        {
            CONTEXT("When subtracting DenseVectorBase from DenseVectorBase:");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());

            typename Vector<DT2_>::ConstElementIterator r(b.begin_elements());
            for (typename Vector<DT1_>::ElementIterator l(a.begin_elements()),
                    l_end(a.end_elements()) ; l != l_end ; ++l)
            {
                *l -= *r;
                ++r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static SparseVector<DT1_> & value(SparseVector<DT1_> & a, const SparseVector<DT2_> & b)
        {
            CONTEXT("When subtracting SparseVector from SparseVector:");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());

            typename Vector<DT1_>::ElementIterator l(a.begin_non_zero_elements());
            for (typename Vector<DT2_>::ConstElementIterator r(b.begin_non_zero_elements()),
                    r_end(b.end_non_zero_elements()) ; r != r_end ; )
            {
                if (r.index() < l.index())
                {
                    a[r.index()] = -(*r);
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
            ///\todo: perhaps sparsify - i.e. substraction of 7 and 7 possible.
        }

        template <typename DT1_, typename DT2_>
        static DenseVectorBase<DT1_> & value(DenseVectorBase<DT1_> & a, const SparseVector<DT2_> & b)
        {
            CONTEXT("When subtracting SparseVector from DenseVectorBase:");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());

            for (typename Vector<DT2_>::ConstElementIterator r(b.begin_non_zero_elements()),
                    r_end(b.end_non_zero_elements()) ; r != r_end ; ++r)
            {
                a[r.index()] -= *r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static DenseVectorBase<DT1_> & value(SparseVector<DT1_> & a, const DenseVectorBase<DT2_> & b)
        {
            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());

            typename Vector<DT1_>::ConstElementIterator l(a.begin_elements());
            for (typename Vector<DT2_>::ElementIterator r(b.begin_elements()),
                    r_end(b.end_elements()) ; r != r_end ; ++r)
            {
                *r = (-1) * *r + *l;
            }

            return b;
        }

        /// \}
    };

    /**
     * \brief Difference of two entities.
     *
     * Difference is the class template for the subtraction operation
     * \f[
     *     \texttt{Difference}(a, b): \quad r \leftarrow a - b,
     * \f]
     * which yields r, the difference of entities a and b.
     *
     * Usually, the return value is the minuend a after modification. However,
     * there are signatures for which b is returned. For these cases a short
     * notice is added.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <>
    struct Difference<tags::CPU::SSE>
    {
        /**
         * \name Vector differences
         * \{
         *
         * \brief Returns the the difference of two given vectors.
         *
         * \param a The vector that is the left-hand minuend of the operation.
         * \param b The vector that is the right-hand subtrahend of the operation.
         *
         * \retval a Will normally modify the minuend a and return it. Only for
         *           Difference(SparseVector, DenseVector) the subtrahend b is modified and returned.
         *
         * \exception VectorSizeDoesNotMatch is thrown if the two vectors don't have the same size.
         */

        static DenseVector<float> & value(DenseVector<float> & a, const DenseVector<float> & b);

        static DenseVector<double> & value(DenseVector<double> & a, const DenseVector<double> & b);

        /// \}
    };
}
#endif
