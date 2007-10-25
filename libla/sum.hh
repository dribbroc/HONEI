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

#ifndef LIBLA_GUARD_SUM_HH
#define LIBLA_GUARD_SUM_HH 1

#include <libla/banded_matrix.hh>
#include <libla/dense_matrix.hh>
#include <libla/matrix_error.hh>
#include <libla/sparse_matrix.hh>
#include <libutil/tags.hh>

#include <algorithm>

namespace honei
{
    template <typename Tag_ = tags::CPU> struct Sum;

    /**
     * \brief Sum of two entities
     *
     * Sum is the class template for the addition operation
     * \f[
     *     \texttt{Sum}(a, b): \quad r \leftarrow a + b,
     * \f]
     * which yields r, the sum of entities a and b.
     *
     * The return value is the summand a after modification.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <> struct Sum<tags::CPU>
    {
        /**
         * \name Sums of two matrices
         * \{
         *
         * \brief Returns the the sum of two given entities.
         *
         * \param a The entity that is the left-hand summand of the operation.
         * \param b The entity that is the right-hand summand of the operation.
         *
         * \retval a Will modify the summand a and return it.
         *
         * \exception MatrixSizeDoesNotMatch is thrown if two banded matrices do not have the same size.
         * \exception MatrixRowsDoNotMatch is thrown if two matrices do not have the same number of rows.
         * \exception MatrixColumnsDoNotMatch is thrown if two matrices do not have the same number of columns.
         * \exception MatrixIsNotSquare is thrown if a row access matrix's number of rows does not equal its number of columns.
         * \exception VectorSizeDoesNotMatch is thrown if the two vectors don't have the same size.
         */

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

            typename MutableMatrix<DT1_>::ElementIterator l(a.begin_non_zero_elements());
            for (typename Matrix<DT2_>::ConstElementIterator r(b.begin_non_zero_elements()),
                    r_end(b.end_non_zero_elements()) ; r != r_end ; )
            {
                if (r.index() < l.index())
                {
                    a[r.row()][r.column()] = (*r);
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
                    Sum<>::value(*l, *r);
                else
                    a.band(r.index()) = r->copy();
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
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
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

            if (a.columns() != a.rows())
            {
                throw MatrixIsNotSquare(a.rows(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            typename MutableMatrix<DT1_>::ElementIterator l(a.begin_non_zero_elements());
            for (typename Matrix<DT2_>::ConstElementIterator r(b.begin_elements()), r_end(b.end_elements()) ;
                    r != r_end ; ++r)
            {
                if (r.index() < l.index())
                {
                    a[r.row()][r.column()] = *r;
                }
                else // l.index() < r.index() not possible
                {
                *l += *r;
                ++l;
                }
            }

            return a;
        }

        /// \}

        /**
         * \name Sums of scalar and matrix
         * \{
         *
         * Returns the matrix after adding a scalar to every element.
         *
         * \param a The DenseMatrix to be used.
         * \param b The scalar to be added.
         *
         * \retval a The referenced matrix is changed by adding the given scalar to each of its elements.
         */

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const DT2_ b)
        {
            CONTEXT("When adding scalar to DenseMatrix:");

            for (typename MutableMatrix<DT1_>::ElementIterator l(a.begin_elements()),
                    l_end(a.end_elements()) ; l != l_end ; ++l)
            {
                *l += b;
            }

            return a;
        }

        /// \}

        /**
         * \name Sums of two vectors
         * \{
         *
         * \brief Returns the the sum of two given vectors.
         *
         * \param a The vector that is the left-hand summand of the operation.
         * \param b The vector that is the right-hand summand of the operation.
         *
         * \retval Will modify the summand a and return it.
         *
         * \exception VectorSizeDoesNotMatch is thrown if the two vectors don't have the same size.
         */

        template <typename DT1_, typename DT2_>
        static DenseVector<DT1_> & value(DenseVector<DT1_> & a, const DenseVector<DT2_> & b)
        {
            CONTEXT("When adding DenseVector to DenseVector:");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());

            typename Vector<DT2_>::ConstElementIterator r(b.begin_elements());
            for (typename Vector<DT1_>::ElementIterator l(a.begin_elements()),
                    l_end(a.end_elements()) ; l != l_end ; ++l)
            {
                *l += *r;
                ++r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static SparseVector<DT1_> & value(SparseVector<DT1_> & a, const SparseVector<DT2_> & b)
        {
            CONTEXT("When adding SparseVector to SparseVector:");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());

            typename Vector<DT1_>::ElementIterator l(a.begin_non_zero_elements());
            for (typename Vector<DT2_>::ConstElementIterator r(b.begin_non_zero_elements()),
                    r_end(b.end_non_zero_elements()) ; r != r_end ; )
            {
                if (r.index() < l.index())
                {
                    a[r.index()] = *r;
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
            ///\todo: perhaps sparsify - i.e. addition of -7 and 7 possible
            return a;
        }

        template <typename DT1_, typename DT2_>
        static DenseVector<DT1_> & value(DenseVector<DT1_> & a, const SparseVector<DT2_> & b)
        {
            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());

            for (typename Vector<DT2_>::ConstElementIterator r(b.begin_non_zero_elements()),
                    r_end(b.end_non_zero_elements()) ; r != r_end ; ++r )
            {
                a[r.index()] += *r;
            }

            return a;
        }

        /// \}

        /**
         * \name Sums of scalar and vector
         * \{
         *
         * Returns the vector after adding a scalar to every element.
         *
         * \param x The Vector that shall be added.
         * \param a The scalar that shall be added.
         *
         * \retval x Will return x after modification.
         */

        template <typename DT_>
        static DenseVector<DT_> & value(const DT_ a, DenseVector<DT_> & x)
        {
            CONTEXT("When adding scalar to DenseVector:");

            for (typename Vector<DT_>::ElementIterator l(x.begin_elements()),
                    l_end(x.end_elements()) ; l != l_end ; ++l)
            {
                *l += a;
            }

            return x;
        }

        /// \}
    };

    /**
     * \brief Sum of two entities
     *
     * Sum is the class template for the addition operation
     * \f[
     *     \texttt{Sum}(a, b): \quad r \leftarrow a + b,
     * \f]
     * which yields r, the sum of entities a and b.
     *
     * The return value is the summand a after modification.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <> struct Sum<tags::Cell>
    {
        /**
         * \name Sums of two vectors
         * \{
         *
         * \brief Returns the the sum of two given vectors.
         *
         * \param a The vector that is the left-hand summand of the operation.
         * \param b The vector that is the right-hand summand of the operation.
         *
         * \retval Will modify the summand a and return it.
         *
         * \exception VectorSizeDoesNotMatch is thrown if the two vectors don't have the same size.
         */

        static DenseVector<float> & value(DenseVector<float> & a, const DenseVector<float> & b);
        static DenseVector<float> & value(DenseVector<float> & a, const SparseVector<float> & b);

        /// \}
    };

    /**
     * \brief Sum of two entities
     *
     * Sum is the class template for the addition operation
     * \f[
     *     \texttt{Sum}(a, b): \quad r \leftarrow a + b,
     * \f]
     * which yields r, the sum of entities a and b.
     *
     * The return value is the summand a after modification.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <> struct Sum<tags::CPU::SSE>
    {
        /**
         * \name Sums of two vectors
         * \{
         *
         * \brief Returns the the sum of two given vectors.
         *
         * \param a The vector that is the left-hand summand of the operation.
         * \param b The vector that is the right-hand summand of the operation.
         *
         * \retval Will modify the summand a and return it.
         *
         * \exception VectorSizeDoesNotMatch is thrown if the two vectors don't have the same size.
         */

        static DenseVector<float> & value(DenseVector<float> & a, const DenseVector<float> & b);

        static DenseVector<double> & value(DenseVector<double> & a, const DenseVector<double> & b);

        /// \}
    };
}

#endif
