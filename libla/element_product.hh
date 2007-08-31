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

#ifndef LIBLA_GUARD_ELEMENT_PRODUCT_HH
#define LIBLA_GUARD_ELEMENT_PRODUCT_HH 1

#include <libla/dense_matrix.hh>
#include <libla/dense_vector.hh>
#include <libla/banded_matrix.hh>
#include <libla/sparse_matrix.hh>
#include <libla/sparse_vector.hh>
#include <libla/vector_error.hh>
#include <libutil/tags.hh>

/**
 * \file
 *
 * Templatized definitions of elementwise vector products.
 *
 * \ingroup grpoperations
 */
namespace honei
{
    /**
     * \brief Multiplication of the elements of the given entities.
     *
     * ElementProduct is the template for the multiplication of the elements 
     * \f[
     *     ElementProduct(a,b): \quad a \leftarrow a[i]*b[i],
     * \f]
     *
     * of a given entity.
     *
     *
     * \ingroup grpoperations
     * \ingroup grpmatrixoperations
     * \ingroup grpvectoroperations
     */
    template <typename Tag_ = tags::CPU> struct ElementProduct
    {
        /**
         * \brief Returns the the result of elementwise multiplication of two entities.
         *
         * \param a Entity that is a factor of the operation.
         * \param b idem 
         *
         * \retval a Will modify the factor a and return it.
         *
         * \exception MatrixSizeDoesNotMatch is thrown if two banded matrices do not have the same size.
         * \exception MatrixRowsDoNotMatch is thrown if two matrices do not have the same number of rows.
         * \exception MatrixColumnsDoNotMatch is thrown if two matrices do not have the same number of columns.
         * \exception MatrixIsNotSquare is thrown if a row access matrix's number of rows does not equal its number of columns.
         */

        template <typename DT1_, typename DT2_>
        static DenseVector<DT1_> & value(DenseVector<DT1_> & left, const DenseVector<DT2_> & right)
        {
            CONTEXT("When calculating the product of DenseVectors elements");

            if (left.size() != right.size())
                throw VectorSizeDoesNotMatch(right.size(), left.size());

            typename Vector<DT2_>::ConstElementIterator r(right.begin_elements());
            for (typename Vector<DT1_>::ElementIterator l(left.begin_elements()),
                    l_end(left.end_elements()) ; l != l_end ; ++l)
            {
                *l *= *r;
                ++r;
            }

            return left;
        }

        template <typename DT1_, typename DT2_>
        static SparseVector<DT1_> value(SparseVector<DT1_> & left, const SparseVector<DT2_> & right)
        {
            CONTEXT("When calculating the product of SparseVector elements");

            if (left.size() != right.size())
                throw VectorSizeDoesNotMatch(right.size(), left.size());

            typename Vector<DT2_>::ConstElementIterator r(right.begin_non_zero_elements());
            for (typename Vector<DT1_>::ElementIterator l(left.begin_non_zero_elements()),
                    l_end(left.end_non_zero_elements()) ; l != l_end ;)
            {
                if (r.index() < l.index())
                {
                    ++r;
                }
                else if (l.index() < r.index())
                {
                    *l = DT1_(0);
                    ++l;
                }
                else
                {
                    *l *= *r;
                    ++l; ++r;
                }
            }
            return left;
            ///\todo: perhaps sparsify - in case l.index < r.index Write of 0 possible.
        }

        template <typename DT1_, typename DT2_>
        static SparseVector<DT1_> & value(SparseVector<DT1_> & left, const DenseVector<DT2_> & right)
        {
            CONTEXT("When calculating the product of SparseVector and DenseVector elements");

            if (left.size() != right.size())
                throw VectorSizeDoesNotMatch(right.size(), left.size());

            for (typename Vector<DT1_>::ElementIterator l(left.begin_non_zero_elements()),
                    l_end(left.end_non_zero_elements()) ; l != l_end ; ++l )
            {
                *l *= right[l.index()];
            }
            return left;
            ///\todo: perhaps sparsify - if *right[l.index()] == 0 -> write of zero.
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & left, const DenseMatrix<DT2_> & right)
        {
            CONTEXT("When calculating the product of DenseMatrix elements");

            if (left.columns() != right.columns())
            {
                throw MatrixColumnsDoNotMatch(right.columns(), left.columns());
            }

            if (left.rows() != right.rows())
            {
                throw MatrixRowsDoNotMatch(right.rows(), left.rows());
            }

            typename Matrix<DT2_>::ConstElementIterator r(right.begin_elements());
            for (typename MutableMatrix<DT1_>::ElementIterator l(left.begin_elements()),
                    l_end(left.end_elements()) ; l != l_end ; ++l)
            {
                *l *= *r;
                ++r;
            }

            return left;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & left, const DenseMatrix<DT2_> & right)
        {
            CONTEXT("When calculating the product of SparseMatrix and DenseMatrix elements");

            if (left.columns() != right.columns())
            {
                throw MatrixColumnsDoNotMatch(right.columns(), left.columns());
            }

            if (left.rows() != right.rows())
            {
                throw MatrixRowsDoNotMatch(right.rows(), left.rows());
            }

            for (typename MutableMatrix<DT1_>::ElementIterator l(left.begin_non_zero_elements()),
                    l_end(left.end_non_zero_elements()); l != l_end ; ++l)
            {
                *l *= right[l.row()][l.column()];
            }

            return left; ///\todo: perhaps sparsify, dense_matrix[row][col] may be zero.
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & left, const SparseMatrix<DT2_> & right)
        {
            CONTEXT("When calculating the product of SparseMatrix elements");

            if (left.columns() != right.columns())
            {
                throw MatrixColumnsDoNotMatch(right.columns(), left.columns());
            }

            if (left.rows() != right.rows())
            {
                throw MatrixRowsDoNotMatch(right.rows(), left.rows());
            }

            typename Matrix<DT2_>::ConstElementIterator r(right.begin_non_zero_elements());
            for (typename MutableMatrix<DT1_>::ElementIterator l(left.begin_non_zero_elements()),
                    l_end(left.end_non_zero_elements()) ; l != l_end ; )
            {
                if (l.index() < r.index())
                {
                    *l = DT1_(0);
                    ++l;
                }

                else if (r.index() < l.index())
                {
                    ++r;
                }

                else
                {
                    *l *= *r;
                    ++l; ++r;
                }
            }
            ///\todo: perhaps sparsify - in case l.index < r.index set to zero possible.
            return left;
        }

        template <typename DT1_, typename DT2_>
        static BandedMatrix<DT1_> & value(BandedMatrix<DT1_> & left, const BandedMatrix<DT2_> & right)
        {
            CONTEXT("When calculating the product of SparseVector elements");

            if (left.rows() != right.rows())
            {
                throw MatrixSizeDoesNotMatch(right.rows(), left.rows());
            }

            typename BandedMatrix<DT2_>::ConstVectorIterator r(right.begin_bands());
            for (typename BandedMatrix<DT1_>::VectorIterator l(left.begin_bands()),
                    l_end(left.end_bands()) ; l != l_end ; ++l)
            {
                if (! r.exists()) //won't check l.exists() here, cause it is always created by Iterator.
                {
                    ++r;
                    continue;
                }

                *l = ElementProduct<>::value(*l, *r);
                ++r;
            }

            return left;
        }

        template <typename DT1_, typename DT2_>
        static BandedMatrix<DT1_> & value(BandedMatrix<DT1_> & left, const DenseMatrix<DT2_> & right)
        {
            CONTEXT("When calculating the product of SparseMatrix and BandedMatrix elements");

            if (left.columns() != right.columns())
            {
                throw MatrixColumnsDoNotMatch(right.columns(), left.columns());
            }

            if (left.rows() != right.rows())
            {
                throw MatrixRowsDoNotMatch(right.rows(), left.rows());
            }

            typename Matrix<DT2_>::ConstElementIterator r(right.begin_elements());
            for (typename MutableMatrix<DT1_>::ElementIterator l(left.begin_elements()),
                    l_end(left.end_elements()) ; l != l_end ; ++l)
            {
                *l *= *r;
                ++r;
            }

            return left;
        }

        template <typename DT1_, typename DT2_>
        static BandedMatrix<DT1_> & value(BandedMatrix<DT1_> & left, const SparseMatrix<DT2_> & right)
        {
            CONTEXT("When calculating the product of BandedMatrix and a SparseMatrix elements");

            if (left.columns() != right.columns())
            {
                throw MatrixColumnsDoNotMatch(right.columns(), left.columns());
            }

            if (left.rows() != right.rows())
            {
                throw MatrixRowsDoNotMatch(right.rows(), left.rows());
            }

            typename MutableMatrix<DT1_>::ElementIterator l(left.begin_elements());
            for (typename Matrix<DT2_>::ConstElementIterator r(right.begin_non_zero_elements()),
                    r_end(right.end_non_zero_elements()); r != r_end ; )
            {
                while (l.index() < r.index())
                {
                    *l = DT1_(0);
                    ++l;
                }

                *l *= *r;
                ++l; ++r;
            }

            for (typename MutableMatrix<DT1_>::ElementIterator l_end(left.end_elements()) ;
                l != l_end ; ++l)
            {
                *l = DT1_(0);
            }
            /// \todo Left is complete dense at this point.
            return left;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & left, const BandedMatrix<DT2_> & right)
        {
            CONTEXT("When calculating the product of DenseMatrix and BandedMatrix elements");

            if (left.columns() != right.columns())
            {
                throw MatrixColumnsDoNotMatch(right.columns(), left.columns());
            }

            if (left.rows() != right.rows())
            {
                throw MatrixRowsDoNotMatch(right.rows(), left.rows());
            }

            typename Matrix<DT2_>::ConstElementIterator r(right.begin_elements());
            for (typename MutableMatrix<DT1_>::ElementIterator l(left.begin_elements()),
                    l_end(left.end_elements()) ; l != l_end ; ++l)
            {
                *l *= *r;
                ++r;
            }

            return left;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & left, const BandedMatrix<DT2_> & right)
        {
            CONTEXT("When calculating the product of SparseMatrix and BandedMatrix elements");

            if (left.columns() != right.columns())
            {
                throw MatrixColumnsDoNotMatch(right.columns(), left.columns());
            }

            if (left.rows() != right.rows())
            {
                throw MatrixRowsDoNotMatch(right.rows(), left.rows());
            }

            typename Matrix<DT2_>::ConstElementIterator r(right.begin_elements());
            for (typename MutableMatrix<DT1_>::ElementIterator l(left.begin_non_zero_elements()),
                    l_end(left.end_non_zero_elements()) ; l != l_end ; )
            {
                while (r.index() < l.index())
                    ++r;

                *l *= *r;
                ++r; ++l;
            }

            return left;
        }
    };
}
#endif
