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

#ifndef LIBLA_GUARD_MATRIX_ELEMENTWISE_PRODUCT_HH
#define LIBLA_GUARD_MATRIX_ELEMENTWISE_PRODUCT_HH 1

#include <libutil/tags.hh>
#include <libla/dense_matrix.hh>
#include <libla/banded_matrix.hh>
#include <libla/sparse_matrix.hh>
#include <libla/matrix_error.hh>
#include <libla/vector_elementwise_product.hh>

/**
 * \file
 *
 * Templatized definitions of elementwise matrix products.
 *
 * \ingroup grpmatrixoperations
 **/
namespace pg512
{

    /**
     * MatrixElementwiseProduct is the class template for multiplying two matrices elementwise.
     * \brief The first referenced matrix is changed under this operation.
     * \ingroup grpmatrixoperations
     **/
    template <typename Tag_ = tags::CPU> struct MatrixElementwiseProduct
    {
        /**
         * Returns the resulting matrix after multiplying two DenseMatrix instances elementwise.
         * \param left Reference to a DenseMatrix. Its return type is used for the result matrix.
         * \param right Reference to a DenseMatrix.
         **/
        template <typename DataType1_, typename DataType2_> static DenseMatrix<DataType1_> & value(DenseMatrix<DataType1_> & left, const DenseMatrix<DataType2_> & right)
        {
            if (left.columns() != right.columns())
            {
                throw MatrixColumnsDoNotMatch(right.columns(), left.columns());
            }

            if (left.rows() != right.rows())
            {
                throw MatrixRowsDoNotMatch(right.rows(), left.rows());
            }

            typename Matrix<DataType2_>::ConstElementIterator r(right.begin_elements());
            for (typename MutableMatrix<DataType1_>::ElementIterator l(left.begin_elements()),
                    l_end(left.end_elements()) ; l != l_end ; ++l)
            {
                *l *= *r;
                ++r;
            }

            return left;
        }

		/**
         * Returns the resulting matrix after multiplying a sparse and a dense matrix instance elementwise.
         * \param left Reference to a SparseMatrix. Its return type is used for the result matrix.
         * \param right Reference to a DenseMatrix.
         **/
        template <typename DataType1_, typename DataType2_> static SparseMatrix<DataType1_> & value(SparseMatrix<DataType1_> & left, const DenseMatrix<DataType2_> & right)
        {
            if (left.columns() != right.columns())
            {
                throw MatrixColumnsDoNotMatch(right.columns(), left.columns());
            }

            if (left.rows() != right.rows())
            {
                throw MatrixRowsDoNotMatch(right.rows(), left.rows());
            }

            for (typename MutableMatrix<DataType1_>::ElementIterator l(left.begin_non_zero_elements()),
                    l_end(left.end_non_zero_elements()); l != l_end ; ++l)
            {
                *l *= right[l.row()][l.column()];
            }

            return left; ///\todo: perhaps sparsify, dense_matrix[row][col] may be zero.
        }

		/**
         * Returns the resulting matrix after multiplying two sparse matrix instances elementwise.
         * \param left Reference to a SparseMatrix. Its return type is used for the result matrix.
         * \param right Reference to a SparseMatrix
         **/
        template <typename DataType1_, typename DataType2_> static SparseMatrix<DataType1_> & value(SparseMatrix<DataType1_> & left, const SparseMatrix<DataType2_> & right)
        {
            if (left.columns() != right.columns())
            {
                throw MatrixColumnsDoNotMatch(right.columns(), left.columns());
            }

            if (left.rows() != right.rows())
            {
                throw MatrixRowsDoNotMatch(right.rows(), left.rows());
            }

            typename Matrix<DataType2_>::ConstElementIterator r(right.begin_non_zero_elements());
            for (typename MutableMatrix<DataType1_>::ElementIterator l(left.begin_non_zero_elements()),
                    l_end(left.end_non_zero_elements()) ; l != l_end ; )
            {
				if (l.index() < r.index())
				{
					*l = DataType1_(0);
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

        /**
         * Returns the resulting matrix after multiplying two Banded Matrix instances elementwise.
         * \param left Reference to a Banded Matrix. Its return type is used for the result matrix.
         * \param right Reference to a Banded Matrix.
         **/
        template <typename DataType1_, typename DataType2_> static BandedMatrix<DataType1_> & value(BandedMatrix<DataType1_> & left, const BandedMatrix<DataType2_> & right)
        {
            if (left.rows() != right.rows())
            {
                throw MatrixSizeDoesNotMatch(right.rows(), left.rows());
            }

            typename BandedMatrix<DataType2_>::ConstVectorIterator r(right.begin_bands());
            for (typename BandedMatrix<DataType1_>::VectorIterator l(left.begin_bands()),
                    l_end(left.end_bands()) ; l != l_end ; ++l)
            {
                if (! r.exists()) //won't check l.exists() here, cause it is always created by Iterator.
                {
                    ++r;
                    continue;
                }

                *l = VectorElementwiseProduct<>::value(*l, *r);
                ++r;
            }

            return left;
        }

        /**
         * Returns the resulting matrix after multiplying a DenseMatrix instance to a BandedMatrix instance elementwise.
         * \param left Reference to a BandedMatrix. Its return type is used for the result matrix.
         * \param right Reference to a DenseMatrix.
         **/
        template <typename DataType1_, typename DataType2_> static BandedMatrix<DataType1_> & value(BandedMatrix<DataType1_> & left, const DenseMatrix<DataType2_> & right)
        {
            if (left.columns() != right.columns())
            {
                throw MatrixColumnsDoNotMatch(right.columns(), left.columns());
            }

            if (left.rows() != right.rows())
            {
                throw MatrixRowsDoNotMatch(right.rows(), left.rows());
            }

            typename Matrix<DataType2_>::ConstElementIterator r(right.begin_elements());
            for (typename MutableMatrix<DataType1_>::ElementIterator l(left.begin_elements()),
                    l_end(left.end_elements()) ; l != l_end ; ++l)
            {
                *l *= *r;
                ++r;
            }

            return left;
        }

        /**
         * Returns the resulting matrix after multiplying a SparseMatrix instance to a BandedMatrix instance elementwise.
         * \param left Reference to a BandedMatrix. Its return type is used for the result matrix.
         * \param right Reference to a SparseMatrix.
         **/
        template <typename DataType1_, typename DataType2_> static BandedMatrix<DataType1_> & value(BandedMatrix<DataType1_> & left, const SparseMatrix<DataType2_> & right)
        {
            if (left.columns() != right.columns())
            {
                throw MatrixColumnsDoNotMatch(right.columns(), left.columns());
            }

            if (left.rows() != right.rows())
            {
                throw MatrixRowsDoNotMatch(right.rows(), left.rows());
            }

            typename MutableMatrix<DataType1_>::ElementIterator l(left.begin_elements());
            for (typename Matrix<DataType2_>::ConstElementIterator r(right.begin_non_zero_elements()),
                    r_end(right.end_non_zero_elements()); r != r_end ; )
            {
                while (l.index() < r.index())
                {
					*l = DataType1_(0);
                    ++l;
                }

                *l *= *r;
                ++l; ++r;
            }

            return left;
        }

        /**
         * Returns the resulting matrix after multiplying a BandedMatrix instance to a DenseMatrix instance elementwise.
         * \param left Reference to a DenseMatrix. Its return type is used for the result matrix.
         * \param right Reference to a BandedMatrix.
         **/
        template <typename DataType1_, typename DataType2_> static DenseMatrix<DataType1_> & value(DenseMatrix<DataType1_> & left, const BandedMatrix<DataType2_> & right)
        {
            if (left.columns() != right.columns())
            {
                throw MatrixColumnsDoNotMatch(right.columns(), left.columns());
            }

            if (left.rows() != right.rows())
            {
                throw MatrixRowsDoNotMatch(right.rows(), left.rows());
            }

            typename Matrix<DataType2_>::ConstElementIterator r(right.begin_elements());
            for (typename MutableMatrix<DataType1_>::ElementIterator l(left.begin_elements()),
                    l_end(left.end_elements()) ; l != l_end ; ++l)
            {
                *l *= *r;
                ++r;
            }

            return left;
        }

        /**
         * Returns the resulting matrix after multiplying a BandedMatrix instance to a SparseMatrix instance elementwise.
         * \param left Reference to a SparseMatrix. Its return type is used for the result matrix.
         * \param right Reference to a BandedMatrix.
         **/
        template <typename DataType1_, typename DataType2_> static SparseMatrix<DataType1_> & value(SparseMatrix<DataType1_> & left, const BandedMatrix<DataType2_> & right)
        {
            if (left.columns() != right.columns())
            {
                throw MatrixColumnsDoNotMatch(right.columns(), left.columns());
            }

            if (left.rows() != right.rows())
            {
                throw MatrixRowsDoNotMatch(right.rows(), left.rows());
            }

            typename Matrix<DataType2_>::ConstElementIterator r(right.begin_elements());
            for (typename MutableMatrix<DataType1_>::ElementIterator l(left.begin_non_zero_elements()),
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
