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

/**
 * \file
 *
 * Templatized definitions of matrix differences.
 *
 * \ingroup grpmatrixoperations
 **/
namespace pg512
{

    /**
     * MatrixDifference is the class template for the difference of two matrix.
     * \brief The first referenced matrix is changed under this operation.
     * \ingroup grpmatrixoperations
     **/
    template <typename Tag_ = tags::CPU> struct MatrixDifference
    {
        /**
         * Returns the the resulting matrix of the difference of two dense matrix instances.
         *
         * \param left Reference to dense matrix that will be also used as result matrix.
         * \param right Reference to constant dense matrix to be subtracted.
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
                    l_end(left.end_elements()) ; l != l_end ; ++l, ++r)
            {
                *l -= *r;
            }

            return left;
        }

		/**
         * Returns the the resulting matrix of the difference of a DenseMatrix and a SparseMatrix instance.
         *
         * \param left Reference to dense matrix that will be also used as result matrix.
         * \param right Reference to constant sparse matrix to be subtracted.
         **/
        template <typename DataType1_, typename DataType2_> static DenseMatrix<DataType1_> & value(DenseMatrix<DataType1_> & left, const SparseMatrix<DataType2_> & right)
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
            for (typename MutableMatrix<DataType1_>::ElementIterator l(left.begin_elements()),
                    l_end(left.end_elements()) ; l != l_end ; )
            {
				while (l.index() < r.index() && (l != l_end))
                {
                    ++l;
                }

                *l -= *r;
                ++r;
            }

            return left;
        }

		/**
         * Returns the the resulting matrix of the difference of a sparse matrix instance and a sparse matrix instance.
         *
         * \param left Reference to sparse matrix that will be also used as result matrix.
         * \param right Reference to constant sparse matrix to be subtracted.
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
                if (r.index() < l.index())
                {
                    left[r.row()](r.column) = -(*r);
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
			///\todo: perhaps sparsify - i.e. addition of -7 and 7 possible
            return left;
        }

        /**
         * Returns the the resulting matrix of the difference of two given BandedMatrix instances.
         *
         * \param left Reference to banded matrix that will be also used as result matrix.
         * \param right Reference to constant banded matrix to be subtracted.
         **/
        template <typename DataType1_, typename DataType2_> static BandedMatrix<DataType1_> & value(BandedMatrix<DataType1_> & left, const BandedMatrix<DataType2_> & right)
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
                    l_end(left.end_elements()) ; l != l_end ; ++l, ++r)
            {
                *l -= *r;
            }

            return left;
        }

		///\todo: Find a solution for different order and return types when changing RowAccess to Dense/Sparse, because now this would mean duplicate signatures.

        /**
         * Returns the the resulting matrix of the difference of a given BandedMatrix instance and a given RowAccessMatrix instance.
         *
         * \param left Reference to banded matrix that will be also used as result matrix.
         * \param right Reference to constant row access matrix to be subtracted.
         **/
        template <typename DataType1_, typename DataType2_> static BandedMatrix<DataType1_> & value(BandedMatrix<DataType1_> & left, const RowAccessMatrix<DataType2_> & right)
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
                    l_end(left.end_elements()) ; l != l_end ; ++l, ++r)
            {
                *l -= *r;
            }

            return left;
        }

        /**
         * Returns the the resulting matrix of the difference of a given RowAccessMatrix instance and a given BandedMatrix instance.
         *
         * \param left Reference to row access matrix that will be also used as result matrix.
         * \param right Reference to constant banded matrix to be subtracted.
         **/
        template <typename DataType1_, typename DataType2_> static BandedMatrix<DataType1_> & value(const RowAccessMatrix<DataType2_> & left, BandedMatrix<DataType1_> & right)
        {
            if (left.columns() != right.columns())
            {
                throw MatrixColumnsDoNotMatch(right.columns(), left.columns());
            }

            if (left.rows() != right.rows())
            {
                throw MatrixRowsDoNotMatch(right.rows(), left.rows());
            }

            typename Matrix<DataType1_>::ConstElementIterator l(left.begin_elements());
            for (typename MutableMatrix<DataType2_>::ElementIterator r(right.begin_elements()),
                    r_end(right.end_elements()) ; r != r_end ; ++r, ++l)
            {
                *r = *l - *r;
            }

            return right;
        }

        /**
         * Returns the the resulting matrix of the difference of a given BandedMatrix instance and a given DenseMatrix instance.
         *
         * \param left Reference to constant banded matrix which is the base.
         * \param right Reference to dense matrix to be substracted that will be also used as result matrix.
         **/
        template <typename DataType1_, typename DataType2_> static DenseMatrix<DataType1_> & value(const BandedMatrix<DataType2_> & left, DenseMatrix<DataType1_> & right)
        {
            if (left.columns() != right.columns())
            {
                throw MatrixColumnsDoNotMatch(right.columns(), left.columns());
            }

            if (left.rows() != right.rows())
            {
                throw MatrixRowsDoNotMatch(right.rows(), left.rows());
            }

            typename Matrix<DataType1_>::ConstElementIterator l(left.begin_elements());
            for (typename MutableMatrix<DataType2_>::ElementIterator r(right.begin_elements()),
                    r_end(right.end_elements()) ; r != r_end ; ++r, ++l)
            {
                *r = *l - *r;
            }

            return right;
        }
		
		 /**
         * Returns the the resulting matrix of the difference of a given BandedMatrix instance and a given SparseMatrix instance.
         *
         * \param left Reference to constant banded matrix which is the base.
         * \param right Reference to sparse matrix to be substracted that will be also used as result matrix.
         **/
        template <typename DataType1_, typename DataType2_> static SparseMatrix<DataType1_> & value(const BandedMatrix<DataType2_> & left, SparseMatrix<DataType1_> & right)
        {
            if (left.columns() != right.columns())
            {
                throw MatrixColumnsDoNotMatch(right.columns(), left.columns());
            }

            if (left.rows() != right.rows())
            {
                throw MatrixRowsDoNotMatch(right.rows(), left.rows());
            }

            typename Matrix<DataType1_>::ConstElementIterator l(left.begin_elements()), l_end(left.begin_elements());
            for (typename MutableMatrix<DataType2_>::ElementIterator r(right.begin_non_zero_elements()),
                    r_end(right.end_non_zero_elements()) ; r != r_end ; ++r, ++l)
            {
				while (l.index() < r.index() && (l != l_end))
                {
                    ++l;
                }

				*r = *l - *r;
                ++r;
            }
			///\todo: perhaps sparsify - i.e. addition of -7 and 7 possible
            return right;
        }

    };
}
#endif
