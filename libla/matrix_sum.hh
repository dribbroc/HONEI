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

/**
 * \file
 *
 * Templatized definitions of matrix sums.
 *
 * \ingroup grpmatrixoperations
 **/
namespace pg512
{

    /**
     * MatrixSum is the class template for the sum of two matrix instances.
     * \brief The first referenced matrix is changed under this operation.
     * \ingroup grpmatrixoperations
     **/
    template <typename Tag_ = tags::CPU> struct MatrixSum
    {
        /**
         * Returns the the resulting matrix of the sum of two DenseMatrix instances.
         *
         * \param left Reference to dense matrix that will be also used as result matrix.
         * \param right Reference to constant dense matrix to be added.
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
                *l += *r;
                ++r;
            }

            return left;
        }

		/**
         * Returns the the resulting matrix of the sum of a DenseMatrix and a SparseMatrix instance.
         *
         * \param left Reference to dense matrix that will be also used as result matrix.
         * \param right Reference to constant sparse matrix to be added.
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

                *l += *r;
                ++r;
            }

            return left;
        }

		/**
         * Returns the the resulting matrix of the sum of two given SparseMatrix instances.
         *
         * \param left Reference to sparse matrix that will be also used as result matrix.
         * \param right Reference to constant sparse matrix to be added.
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
                    l_end(left.end_elements()) ; l != l_end ; )
            {
				if (r.index() < l.index())
                {
                    left[r.row()][r.column()] = *r;
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
            return left;
        }

        /**
         * Returns the the resulting matrix of the sum of two given BandedMatrix instances.
         *
         * \param left Reference to banded matrix that will be also used as result matrix.
         * \param right Reference to constant banded matrix to be added.
         **/
        template <typename DataType1_, typename DataType2_> static BandedMatrix<DataType1_> & value(BandedMatrix<DataType1_> & left, const BandedMatrix<DataType2_> & right)
        {
            if (left.rows() != right.rows())
            {
                throw MatrixSizeDoesNotMatch(right.rows(), left.rows()); 
            }

            typename Matrix<DataType2_>::ConstElementIterator r(right.begin_elements());
            for (typename MutableMatrix<DataType1_>::ElementIterator l(left.begin_elements()),
                    l_end(left.end_elements()) ; l != l_end ; ++l)
            {
                *l += *r;
                ++r;
            }

            return left;
        }

        /**
         * Returns the the resulting matrix of the sum of a given BandedMatrix instance and a given DenseMatrix instance.
         *
         * \param left Reference to banded matrix that will be also used as result matrix.
         * \param right Reference to constant dense matrix to be subtracted.
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
                *l += *r;
                ++r;
            }

            return left;
        }

		/**
         * Returns the the resulting matrix of the sum of a given BandedMatrix instance and a given SparseMatrix instance.
         *
         * \param left Reference to banded matrix that will be also used as result matrix.
         * \param right Reference to constant sparse matrix to be subtracted.
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

            typename Matrix<DataType2_>::ConstElementIterator r(right.begin_non_zero_elements());
            for (typename MutableMatrix<DataType1_>::ElementIterator l(left.begin_elements()),
                    l_end(left.end_elements()) ; l != l_end ; )
            {
                while (l.index() < r.index() && (l != l_end))
                {
                    ++l;
                }

                *l += *r;
                ++r;
            }

            return left;
        }

        /**
         * Returns the the resulting matrix of the sum of a given DenseMatrix instance and a given BandedMatrix instance.
         *
         * \param left Reference to dense matrix that will be also used as result matrix.
         * \param right Reference to constant banded matrix to be subtracted.
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
                    l_end(left.end_elements()) ; l != l_end ; )
            {
                *l += *r;
                ++r;
            }

            return left;
        }

		/**
         * Returns the the resulting matrix of the sum of a given SparseMatrix instance and a given BandedMatrix instance.
         *
         * \param left Reference to sparse matrix that will be also used as result matrix.
         * \param right Reference to constant banded matrix to be subtracted.
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

            typename Matrix<DataType2_>::ConstElementIterator r(right.begin_elements()), r_end(right.end_elements());
            for (typename MutableMatrix<DataType1_>::ElementIterator l(left.begin_non_zero_elements()),
                    l_end(left.end_non_zero_elements()) ; l != l_end ; )
            {
				while (r.index() < l.index() && (r != r_end))
                {
                    ++r;
                }

                *l += *r;
                ++l;
            }
			///\todo: perhaps sparsify - i.e. addition of -7 and 7 possible
            return left;
        }

    };
}
#endif
