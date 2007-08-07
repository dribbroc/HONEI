/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Sven Mallach <sven.mallach@uni-dortmund.de>
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_MATRIX_PRODUCT_HH
#define LIBLA_GUARD_MATRIX_PRODUCT_HH 1

#include <libutil/tags.hh>
#include <libla/dense_matrix.hh>
#include <libla/sparse_matrix.hh>
#include <libla/scalar_product.hh>
#include <libla/banded_matrix.hh>
#include <libla/matrix_error.hh>

/**
 * \file
 *
 * Templatized definitions of matrix products.
 *
 * \ingroup grpmatrixoperations
 **/
namespace pg512
{

	///\todo: When DenseMatrix operator[] return type is restored, some cases can be put together to RowAccessMatrix.

    /**
     * MatrixProduct is the class template for multiplying two matrices
     * \brief The two referenced matrices are invariant under this operation.
     * \ingroup grpmatrixoperations
     **/
    template <typename Tag_ = tags::CPU> struct MatrixProduct
    {
        /**
         * Returns the resulting matrix after multiplying two DenseMatrix instances.
         * \param left Reference to a DenseMatrix used as first factor. Its return type is used for the result matrix.
         * \param right Reference to a DenseMatrix used as second factor.
         **/
        template <typename DataType1_, typename DataType2_> static DenseMatrix<DataType1_> value(const DenseMatrix<DataType1_> & left, const DenseMatrix<DataType2_> & right)
        {
            if (left.columns() != right.rows())
                throw MatrixRowsDoNotMatch(right.rows(), left.columns());

            DenseMatrix<DataType1_> result(right.columns(), left.rows());
            typename MutableMatrix<DataType1_>::ElementIterator i(result.begin_elements());

            for (unsigned int s=0 ; s < left.rows() ; ++s)
            {
                const DenseVector<DataType1_> left_row = left[s];
                for (unsigned int t=0; t < right.columns() ; ++t)
                {
                    const DenseVector<DataType2_> right_column = right.column(t);
					*i = ScalarProduct<>::value(left_row, right_column);
					++i;
                }

            }
            return result;
        }


		/**
         * Returns the resulting matrix after multiplying a DenseMatrix and a SparseMatrix instance.
         * \param left Reference to a DenseMatrix used as first factor. Its return type is used for the result matrix.
         * \param right Reference to a SparseMatrix used as second factor.
         **/

        template <typename DataType1_, typename DataType2_> static DenseMatrix<DataType1_> value(const DenseMatrix<DataType1_> & left, const SparseMatrix<DataType2_> & right)
        {
            if (left.columns() != right.rows())
                throw MatrixRowsDoNotMatch(right.rows(), left.columns());

            DenseMatrix<DataType1_> result(right.columns(), left.rows());
            typename MutableMatrix<DataType1_>::ElementIterator i(result.begin_elements());

            ///\todo: Should be optimized !!! (Use NonZeroIterators, less []-access ...)
            for (unsigned int l_row=0 ; l_row < left.rows() ; ++l_row)
            {
                const SparseVector<DataType1_> left_row(left[l_row]);


                for (unsigned int r_column=0; r_column < right.columns() ; ++r_column)
                {
                    typename Vector<DataType1_>::ConstElementIterator l(left_row.begin_elements());
                    for (unsigned int r_row=0; r_row < right.rows() ; ++r_row)
                    {
                        const SparseVector<DataType2_> right_row(right[r_row]);
                        DataType2_ right_value(right_row[r_column]);
                        //result[l_row][r_column] += right_value * *l;
                        *i += right_value * *l;
                        ++l;
                    }
                    ++i;
				}
            }
            return result;
        }

		 /**
         * Returns the resulting matrix after multiplying a sparse and a sparse matrix instance.
         * \param left Reference to a SparseMatrix used as first factor. Its return type is used for the result matrix.
         * \param right Reference to a SparseMatrix used as second factor.
         **/
        template <typename DataType1_, typename DataType2_> static SparseMatrix<DataType1_> value(const SparseMatrix<DataType1_> & left, const SparseMatrix<DataType2_> & right)
        {
            if (left.columns() != right.rows())
                throw MatrixRowsDoNotMatch(right.rows(), left.columns());

            SparseMatrix<DataType1_> result(right.columns(), left.rows());
            typename MutableMatrix<DataType1_>::ElementIterator i(result.begin_elements());
            ///\todo: Should be optimized !!! (Use NonZeroIterators, less []-access ...)
            for (unsigned int l_row=0 ; l_row < left.rows() ; ++l_row)
            {
                const SparseVector<DataType1_> left_row(left[l_row]);


                for (unsigned int r_column=0; r_column < right.columns() ; ++r_column)
                {
                    typename Vector<DataType1_>::ConstElementIterator l(left_row.begin_elements());
                    for (unsigned int r_row=0; r_row < right.rows() ; ++r_row)
                    {
                        const SparseVector<DataType2_> right_row(right[r_row]);
                        DataType2_ right_value(right_row[r_column]);
                        //result[l_row][r_column] += right_value * *l;
                        *i += right_value * *l;
                        ++l;
                    }
                    ++i;
				}
            }
            return result;
        }

		/**
         * Returns the resulting matrix after multiplying a SparseMatrix and a DenseMatrix instance.
         * \param left Reference to a SparseMatrix used as first factor.
         * \param right Reference to a DenseMatrix used as second factor. Its return type is used for the result matrix.
         **/
        template <typename DataType1_, typename DataType2_> static DenseMatrix<DataType1_> value(const SparseMatrix<DataType1_> & left, const DenseMatrix<DataType2_> & right)
        {
            if (left.columns() != right.rows())
                throw MatrixRowsDoNotMatch(right.rows(), left.columns());

            DenseMatrix<DataType1_> result(right.columns(), left.rows());
            typename MutableMatrix<DataType1_>::ElementIterator i(result.begin_elements());

            for (unsigned int s=0 ; s < left.rows() ; ++s)
            {
                const SparseVector<DataType1_> left_row(left[s]);
                for (unsigned int t=0; t < right.columns() ; ++t)
                {
                    const DenseVector<DataType2_> right_column = right.column(t);
					*i = ScalarProduct<>::value(right_column, left_row);
					++i;
                }

            }
            return result;
        }

		/**
         * Returns the resulting matrix after multiplying two BandedMatrix instances.
         * \param left Reference to a BandedMatrix used as first factor. Its return type is used for the result matrix.
         * \param right Reference to a BandedMatrix used as second factor.
         **/
        template <typename DataType1_, typename DataType2_> static BandedMatrix<DataType1_> value(const BandedMatrix<DataType1_> & left, const BandedMatrix<DataType2_> & right)
        {
            if (left.columns() != right.rows())
                throw MatrixRowsDoNotMatch(right.rows(), left.columns());

            BandedMatrix<DataType1_> result(left.rows());
            ///\todo: Implement when BandIterator ready
            return result;
        }

        /**
         * Returns the resulting matrix after multiplying a BandedMatrix instance and a DenseMatrix instance.
         * \param left Reference to a BandedMatrix used as first factor. Its return type is used for the result matrix.
         * \param right Reference to a DenseMatrix used as second factor.
         **/
        template <typename DataType1_, typename DataType2_> static BandedMatrix<DataType1_> value(const BandedMatrix<DataType1_> & left, const DenseMatrix<DataType2_> & right)
        {
            if (left.columns() != right.rows())
                throw MatrixRowsDoNotMatch(right.rows(), left.columns());

            BandedMatrix<DataType1_> result(left.rows());
            ///\todo: Implement when BandIterator ready
            return result;
        }

		/**
         * Returns the resulting matrix after multiplying a BandedMatrix instance and a SparseMatrix instance.
         * \param left Reference to a BandedMatrix used as first factor. Its return type is used for the result matrix.
         * \param right Reference to a SparseMatrix used as second factor.
         **/
        template <typename DataType1_, typename DataType2_> static BandedMatrix<DataType1_> value(const BandedMatrix<DataType1_> & left, const SparseMatrix<DataType2_> & right)
        {
            if (left.columns() != right.rows())
                throw MatrixRowsDoNotMatch(right.rows(), left.columns());

            BandedMatrix<DataType1_> result(left.rows());
            ///\todo: Implement when BandIterator ready
            return result;
        }


        /**
         * Returns the resulting matrix after multiplying a DenseMatrix instance and a BandedMatrix instance.
         * \param left Reference to a DenseMatrix used as first factor. Its return type is used for the result matrix.
         * \param right Reference to a BandedMatrix used as second factor.
         **/
        template <typename DataType1_, typename DataType2_> static DenseMatrix<DataType1_> value(const DenseMatrix<DataType1_> & left, const BandedMatrix<DataType2_> & right)
        {
            if (left.columns() != right.rows())
                throw MatrixRowsDoNotMatch(right.rows(), left.columns());

            DenseMatrix<DataType1_> result(right.columns(), left.rows());
            ///\todo: Implement when BandIterator ready
            return result;
        }

		/**
         * Returns the resulting matrix after multiplying a SparseMatrix instance and a BandedMatrix instance.
         * \param left Reference to a SparseMatrix used as first factor.
         * \param right Reference to a BandedMatrix used as second factor. Its return type is used for the result matrix.
         **/
        template <typename DataType1_, typename DataType2_> static DenseMatrix<DataType1_> value(const SparseMatrix<DataType1_> & left, const BandedMatrix<DataType2_> & right)
        {
            if (left.columns() != right.rows())
                throw MatrixRowsDoNotMatch(right.rows(), left.columns());

            DenseMatrix<DataType1_> result(right.columns(), left.rows());
            ///\todo: Implement when BandIterator ready
            return result;
        }

    };

	/// Use the following algorithm for Cell processor, cause of optimized multiply-accumulate.

	template <> struct MatrixProduct<tags::Cell>
	{
		/**
         * Returns the resulting matrix after multiplying two DenseMatrix instances.
		 * \param left Reference to a DenseMatrix used as first factor. Its return type is used for the result matrix.
		 * \param right Reference to a DenseMatrix used as second factor.
         **/
        template <typename DataType1_, typename DataType2_> static DenseMatrix<DataType1_> value(const DenseMatrix<DataType1_> & left, const DenseMatrix<DataType2_> & right)
        {
            if (left.columns() != right.rows())
                throw MatrixRowsDoNotMatch(right.rows(), left.columns());

            DenseMatrix<DataType1_> result(right.columns(), left.rows());
            typename MutableMatrix<DataType1_>::ElementIterator i(result.begin_elements());

            for (unsigned int s=0 ; s < left.rows() ; ++s)
            {
                const DenseVector<DataType1_> left_row = left[s];
                for (unsigned int t=0; t < right.columns() ; ++t)
                {
					const DenseVector<DataType2_> right_column = right.column(t);
					typename Vector<DataType2_>::ConstElementIterator r(right_column.begin_elements());
					for (typename Vector<DataType1_>::ConstElementIterator l(left_row.begin_elements()),
							l_end(left_row.end_elements()) ; l != l_end ; ++l, ++r)
					{
						*i += (*l) * (*r);
					}
					++i;
				}
            }
            return result;
        }
	};


}
#endif
