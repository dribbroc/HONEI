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

#include <libla/tags.hh>
#include <libla/dense_matrix.hh>
#include <libla/banded_matrix.hh>
#include <libla/matrix_error.hh>

/**
 * \file
 *
 * Templatized definitions of matrix sums.<br/>
 *
 * \ingroup grpmatrixoperations
 **/
namespace pg512
{

    /**
     * MatrixSum is the class template for the sum of two matrix.
     * \brief The first referenced matrix is changed under this operation.
     * \ingroup grpmatrixoperations
     **/
    template <typename Tag_ = tags::CPU> struct MatrixSum
    {
        /**
         * Returns the the resulting matrix of the sum of two given DenseMatrix instances.
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

            if (left.rows() != left.rows())
            {
                throw MatrixRowsDoNotMatch(right.rows(), left.rows());
            }

            typename Matrix<DataType2_>::ConstElementIterator r(right.begin_elements());
            for (typename MutableMatrix<DataType1_>::ElementIterator l(left.begin_elements()),
                    l_end(left.end_elements()) ; l != l_end ; ++l)
            {
                *l += *r;
                r++;
            }

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
            if (left.columns() != right.columns())
            {
                throw MatrixColumnsDoNotMatch(right.columns(), left.columns());
            }

            if (left.rows() != left.rows())
            {
                throw MatrixRowsDoNotMatch(right.rows(), left.rows());
            }

            typename Matrix<DataType2_>::ConstElementIterator r(right.begin_elements());
            for (typename MutableMatrix<DataType1_>::ElementIterator l(left.begin_elements()),
                    l_end(left.end_elements()) ; l != l_end ; ++l)
            {
                *l += *r;
                r++;
            }

            return left;
        }

    };
}
#endif
