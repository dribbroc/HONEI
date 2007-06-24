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

#ifndef LIBLA_GUARD_MATRIX_PRODUCT_HH
#define LIBLA_GUARD_MATRIX_PRODUCT_HH 1

#include <libutil/tags.hh>
#include <libla/dense_matrix.hh>
#include <libla/banded_matrix.hh>
#include <matrix_error.hh>

/**
 * \file
 *
 * Templatized definitions of matrix products.<br/>
 *
 * \ingroup grpmatrixoperations
 **/
namespace pg512
{

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

            DenseMatrix<DataType1_> result(right.columns(), left.rows(), DataType1_(0));
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
