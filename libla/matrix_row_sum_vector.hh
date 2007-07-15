/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_MATRIX_ROW_SUM_VECTOR_HH
#define LIBLA_GUARD_MATRIX_ROW_SUM_VECTOR_HH 1

#include <libla/dense_vector.hh>
#include <libla/matrix.hh>
#include <libla/banded_matrix.hh>
#include <libla/vector_element_sum.hh>

#include <tr1/memory>

namespace pg512
{
    /**
     * A MatrixRowSumVector yields a Vector of which each element holds the sum
     * of the matrix's corresponding row's elements.
      * \ingroup grpmatrixoperations
     **/
    template <typename Tag_ = tags::CPU> struct MatrixRowSumVector
    {
        template <typename DataType_> static DenseVector<DataType_> value(const RowAccessMatrix<DataType_> & matrix)
        {
            DenseVector<DataType_> result(matrix.rows());

            for (typename Vector<DataType_>::ElementIterator i(result.begin_elements()), i_end(result.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = VectorElementSum<DataType_, Tag_>::value(matrix[i.index()]);
            }

            return result;
        }


        template <typename DataType_> static DenseVector<DataType_> value(const BandedMatrix<DataType_> & matrix)
        {
            DenseVector<DataType_> result(matrix.rows(), DataType_(0));
            typename Matrix<DataType_>::ConstElementIterator b(matrix.begin_elements()), b_end(matrix.end_elements());

            for (typename Vector<DataType_>::ElementIterator a(result.begin_elements()), a_end(result.end_elements()); a != a_end; ++a)
            {
                for (b ; b.row() == a.index() && b != b_end; ++b)
                {
                    *a += *b;
                }
            }

            return result;
        }
    };
}

#endif
