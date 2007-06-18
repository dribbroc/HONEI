/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
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

#ifndef LIBLA_GUARD_MATRIX_ROW_SUM_VECTOR_HH
#define LIBLA_GUARD_MATRIX_ROW_SUM_VECTOR_HH 1

#include <libla/dense_vector.hh>
#include <libla/matrix.hh>
#include <libla/vector_element_sum.hh>

#include <tr1/memory>

namespace pg512
{
    /**
     * A MatrixRowSumVector yields a Vector of which each element holds the sum
     * of the matrix's corresponding row's elements.
      * \ingroup grpmatrixoperations
     **/
    template <typename DataType_, typename Tag_ = tags::CPU> struct MatrixRowSumVector
    {
        static std::tr1::shared_ptr<DenseVector<DataType_> > value(const RowAccessMatrix<DataType_> & matrix)
        {
            std::tr1::shared_ptr<DenseVector<DataType_> > result(new DenseVector<DataType_>(matrix.rows()));

            for (typename Vector<DataType_>::ElementIterator i(result->begin_elements()), i_end(result->end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = VectorElementSum<DataType_, Tag_>::value(matrix[i.index()]);
            }

            return result;
        }
    };
}

#endif
