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

#ifndef LIBLA_GUARD_VECTOR_MATRIX_PRODUCT_HH
#define LIBLA_GUARD_VECTOR_MATRIX_PRODUCT_HH 1

#include <libutil/tags.hh>
#include <libla/dense_vector.hh>
#include <libla/sparse_vector.hh>
#include <libla/dense_matrix.hh>
#include <libla/banded_matrix.hh>
#include <libla/matrix_error.hh>
#include <libla/scalar_product.hh>

/**
 * \file
 *
 * Templatized definitions of vector-matrix products.<br/>
 *
 * \ingroup grpvectoroperations
 **/
namespace pg512
{

    /**
     * VectorMatrixProduct is the class template for multiplying a matrix to a vector
     * \brief The referenced vector and matrix are invariant under this operation.
     * \ingroup grpvectoroperations
     * \ingroup grpmatrixoperations
     **/
    template <typename Tag_ = tags::CPU> struct VectorMatrixProduct
    {
        /**
         * Returns the resulting vector after multiplying a matrix to a given DenseVector instance.
         * \param vector DenseVector to be used as factor.
         * \param matrix The DenseMatrix to be used as factor.
         **/
        template <typename DataType1_, typename DataType2_> static DenseVector<DataType1_> value(const DenseVector<DataType1_> & vector, const DenseMatrix<DataType2_> & matrix)
        {
            if (vector.size() != matrix.rows())
                throw MatrixRowsDoNotMatch(matrix.rows(), vector.size());

            DenseVector<DataType1_> result(vector.size());
            typename Vector<DataType1_>::ElementIterator l(result.begin_elements());
            for (unsigned long i=0; i < matrix.columns(); ++i)
            {
                DenseVector<DataType2_> dv = matrix.column(i);
                *l = ScalarProduct<Tag_>::value(vector, dv);
                ++l;
            }

            return result;
        }

        /**
         * Returns the resulting vector after multiplying a matrix to a given SparseVector instance.
         * \param vector SparseVector to be used as factor.
         * \param matrix The DenseMatrix to be used as factor.
         **/
        template <typename DataType1_, typename DataType2_> static SparseVector<DataType1_> value(const SparseVector<DataType1_> & vector, const DenseMatrix<DataType2_> & matrix)
        {
            if (vector.size() != matrix.rows())
                throw MatrixRowsDoNotMatch(matrix.rows(), vector.size());

            SparseVector<DataType1_> result(vector.size());
            typename Vector<DataType1_>::ElementIterator l(result.begin_elements());
            for (unsigned long i=0; i < matrix.columns(); ++i)
            {
                SparseVector<DataType2_> dv = matrix.column(i);
                *l = ScalarProduct<Tag_>::value(vector, dv);
                ++l;
            }

            return result;
        }

    };
}
#endif
