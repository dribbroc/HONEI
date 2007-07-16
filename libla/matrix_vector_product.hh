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

#ifndef LIBLA_GUARD_MATRIX_VECTOR_PRODUCT_HH
#define LIBLA_GUARD_MATRIX_VECTOR_PRODUCT_HH 1

#include <libutil/tags.hh>
#include <libla/dense_vector.hh>
#include <libla/sparse_vector.hh>
#include <libla/dense_matrix.hh>
#include <libla/banded_matrix.hh>
#include <libla/sparse_matrix.hh>
#include <libla/matrix_error.hh>
#include <libla/scalar_product.hh>

/**
 * \file
 *
 * Templatized definitions of matrix-vector products.
 *
 * \ingroup grpmatrixoperations
 * \ingroup grpvectoroperations
 **/
namespace pg512
{

    /**
     * MatrixVectorProduct is the class template for multiplying a matrix to a vector
     * \brief The referenced matrix and vector are invariant under this operation.
     * \ingroup grpmatrixoperations
     * \ingroup grpvectoroperations
     **/
    template <typename Tag_ = tags::CPU> struct MatrixVectorProduct
    {
        /**
         * Returns the resulting vector after multiplying a DenseVector to a given DenseMatrix instance.
         * \param matrix The DenseMatrix to be used as factor.
         * \param vector DenseVector to be used as factor.
         **/
        template <typename DataType1_, typename DataType2_> static DenseVector<DataType1_> value(const DenseMatrix<DataType1_> & matrix, const DenseVector<DataType2_> & vector)
        {
            if (vector.size() != matrix.columns())
                {
                throw MatrixRowsDoNotMatch(matrix.columns(), vector.size());
                }

            DenseVector<DataType1_> result(matrix.rows());
            typename Vector<DataType1_>::ElementIterator l(result.begin_elements());
            for (unsigned long i=0; i < matrix.rows(); ++i)
            {
                DenseVector<DataType1_> dv = matrix[i];
                *l = ScalarProduct<Tag_>::value(vector, dv);
                ++l;
            }

            return result;
        }

		/**
         * Returns the resulting vector after multiplying a SparseVector to a given DenseMatrix instance.
         * \param matrix The DenseMatrix to be used as factor.
         * \param vector SparseVector to be used as factor.
         **/
        template <typename DataType1_, typename DataType2_> static DenseVector<DataType1_> value(const DenseMatrix<DataType1_> & matrix, const SparseVector<DataType2_> & vector)
        {
            if (vector.size() != matrix.columns())
                {
                throw MatrixRowsDoNotMatch(matrix.columns(), vector.size());
                }

            DenseVector<DataType1_> result(matrix.rows());
            typename Vector<DataType1_>::ElementIterator l(result.begin_elements());
            for (unsigned long i=0; i < matrix.rows(); ++i)
            {
                DenseVector<DataType1_> dv = matrix[i];
                *l = ScalarProduct<Tag_>::value(vector, dv);
                ++l;
            }

            return result;
        }

		/**
         * Returns the resulting vector after multiplying a DenseVector to a given SparseMatrix instance.
         * \param matrix The SparseMatrix to be used as factor.
         * \param vector DenseVector to be used as factor.
         **/
        template <typename DataType1_, typename DataType2_> static SparseVector<DataType1_> value(const SparseMatrix<DataType1_> & matrix, const DenseVector<DataType2_> & vector)
        {
            if (vector.size() != matrix.columns())
                {
                throw MatrixRowsDoNotMatch(matrix.columns(), vector.size());
                }

            SparseVector<DataType1_> result(matrix.rows());
            typename Vector<DataType1_>::ElementIterator l(result.begin_elements());
            for (unsigned long i=0; i < matrix.rows(); ++i)
            {
                Vector<DataType1_> & v(matrix[i]);
                *l = ScalarProduct<Tag_>::value(vector, v);
                ++l;
            }
			///\todo: perhaps sparsify
            return result;
        }

        /**
         * Returns the resulting vector after multiplying a SparseVector to a SparseMatrix instance.
         * \param matrix The SparseMatrix to be used as factor.
         * \param vector SparseVector to be used as factor.
         **/
        template <typename DataType1_, typename DataType2_> static SparseVector<DataType1_> value(const SparseMatrix<DataType1_> & matrix, const SparseVector<DataType2_> & vector)
        {
            if (vector.size() != matrix.columns())
                throw MatrixRowsDoNotMatch(matrix.columns(), vector.size());

            SparseVector<DataType1_> result(matrix.rows(), matrix.rows());
            typename Vector<DataType1_>::ElementIterator l(result.begin_elements());
            for (unsigned long i=0; i < matrix.rows(); ++i)
            {
                const Vector<DataType1_> & v(matrix[i]);
                *l = ScalarProduct<Tag_>::value(vector, v);
                ++l;
            }

            return result;
        }

        /**
         * Returns the resulting vector after multiplying a DenseVector to a given BandedMatrix instance.
         * \param matrix The BandedMatrix to be used as factor.
         * \param vector DenseVector to be used as factor.
         **/
        template <typename DataType1_, typename DataType2_> static DenseVector<DataType1_> value(const BandedMatrix<DataType1_> & matrix, const DenseVector<DataType2_> & vector)
        {
            if (vector.size() != matrix.columns())
                throw MatrixRowsDoNotMatch(matrix.columns(), vector.size());

            DenseVector<DataType1_> result(matrix.rows());
            ///\todo: Implement when band-iterator available.
            /*
            typename Vector<DataType1_>::ElementIterator l(result.begin_elements());
            for (unsigned long i=0; i < matrix.columns(); ++i)
            {
                DenseVector<DataType1_> dv = matrix.column(i);
                *l = ScalarProduct<Tag_>::value(vector, dv);
                ++l;
            }
            */
            return result;
        }

        /**
         * Returns the resulting vector after multiplying a SparseVector to a given BandedMatrix instance.
         * \param matrix The BandedMatrix to be used as factor.
         * \param vector SparseVector to be used as factor.
         **/
        template <typename DataType1_, typename DataType2_> static SparseVector<DataType1_> value(const BandedMatrix<DataType1_> & matrix, const SparseVector<DataType2_> & vector)
        {
            if (vector.size() != matrix.columns())
                throw MatrixRowsDoNotMatch(matrix.columns(), vector.size());

            SparseVector<DataType1_> result(matrix.rows(), matrix.rows());
            ///\todo: Implement when band-iterator available.
            /*
            typename Vector<DataType1_>::ElementIterator l(result.begin_elements());
            for (unsigned long i=0; i < matrix.columns(); ++i)
            {
                DenseVector<DataType1_> dv = matrix.column(i);
                *l = ScalarProduct<Tag_>::value(vector, dv);
                ++l;
            }
            */
            return result;
        }



    };
}
#endif
