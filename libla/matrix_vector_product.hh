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
#include <iostream>

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

            SparseVector<DataType1_> result(matrix.rows(),1);
            typename Vector<DataType1_>::ElementIterator l(result.begin_elements());
            for (unsigned long i=0; i < matrix.rows(); ++i)
            {
                SparseVector<DataType1_>  v(matrix[i]);
                *l = ScalarProduct<Tag_>::value(v, vector);
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
            {
                throw MatrixRowsDoNotMatch(matrix.columns(), vector.size());
            }

            SparseVector<DataType1_> result(matrix.rows(), matrix.rows());
            typename SparseVector<DataType1_>::ElementIterator l(result.begin_elements());
            for (unsigned long i=0; i < matrix.rows(); ++i)
            {
                const Vector<DataType1_> & v(matrix[i]);
                *l = ScalarProduct<Tag_>::value(vector, v);
                ++l;
            }
			///\todo: perhaps sparsify
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
            {
                throw MatrixRowsDoNotMatch(matrix.columns(), vector.size());
            }

            DenseVector<DataType1_> result(matrix.rows(), DataType1_(0));
            for (typename BandedMatrix<DataType1_>::ConstVectorIterator vi(matrix.begin_bands()), vi_end(matrix.end_bands()) ;  vi != vi_end ; ++vi)
            {
                DenseVector<DataType1_> band = *vi;
                typename Vector<DataType2_>::ConstElementIterator j(vector.begin_elements()), j_end(vector.end_elements());
                typename Vector<DataType1_>::ElementIterator r(result.begin_elements()), r_end(result.end_elements());

                int middle_index = matrix.rows() -1;
                // If we are above or on the diagonal band, we start at Element 0 and go on until Element band_size-band_index.
                if (vi.index() >= middle_index)
                {
                    for (unsigned int a=0 ; a < (vi.index()-middle_index) && j != j_end ; ++a)
                    {
                        ++j; // Get the right position in vector.
                    }
                    unsigned long end = band.size() - (vi.index() - middle_index); //Calculation of the element-index to stop in iteration!
                    for(typename Vector<DataType1_>::ConstElementIterator b(band.begin_elements()), b_end(band.element_at(end)) ; b != b_end ; ++b)
                    {
                        *r += *b * *j;
                        if (j != j_end && r != r_end)
                        {
                            ++r;
                            ++j;
                        }
                    }
                }
                // If we are below the diagonal band, we start at Element index and go on until the last element.
                else
                {
                    unsigned long start = middle_index - vi.index(); //Calculation of the element-index to start in iteration!
                    for (unsigned int a=0; a < middle_index - vi.index() && r != r_end; ++a)
                    {
                        ++r; // Get the right position in result vector.
                    }
                    for(typename Vector<DataType1_>::ConstElementIterator b(band.element_at(start)), b_end(band.end_elements()) ; b != b_end ; ++b)
                    {
                        *r += *b * *j;
                        if (j != j_end && r != r_end)
                        {
                            ++j;
                            ++r;
                        }
                    }
                }
            }

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
            {
                throw MatrixRowsDoNotMatch(matrix.columns(), vector.size());
            }

            SparseVector<DataType1_> result(vector.size(), vector.capacity());
            for (typename BandedMatrix<DataType1_>::ConstVectorIterator vi(matrix.begin_bands()), vi_end(matrix.end_bands()) ;  vi != vi_end ; ++vi)
            {
                DenseVector<DataType1_> band = *vi;
                typename Vector<DataType1_>::ConstElementIterator b(band.begin_elements()), b_end(band.end_elements());;
                typename Vector<DataType1_>::ElementIterator r(result.begin_elements());
                int middle_index = matrix.rows() -1;

                // If we are above or on the diagonal band, we first search the first non-zero element, then we make positions in band and result vector meet it.
                if (vi.index() >= middle_index)
                {
                    for (typename Vector<DataType2_>::ConstElementIterator j(vector.begin_non_zero_elements()), j_end(vector.end_non_zero_elements()) ; j != j_end ; ++j)
                    {
                        while (j.index() < vi.index()-middle_index && j != j_end)
                        {
                            ++j;
                        }

                        while (b.index() < j.index() - (vi.index()-middle_index) && b != b_end)
                        {
                            ++b;
                            ++r;
                        }

                        *r += *b * *j;
                    }
                }
                // If we are below the diagonal band we correct the position in band an result vector.
                else
                {
                    for (typename Vector<DataType2_>::ConstElementIterator j(vector.begin_non_zero_elements()), j_end(vector.end_non_zero_elements()) ; j != j_end ; ++j)
                    {
                        while (b.index() < j.index() - (vi.index() - middle_index) && b != b_end) //normally + but we get a negative index where we want to work with a positive one
                        {
                            ++b;
                            ++r;
                        }

                        *r += *b * *j;
                    }
                }
            }
			///\todo: perhaps sparsify (*b can be zero)
            return result;
        }
    };
}
#endif
