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

#ifndef LIBLA_GUARD_PRODUCT_HH
#define LIBLA_GUARD_PRODUCT_HH 1

#include <libla/banded_matrix.hh>
#include <libla/dense_vector.hh>
#include <libla/dense_matrix.hh>
#include <libla/dot_product.hh>
#include <libla/element_product.hh>
#include <libla/matrix_error.hh>
#include <libla/sparse_matrix.hh>
#include <libla/sparse_vector.hh>
#include <libla/sum.hh>
#include <libutil/tags.hh>

/**
 * \file
 *
 * Templatized definitions of operation Product.
 *
 * \ingroup grpoperations
 */
namespace honei
{
    /**
     * \brief Product of two entities.
     *
     * MatrixProduct is the class template for the product operation
     * \f[
     *     MatrixProduct(a, b): \quad c \leftarrow a * b,
     * \f]
     * which yields c, the product of entities a and b.
     *
     * The referenced containers are invariant under this operation.
     * In every case, a new object is created and returned.
     *
     * \ingroup grpoperations
     * \ingroup grpmatrixoperations
     * \ingroup grpvectoroperations
     */
    template <typename Tag_ = tags::CPU> struct Product
    {
        /**
         * \brief Returns the the product of two given entities.
         *
         * \param a The entity that is the first factor of the operation.
         * \param b The matrix that is the second factor of the operation.
         *
         * \retval c Will create a new object with Datatype of the first factor and return it.
         *
         * \exception MatrixSizeDoesNotMatch is thrown if two banded matrices do not have the same size.
         * \exception MatrixRowsDoNotMatch is thrown if two matrices do not have the same number of rows.
         * \exception MatrixColumnsDoNotMatch is thrown if two matrices do not have the same number of columns.
         * \exception MatrixIsNotSquare is thrown if a row access matrix's number of rows does not equal its number of columns.
         * \exception VectorSizeDoesNotMatch is thrown if two vectors do not have the same size.
         */

        /// \{

        template <typename DT1_, typename DT2_>
        static DenseVector<DT1_> value(const DenseMatrix<DT1_> & matrix, const DenseVector<DT2_> & vector)
        {
            CONTEXT("When multiplying DenseMatrix with DenseVector:");

            if (vector.size() != matrix.columns())
            {
                throw MatrixRowsDoNotMatch(matrix.columns(), vector.size());
            }

            DenseVector<DT1_> result(matrix.rows());
            typename Vector<DT1_>::ElementIterator l(result.begin_elements());
            for (unsigned long i=0; i < matrix.rows(); ++i)
            {
                DenseVector<DT1_> dv = matrix[i];
                *l = DotProduct<Tag_>::value(vector, dv);
                ++l;
            }

            return result;
        }

        template <typename DT1_, typename DT2_>
        static DenseVector<DT1_> value(const DenseMatrix<DT1_> & matrix, const SparseVector<DT2_> & vector)
        {
            CONTEXT("When multiplying DenseMatrix with SparseVector:");

            if (vector.size() != matrix.columns())
            {
                throw MatrixRowsDoNotMatch(matrix.columns(), vector.size());
            }

            DenseVector<DT1_> result(matrix.rows());
            typename Vector<DT1_>::ElementIterator l(result.begin_elements());
            for (unsigned long i=0; i < matrix.rows(); ++i)
            {
                DenseVector<DT1_> dv = matrix[i];
                *l = DotProduct<Tag_>::value(vector, dv);
                ++l;
            }

            return result;
        }

        template <typename DT1_, typename DT2_>
        static SparseVector<DT1_> value(const SparseMatrix<DT1_> & matrix, const DenseVector<DT2_> & vector)
        {
            CONTEXT("When multiplying SparseMatrix with DenseVector:");

            if (vector.size() != matrix.columns())
            {
                throw MatrixRowsDoNotMatch(matrix.columns(), vector.size());
            }

            SparseVector<DT1_> result(matrix.rows(),1);
            typename Vector<DT1_>::ElementIterator l(result.begin_elements());
            for (unsigned long i=0; i < matrix.rows(); ++i)
            {
                SparseVector<DT1_>  v(matrix[i]);
                *l = DotProduct<Tag_>::value(v, vector);
                ++l;
            }
            ///\todo: perhaps sparsify
            return result;
        }

        template <typename DT1_, typename DT2_>
        static SparseVector<DT1_> value(const SparseMatrix<DT1_> & matrix, const SparseVector<DT2_> & vector)
        {
            CONTEXT("When multiplying SparseMatrix with SparseVector:");

            if (vector.size() != matrix.columns())
            {
                throw MatrixRowsDoNotMatch(matrix.columns(), vector.size());
            }

            SparseVector<DT1_> result(matrix.rows(), matrix.rows());
            typename SparseVector<DT1_>::ElementIterator l(result.begin_elements());
            for (unsigned long i=0; i < matrix.rows(); ++i)
            {
                const Vector<DT1_> & v(matrix[i]);
                *l = DotProduct<Tag_>::value(vector, v);
                ++l;
            }
            ///\todo: perhaps sparsify
            return result;
        }

        template <typename DT1_, typename DT2_>
        static DenseVector<DT1_> value(const BandedMatrix<DT1_> & matrix, const DenseVector<DT2_> & vector)
        {
            CONTEXT("When multiplying BandedMatrix with DenseVector:");

            if (vector.size() != matrix.columns())
            {
                throw MatrixRowsDoNotMatch(matrix.columns(), vector.size());
            }

            DenseVector<DT1_> result(matrix.rows(), DT1_(0));
            for (typename BandedMatrix<DT1_>::ConstVectorIterator vi(matrix.begin_bands()), vi_end(matrix.end_bands()) ;
                    vi != vi_end ; ++vi)
            {
                typename Vector<DT2_>::ConstElementIterator j(vector.begin_elements()), j_end(vector.end_elements());
                typename Vector<DT1_>::ElementIterator r(result.begin_elements()), r_end(result.end_elements());

                int middle_index = matrix.rows() -1;
                // If we are above or on the diagonal band, we start at Element 0 and go on until Element band_size-band_index.
                if (vi.index() >= middle_index)
                {
                    for (unsigned int a=0 ; a < (vi.index()-middle_index) && j != j_end ; ++a)
                    {
                        ++j; // Get the right position in vector.
                    }

                    //Calculation of the element-index to stop in iteration!
                    unsigned long end = vi->size() - (vi.index() - middle_index);
                    for(typename Vector<DT1_>::ConstElementIterator b(vi->begin_elements()),
                            b_end(vi->element_at(end)) ; b != b_end ; ++b)
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
                    for (unsigned int a = 0; a < middle_index - vi.index() && r != r_end; ++a)
                    {
                        ++r; // Get the right position in result vector.
                    }
                    for(typename Vector<DT1_>::ConstElementIterator b(vi->element_at(start)),
                            b_end(vi->end_elements()) ; b != b_end ; ++b)
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

        template <typename DT1_, typename DT2_>
        static SparseVector<DT1_> value(const BandedMatrix<DT1_> & matrix, const SparseVector<DT2_> & vector)
        {
            CONTEXT("When multiplying BandedMatrix with SparseVector:");

            if (vector.size() != matrix.columns())
            {
                throw MatrixRowsDoNotMatch(matrix.columns(), vector.size());
            }

            SparseVector<DT1_> result(vector.size(), vector.capacity());
            for (typename BandedMatrix<DT1_>::ConstVectorIterator vi(matrix.begin_bands()),
                    vi_end(matrix.end_bands()) ;  vi != vi_end ; ++vi)
            {
                unsigned long middle_index = matrix.rows() -1;

                // If we are above or on the diagonal band, we first search the first non-zero element, then
                // we make positions in band and result vector meet it.
                if (vi.index() >= middle_index)
                {
                    unsigned long move_index = vi.index() - middle_index;
                    unsigned long end = vi->size() - move_index; //Calculation of the element-index to stop in iteration!

                    for (typename Vector<DT2_>::ConstElementIterator j(vector.begin_non_zero_elements()),
                            j_end(vector.end_non_zero_elements()) ; j != j_end ; ++j)
                    {
                        if (j.index() < move_index)
                        {
                            continue;
                        }

                        typename Vector<DT1_>::ConstElementIterator b(vi->begin_elements()), b_end(vi->element_at(end));
                        typename Vector<DT1_>::ElementIterator r(result.begin_elements());

                        while (b.index() < (j.index() - move_index) && b != b_end)
                        {
                            ++b;
                            ++r;
                        }

                        if (b != b_end)
                        {
                            *r += *b * *j;
                        }
                    }
                }
                // If we are below the diagonal band we correct the position in band and result vector.
                else
                {
                    unsigned long move_index = middle_index - vi.index(); // is also = the element index to start in iteration for b
                    for (typename Vector<DT2_>::ConstElementIterator j(vector.begin_non_zero_elements()),
                            j_end(vector.end_non_zero_elements()) ; j != j_end ; ++j)
                    {
                        typename Vector<DT1_>::ConstElementIterator b(vi->element_at(move_index)), b_end(vi->end_elements());
                        typename Vector<DT1_>::ElementIterator r(result.element_at(move_index));
                        while (b.index() < (j.index() + move_index) && b != b_end) // Need a positive index here, so + is used!
                        {
                            ++b;
                            ++r;
                        }

                        if (b != b_end)
                        {
                            *r += *b * *j;
                        }
                    }
                }
            }
            ///\todo: perhaps sparsify (*b can be zero)
            return result;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> value(const DenseMatrix<DT1_> & left, const DenseMatrix<DT2_> & right)
        {
            CONTEXT("When multiplying DenseMatrix with DenseMatrix:");

            if (left.columns() != right.rows())
                throw MatrixRowsDoNotMatch(right.rows(), left.columns());

            DenseMatrix<DT1_> result(right.columns(), left.rows());
            typename MutableMatrix<DT1_>::ElementIterator i(result.begin_elements());

            for (unsigned int s=0 ; s < left.rows() ; ++s)
            {
                const DenseVector<DT1_> left_row = left[s];
                for (unsigned int t=0; t < right.columns() ; ++t)
                {
                    const DenseVector<DT2_> right_column = right.column(t);
                    *i = DotProduct<>::value(left_row, right_column);
                    ++i;
                }

            }
            return result;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> value(const DenseMatrix<DT1_> & left, const SparseMatrix<DT2_> & right)
        {
            CONTEXT("When multiplying DenseMatrix with SparseMatrix:");

            if (left.columns() != right.rows())
                throw MatrixRowsDoNotMatch(right.rows(), left.columns());

            DenseMatrix<DT1_> result(right.columns(), left.rows(), DT1_(0));
            typename MutableMatrix<DT1_>::ElementIterator i(result.begin_elements());

            ///\todo: Should be optimized !!! (Use NonZeroIterators, less []-access ...)
            for (unsigned int l_row=0 ; l_row < left.rows() ; ++l_row)
            {
                const DenseVector<DT1_> left_row(left[l_row]);
                for (unsigned int r_column=0; r_column < right.columns() ; ++r_column)
                {
                    typename Vector<DT1_>::ConstElementIterator l(left_row.begin_elements());
                    for (unsigned int r_row=0; r_row < right.rows() ; ++r_row)
                    {
                        const SparseVector<DT2_> right_row(right[r_row]);
                        DT2_ right_value(right_row[r_column]);
                        //result[l_row][r_column] += right_value * *l;
                        *i += right_value * *l;
                        ++l;
                    }
                    ++i;
                }
            }
            return result;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> value(const SparseMatrix<DT1_> & left, const SparseMatrix<DT2_> & right)
        {
            CONTEXT("When multiplying SparseMatrix with SparseMatrix:");

            if (left.columns() != right.rows())
                throw MatrixRowsDoNotMatch(right.rows(), left.columns());

            SparseMatrix<DT1_> result(right.columns(), left.rows());
            typename MutableMatrix<DT1_>::ElementIterator i(result.begin_elements());
            ///\todo: Should be optimized !!! (Use NonZeroIterators, less []-access ...)
            for (unsigned int l_row=0 ; l_row < left.rows() ; ++l_row)
            {
                const SparseVector<DT1_> left_row(left[l_row]);
                for (unsigned int r_column=0; r_column < right.columns() ; ++r_column)
                {
                    typename Vector<DT1_>::ConstElementIterator l(left_row.begin_elements());
                    for (unsigned int r_row=0; r_row < right.rows() ; ++r_row)
                    {
                        const SparseVector<DT2_> right_row(right[r_row]);
                        DT2_ right_value(right_row[r_column]);
                        //result[l_row][r_column] += right_value * *l;
                        *i += right_value * *l;
                        ++l;
                    }
                    ++i;
                }
            }
            return result;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> value(const SparseMatrix<DT1_> & left, const DenseMatrix<DT2_> & right)
        {
            CONTEXT("when multiplying SparseMatrix with DenseMatrix:");

            if (left.columns() != right.rows())
                throw MatrixRowsDoNotMatch(right.rows(), left.columns());

            DenseMatrix<DT1_> result(right.columns(), left.rows());
            typename MutableMatrix<DT1_>::ElementIterator i(result.begin_elements());

            for (unsigned int s=0 ; s < left.rows() ; ++s)
            {
                const SparseVector<DT1_> left_row(left[s]);
                for (unsigned int t=0; t < right.columns() ; ++t)
                {
                    const DenseVector<DT2_> right_column = right.column(t);
                    *i = DotProduct<>::value(right_column, left_row);
                    ++i;
                }

            }
            return result;
        }

        template <typename DT1_, typename DT2_>
        static BandedMatrix<DT1_> value(const BandedMatrix<DT1_> & left, const BandedMatrix<DT2_> & right)
        {
            CONTEXT("When multiplying BandedMatrix with BandedMatrix:");

            if (left.columns() != right.rows())
                throw MatrixRowsDoNotMatch(right.rows(), left.columns());

            BandedMatrix<DT1_> result(left.rows());
            ///\todo: Implement when BandIterator ready
            return result;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT2_> value(const BandedMatrix<DT1_> & left, const DenseMatrix<DT2_> & right)
        {
            CONTEXT("When multiplying BandedMatrix with DenseMatrix:");

            if (left.columns() != right.rows())
                throw MatrixRowsDoNotMatch(right.rows(), left.columns());

            if (right.columns() != left.rows())
                throw MatrixRowsDoNotMatch(left.rows(), right.columns());

            DenseMatrix<DT2_> result(right.columns(), right.rows(), DT2_(0));
            unsigned long middle_index = left.size() -1;

            for (typename BandedMatrix<DT1_>::ConstVectorIterator vi(left.begin_bands()),
                    vi_end(left.end_bands()) ; vi != vi_end ; ++vi)
            {
                if (vi.index() == middle_index) // Are we on diagonal?
                {
                    for (unsigned int s = 0 ; s < right.columns() ; ++s)
                    {
                        Sum<>::value(result.column(s), ElementProduct<>::value(*(right.column(s).copy()), *vi));
                    }
                }
                else if (vi.index() > middle_index) // Are we above?
                {
                    for (unsigned int s = 0 ; s < right.columns() ; ++s)
                    {
                        DenseVector<DT2_> temp(right.rows(), DT2_(0)); // Temporary container for efficient calculation of elementwise vector product.
                        unsigned long real_index = vi.index() - middle_index;
                        typename Vector<DT2_>::ConstElementIterator a(right.column(s).element_at(real_index)),
                                 b(vi->begin_elements());
                        unsigned long end = temp.size();
                        end-= real_index;

                        for (typename Vector<DT2_>::ElementIterator x(temp.begin_elements()),
                                x_end(temp.element_at(end)) ; x != x_end ; ++x)
                        {
                            *x = *a * *b;
                            ++a; ++b;
                        }
                        Sum<>::value(result.column(s), temp);
                    }
                }
                else // We are below.
                {
                    for (unsigned int s = 0 ; s < right.columns() ; ++s)
                    {
                        // Temporary container for efficient calculation of elementwise vector product.
                        DenseVector<DT2_> temp(right.rows(), DT2_(0));
                        unsigned long real_index = middle_index - vi.index();
                        typename Vector<DT2_>::ConstElementIterator a(right.column(s).begin_elements()),
                                 b(vi->element_at(real_index));
                        for (typename Vector<DT2_>::ElementIterator x(temp.element_at(real_index)),
                                x_end(temp.end_elements()) ; x != x_end ; ++x)
                        {
                            *x = *a * *b;
                            ++a; ++b;
                        }
                        Sum<>::value(result.column(s), temp);
                    }
                }
            }
            return result;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> value(const BandedMatrix<DT1_> & left, const SparseMatrix<DT2_> & right)
        {
            CONTEXT("When multiplying BandedMatrix with SparseMatrix:");

            if (left.columns() != right.rows())
                throw MatrixRowsDoNotMatch(right.rows(), left.columns());

            if (right.columns() != left.rows())
                throw MatrixRowsDoNotMatch(left.rows(), right.columns());

            DenseMatrix<DT1_> result(left.columns(), left.rows(), DT1_(0));
            ///\todo: Will be implemented soon.
            return result;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> value(const DenseMatrix<DT1_> & left, const BandedMatrix<DT2_> & right)
        {
            CONTEXT("When multiplying DenseMatrix with BandedMatrix:");

            if (left.columns() != right.rows())
                throw MatrixRowsDoNotMatch(right.rows(), left.columns());

            DenseMatrix<DT1_> result(right.columns(), left.rows(), DT1_(0));
            ///\todo: Will be implemented soon.
            return result;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> value(const SparseMatrix<DT1_> & left, const BandedMatrix<DT2_> & right)
        {
            CONTEXT("When multiplying SparseMatrix with BandedMatrix:");

            if (left.columns() != right.rows())
                throw MatrixRowsDoNotMatch(right.rows(), left.columns());

            DenseMatrix<DT1_> result(right.columns(), left.rows(), DT1_(0));
            ///\todo: Will be implemented soon.
            return result;
        }

    };

#if 0
        /// Use the following algorithm for Cell processor, cause of optimized multiply-accumulate.

        template <> struct MatrixProduct<tags::Cell>
        {
            template <typename DT1_, typename DT2_> static DenseMatrix<DT1_> value(const DenseMatrix<DT1_> & left, const DenseMatrix<DT2_> & right)
            {
                if (left.columns() != right.rows())
                    throw MatrixRowsDoNotMatch(right.rows(), left.columns());

                DenseMatrix<DT1_> result(right.columns(), left.rows());
                typename MutableMatrix<DT1_>::ElementIterator i(result.begin_elements());

                for (unsigned int s=0 ; s < left.rows() ; ++s)
                {
                    const DenseVector<DT1_> left_row = left[s];
                    for (unsigned int t=0; t < right.columns() ; ++t)
                    {
                        const DenseVector<DT2_> right_column = right.column(t);
                        typename Vector<DT2_>::ConstElementIterator r(right_column.begin_elements());
                        for (typename Vector<DT1_>::ConstElementIterator l(left_row.begin_elements()),
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
#endif
}
#endif
