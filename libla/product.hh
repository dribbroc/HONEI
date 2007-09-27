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

#include <cmath>

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
     *     Product(a, b): \quad c \leftarrow a * b,
     * \f]
     * which yields c, the product of entities a and b.
     *
     * The referenced containers are invariant under this operation.
     * In every case, a new object is created and returned.
     *
     * \ingroup grpmatrixoperations
     * \ingroup grpvectoroperations
     */
    template <typename Tag_ = tags::CPU> struct Product
    {
        /**
         * \brief Returns the the product of two given entities.
         *
         * \param a The entity that is the first factor of the operation.
         * \param b The entity that is the second factor of the operation.
         *
         * \retval c Will create a new entity with Datatype of the first factor and return it.
         *
         * \exception MatrixSizeDoesNotMatch is thrown if two banded matrices do not have the same size.
         * \exception MatrixRowsDoNotMatch is thrown if two matrices do not have the same number of rows.
         * \exception MatrixColumnsDoNotMatch is thrown if two matrices do not have the same number of columns.
         * \exception MatrixIsNotSquare is thrown if a row access matrix's number of rows does not equal its number of columns.
         * \exception VectorSizeDoesNotMatch is thrown if two vectors do not have the same size.
         */

        /// \{

        template <typename DT1_, typename DT2_>
        static DenseVector<DT1_> value(const DenseMatrix<DT1_> & a, const DenseVector<DT2_> & b)
        {
            CONTEXT("When multiplying DenseMatrix with DenseVector:");

            if (b.size() != a.columns())
            {
                throw MatrixRowsDoNotMatch(a.columns(), b.size());
            }

            DenseVector<DT1_> result(a.rows());
            typename Vector<DT1_>::ElementIterator l(result.begin_elements());
            for (unsigned long i=0; i < a.rows(); ++i)
            {
                DenseVector<DT1_> dv(a[i]);
                *l = DotProduct<Tag_>::value(b, dv);
                ++l;
            }

            return result;
        }

        template <typename DT1_, typename DT2_>
        static DenseVector<DT1_> value(const DenseMatrix<DT1_> & a, const SparseVector<DT2_> & b)
        {
            CONTEXT("When multiplying DenseMatrix with SparseVector:");

            if (b.size() != a.columns())
            {
                throw MatrixRowsDoNotMatch(a.columns(), b.size());
            }

            DenseVector<DT1_> result(a.rows());
            typename Vector<DT1_>::ElementIterator l(result.begin_elements());
            for (unsigned long i=0; i < a.rows(); ++i)
            {
                DenseVector<DT1_> dv(a[i]);
                *l = DotProduct<Tag_>::value(b, dv);
                ++l;
            }

            return result;
        }

        template <typename DT1_, typename DT2_>
        static SparseVector<DT1_> value(const SparseMatrix<DT1_> & a, const DenseVector<DT2_> & b)
        {
            CONTEXT("When multiplying SparseMatrix with DenseVector:");

            if (b.size() != a.columns())
            {
                throw MatrixRowsDoNotMatch(a.columns(), b.size());
            }

            SparseVector<DT1_> result(a.rows(),1);
            typename Vector<DT1_>::ElementIterator l(result.begin_elements());
            for (unsigned long i=0; i < a.rows(); ++i)
            {
                SparseVector<DT1_>  v(a[i]);
                *l = DotProduct<Tag_>::value(v, b);
                ++l;
            }
            ///\todo: perhaps sparsify
            return result;
        }

        template <typename DT1_, typename DT2_>
        static SparseVector<DT1_> value(const SparseMatrix<DT1_> & a, const SparseVector<DT2_> & b)
        {
            CONTEXT("When multiplying SparseMatrix with SparseVector:");

            if (b.size() != a.columns())
            {
                throw MatrixRowsDoNotMatch(a.columns(), b.size());
            }

            SparseVector<DT1_> result(a.rows(), a.rows());
            typename SparseVector<DT1_>::ElementIterator l(result.begin_elements());
            for (unsigned long i=0; i < a.rows(); ++i)
            {
                const Vector<DT1_> & v(a[i]);
                *l = DotProduct<Tag_>::value(b, v);
                ++l;
            }
            ///\todo: perhaps sparsify
            return result;
        }

        template <typename DT1_, typename DT2_>
        static DenseVector<DT1_> value(const BandedMatrix<DT1_> & a, const DenseVector<DT2_> & b)
        {
            CONTEXT("When multiplying BandedMatrix with DenseVector:");

            if (b.size() != a.columns())
            {
                throw MatrixRowsDoNotMatch(a.columns(), b.size());
            }

            DenseVector<DT1_> result(a.rows(), DT1_(0));
            int middle_index(a.rows() -1);

            // If we are above or on the diagonal band, we start at Element 0 and go on until Element band_size-band_index.
            for (typename BandedMatrix<DT1_>::ConstVectorIterator vi(a.band_at(middle_index)), vi_end(a.end_bands()) ;
                    vi != vi_end ; ++vi)
            {
                typename Vector<DT2_>::ConstElementIterator j(b.begin_elements()), j_end(b.end_elements());
                typename Vector<DT1_>::ElementIterator r(result.begin_elements()), r_end(result.end_elements());

                for (unsigned int i(0) ; i < (vi.index()-middle_index) && j != j_end ; ++i)
                {
                    ++j; // Get the right position in b.
                }

                //Calculation of the element-index to stop in iteration!
                unsigned long end(vi->size() - (vi.index() - middle_index));
                for(typename Vector<DT1_>::ConstElementIterator c(vi->begin_elements()),
                        c_end(vi->element_at(end)) ; c != c_end ; ++c)
                {
                    *r += *c * *j;
                    if (j != j_end && r != r_end)
                    {
                        ++r;
                        ++j;
                    }
                }
            }

            // If we are below the diagonal band, we start at Element index and go on until the last element.
            for (typename BandedMatrix<DT1_>::ConstVectorIterator vi(a.begin_bands()),
                    vi_end(a.band_at(middle_index)) ; vi != vi_end ; ++vi)
            {
                typename Vector<DT2_>::ConstElementIterator j(b.begin_elements()), j_end(b.end_elements());
                typename Vector<DT1_>::ElementIterator r(result.begin_elements()), r_end(result.end_elements());
                unsigned long start(middle_index - vi.index()); //Calculation of the element-index to start in iteration!
                for (unsigned int a(0); a < middle_index - vi.index() && r != r_end; ++a)
                {
                    ++r; // Get the right position in result b.
                }
                for(typename Vector<DT1_>::ConstElementIterator c(vi->element_at(start)),
                        c_end(vi->end_elements()) ; c != c_end ; ++c)
                {
                    *r += *c * *j;
                    if (j != j_end && r != r_end)
                    {
                        ++j;
                        ++r;
                    }
                }
            }

            return result;
        }

        template <typename DT1_, typename DT2_>
        static SparseVector<DT1_> value(const BandedMatrix<DT1_> & a, const SparseVector<DT2_> & b)
        {
            CONTEXT("When multiplying BandedMatrix with SparseVector:");

            if (b.size() != a.columns())
            {
                throw MatrixRowsDoNotMatch(a.columns(), b.size());
            }

            SparseVector<DT1_> result(b.size(), b.capacity());
            unsigned long middle_index(a.rows() -1);

            for (typename BandedMatrix<DT1_>::ConstVectorIterator vi(a.begin_bands()),
                    vi_end(a.band_at(middle_index)) ; vi != vi_end; ++vi)
            {
                // If we are below the diagonal band we correct the position in band and result b.
                unsigned long move_index(middle_index - vi.index()); // is also = the element index to start in iteration for b
                for (typename Vector<DT2_>::ConstElementIterator j(b.begin_non_zero_elements()),
                        j_end(b.end_non_zero_elements()) ; j != j_end ; ++j)
                {
                    typename Vector<DT1_>::ConstElementIterator c(vi->element_at(move_index)), c_end(vi->end_elements());
                    typename Vector<DT1_>::ElementIterator r(result.element_at(move_index));
                    while (c.index() < (j.index() + move_index) && c != c_end) // Need a positive index here, so + is used!
                    {
                        ++c;
                        ++r;
                    }

                    if (c != c_end)
                    {
                        *r += *c * *j;
                    }
                }
            }

            for (typename BandedMatrix<DT1_>::ConstVectorIterator vi(a.band_at(middle_index)),
                    vi_end(a.end_bands()) ;  vi != vi_end ; ++vi)
            {
                // If we are above or on the diagonal band, we first search the first non-zero element, then
                // we make positions in band and result b meet it.
                unsigned long move_index(vi.index() - middle_index);
                unsigned long end(vi->size() - move_index); //Calculation of the element-index to stop in iteration!

                for (typename Vector<DT2_>::ConstElementIterator j(b.begin_non_zero_elements()),
                        j_end(b.end_non_zero_elements()) ; j != j_end ; ++j)
                {
                    if (j.index() < move_index)
                    {
                        continue;
                    }

                    typename Vector<DT1_>::ConstElementIterator c(vi->begin_elements()), c_end(vi->element_at(end));
                    typename Vector<DT1_>::ElementIterator r(result.begin_elements());

                    while (c.index() < (j.index() - move_index) && c != c_end)
                    {
                        ++c;
                        ++r;
                    }

                    if (c != c_end)
                    {
                        *r += *c * *j;
                    }
                }
            }

            ///\todo: perhaps sparsify (*c can be zero)
            return result;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> value(const DenseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
        {
            CONTEXT("When multiplying DenseMatrix with DenseMatrix:");

            if (a.columns() != b.rows())
                throw MatrixRowsDoNotMatch(b.rows(), a.columns());

            DenseMatrix<DT1_> result(b.columns(), a.rows());
            typename MutableMatrix<DT1_>::ElementIterator i(result.begin_elements());

            for (unsigned int s(0) ; s < a.rows() ; ++s)
            {
                const DenseVector<DT1_> a_row(a[s]);
                for (unsigned int t(0); t < b.columns() ; ++t)
                {
                    const DenseVector<DT2_> b_column(b.column(t));
                    *i = DotProduct<>::value(a_row, b_column);
                    ++i;
                }

            }
            return result;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> value(const DenseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
        {
            CONTEXT("When multiplying DenseMatrix with SparseMatrix:");

            if (a.columns() != b.rows())
                throw MatrixRowsDoNotMatch(b.rows(), a.columns());

            DenseMatrix<DT1_> result(b.columns(), a.rows(), DT1_(0));
            typename MutableMatrix<DT1_>::ElementIterator i(result.begin_elements());

            ///\todo: Should be optimized !!! (Use NonZeroIterators, less []-access ...)
            for (unsigned int l_row(0) ; l_row < a.rows() ; ++l_row)
            {
                const DenseVector<DT1_> a_row(a[l_row]);
                for (unsigned int r_column(0); r_column < b.columns() ; ++r_column)
                {
                    typename Vector<DT1_>::ConstElementIterator l(a_row.begin_elements());
                    for (unsigned int r_row(0); r_row < b.rows() ; ++r_row)
                    {
                        const SparseVector<DT2_> b_row(b[r_row]);
                        DT2_ b_value(b_row[r_column]);
                        //result[l_row][r_column] += b_value * *l;
                        *i += b_value * *l;
                        ++l;
                    }
                    ++i;
                }
            }
            return result;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> value(const SparseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
        {
            CONTEXT("When multiplying SparseMatrix with SparseMatrix:");

            if (a.columns() != b.rows())
                throw MatrixRowsDoNotMatch(b.rows(), a.columns());

            SparseMatrix<DT1_> result(b.columns(), a.rows());
            typename MutableMatrix<DT1_>::ElementIterator i(result.begin_elements());
            ///\todo: Should be optimized !!! (Use NonZeroIterators, less []-access ...)
            for (unsigned int l_row(0) ; l_row < a.rows() ; ++l_row)
            {
                const SparseVector<DT1_> a_row(a[l_row]);
                for (unsigned int r_column(0); r_column < b.columns() ; ++r_column)
                {
                    typename Vector<DT1_>::ConstElementIterator l(a_row.begin_elements());
                    for (unsigned int r_row(0); r_row < b.rows() ; ++r_row)
                    {
                        const SparseVector<DT2_> b_row(b[r_row]);
                        DT2_ b_value(b_row[r_column]);
                        //result[l_row][r_column] += b_value * *l;
                        *i += b_value * *l;
                        ++l;
                    }
                    ++i;
                }
            }
            return result;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> value(const SparseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
        {
            CONTEXT("when multiplying SparseMatrix with DenseMatrix:");

            if (a.columns() != b.rows())
                throw MatrixRowsDoNotMatch(b.rows(), a.columns());

            DenseMatrix<DT1_> result(b.columns(), a.rows());
            typename MutableMatrix<DT1_>::ElementIterator i(result.begin_elements());

            for (unsigned int s(0) ; s < a.rows() ; ++s)
            {
                const SparseVector<DT1_> a_row(a[s]);
                for (unsigned int t(0); t < b.columns() ; ++t)
                {
                    const DenseVector<DT2_> b_column(b.column(t));
                    *i = DotProduct<>::value(b_column, a_row);
                    ++i;
                }

            }
            return result;
        }

        template <typename DT1_, typename DT2_>
        static BandedMatrix<DT1_> value(const BandedMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
        {
            CONTEXT("When multiplying BandedMatrix with BandedMatrix:");

            if (a.size() != b.size())
                throw MatrixSizeDoesNotMatch(b.size(), a.size());

            BandedMatrix<DT1_> result(a.size());
            unsigned long diag_index(a.size() -1);

            // Lower part of a
            for (typename BandedMatrix<DT1_>::ConstVectorIterator vi(a.begin_bands()),
                    vi_end(a.band_at(diag_index)) ; vi != vi_end ; ++vi)
            {
                if (! vi.exists())
                    continue;

                // We start at the zero_based band_index of b that is equal to abs(band_index(a) - diag_index)
                    unsigned long iteration_start(abs(vi.index() - diag_index));
                    for(typename BandedMatrix<DT1_>::ConstVectorIterator vj(b.band_at(iteration_start)),
                            vj_end(b.end_bands()) ; vj != vj_end ; ++vj)
                    {
                        if (vj.index() == diag_index) // We are on diagonal of b
                        {
                            signed long result_band_index(vi.index() - diag_index); //index based on diag = 0
                            typename Vector<DT1_>::ConstElementIterator i(vi->element_at(abs(result_band_index)));
                            typename Vector<DT2_>::ConstElementIterator j(vj->begin_elements());
                            for (typename Vector<DT1_>::ElementIterator
                                    r(result.band(result_band_index).element_at(abs(result_band_index))),
                                    r_end(result.band(result_band_index).end_elements()) ; r != r_end ;
                                    ++j, ++i, ++r)
                            {
                                *r += *i * *j;
                            }
                        }
                        else if (vj.index() > diag_index) // We are above diagonal of b
                        {
                            signed long diag_based_index_a(vi.index() - diag_index);
                            signed long diag_based_index_b(vj.index() - diag_index);
                            signed long result_band_index(diag_based_index_a + diag_based_index_b);
                            unsigned long shift(diag_index - vi.index());
                            typename Vector<DT1_>::ConstElementIterator i(vi->element_at(shift));
                            typename Vector<DT2_>::ConstElementIterator j(vj->begin_elements());
                            long res_end(result.size());
                            if (result_band_index > 0)
                                res_end -= result_band_index;

                            for(typename Vector<DT1_>::ElementIterator r(result.band(result_band_index).element_at(shift)), r_end(result.band(result_band_index).element_at(res_end)) ; r != r_end ; ++j, ++i, ++r)
                            {
                                *r += *i * *j;
                            }
                        }

                        else if (vj.index() < diag_index) // We are below diagonal of b
                        {
                            signed long diag_based_index_a(vi.index() - diag_index);
                            signed long diag_based_index_b(vj.index() - diag_index);
                            signed long result_band_index(diag_based_index_a + diag_based_index_b);
                            unsigned long shift(diag_index - vj.index());
                            typename Vector<DT1_>::ConstElementIterator i(vi->element_at(abs(result_band_index)));
                            typename Vector<DT2_>::ConstElementIterator j(vj->element_at(shift));
                            for(typename Vector<DT1_>::ElementIterator r(result.band(result_band_index).element_at(abs(result_band_index))), r_end(result.band(result_band_index).end_elements()) ; r != r_end ; ++j, ++i, ++r)
                            {
                                *r += *i * *j;
                            }
                        }
                    }
                }

            // Diagonal of a
            typename BandedMatrix<DT1_>::ConstVectorIterator vi(a.band_at(diag_index));

            if (vi.exists())
            {
                // Lower part of b
                for (typename BandedMatrix<DT2_>::ConstVectorIterator vj(b.begin_bands()),
                        vj_end(b.band_at(diag_index)) ; vj != vj_end ; ++vj)
                {
                    if (! vj.exists())
                        continue;

                    signed long result_band_index(vj.index() - diag_index);  //index based on diag = 0
                    typename Vector<DT1_>::ElementIterator r(result.band(result_band_index).element_at(abs(result_band_index)));
                    typename Vector<DT1_>::ConstElementIterator i(vi->element_at(abs(result_band_index)));
                    for(typename Vector<DT2_>::ConstElementIterator j(vj->element_at(abs(result_band_index))),
                            j_end(vj->end_elements()) ; j != j_end ; ++j, ++i, ++r)
                    {
                        *r += *i * *j;
                    }
                }

                // Diagonal of b
                for(typename BandedMatrix<DT2_>::ConstVectorIterator vj(b.band_at(diag_index)),
                        vj_end(b.band_at(diag_index + 1)) ; vj != vj_end ; ++vj)
                {
                        if (! vj.exists())
                            continue;

                        DenseVector<DT1_> temp(vi->copy());
                        result.band(0) = Sum<>::value(ElementProduct<>::value(temp, *vj), result.band(0));
                }

                // Upper part of b
                for(typename BandedMatrix<DT2_>::ConstVectorIterator vj(b.band_at(diag_index+1)),
                        vj_end(b.end_bands()) ; vj != vj_end ; ++vj)
                {
                    signed long result_band_index(vj.index() - diag_index); //index based on diag = 0
                    typename Vector<DT1_>::ElementIterator r(result.band(result_band_index).begin_elements());
                    typename Vector<DT1_>::ConstElementIterator i(vi->begin_elements());
                    for(typename Vector<DT2_>::ConstElementIterator j(vj->begin_elements()),
                            j_end(vj->element_at(vj->size() - result_band_index)) ; j != j_end ; ++j, ++i, ++r)
                    {
                        *r += *i * *j;
                    }
                }
            }

            // Upper part of a
            for(typename BandedMatrix<DT1_>::ConstVectorIterator vi(a.band_at(diag_index+1)),
                    vi_end(a.end_bands()) ; vi != vi_end ; ++vi)
            {
                if (! vi.exists())
                    continue;

                    // We only go on until band_index of b is equal to [numbers of bands] - actual band_index of a (zero_based).
                    unsigned long iteration_end((2*b.size()-1) - (vi.index() - diag_index));
                    for(typename BandedMatrix<DT1_>::ConstVectorIterator vj(b.begin_bands()), vj_end(b.band_at(iteration_end)) ; vj != vj_end ; ++vj)
                    {
                        if (vj.index() == diag_index) // We are on diagonal of b
                        {
                            signed long result_band_index(vi.index() - diag_index); //index based on diag = 0
                            typename Vector<DT1_>::ConstElementIterator i(vi->begin_elements());
                            typename Vector<DT2_>::ConstElementIterator j(vj->element_at(result_band_index));
                            for(typename Vector<DT1_>::ElementIterator r(result.band(result_band_index).begin_elements()),
                                    r_end(result.band(result_band_index).element_at(result.size() - result_band_index)) ; r != r_end ; ++j, ++i, ++r)
                            {
                                *r += *i * *j;
                            }
                        }
                        else if (vj.index() > diag_index) // We are above diagonal of b
                        {
                            signed long diag_based_index_a(vi.index() - diag_index);
                            signed long diag_based_index_b(vj.index() - diag_index);
                            signed long result_band_index(diag_based_index_a + diag_based_index_b);
                            typename Vector<DT1_>::ConstElementIterator i(vi->begin_elements());
                            unsigned long shift(vi.index() - diag_index);
                            typename Vector<DT2_>::ConstElementIterator j(vj->element_at(shift));
                            for(typename Vector<DT1_>::ElementIterator r(result.band(result_band_index).begin_elements()),
                                    r_end(result.band(result_band_index).element_at(result.size() - result_band_index)) ; r != r_end ; ++j, ++i, ++r)
                            {
                                *r += *i * *j;
                            }
                        }

                        else if (vj.index() < diag_index) // We are below diagonal of b
                        {
                            signed long diag_based_index_a(vi.index() - diag_index);
                            signed long diag_based_index_b(vj.index() - diag_index);
                            signed long result_band_index(diag_based_index_a + diag_based_index_b);

                            long res_start(0);
                            if (result_band_index < 0)
                                res_start = abs(result_band_index);

                            typename Vector<DT1_>::ConstElementIterator i(vi->element_at(res_start));
                            long vj_start(vi.index() - diag_index);
                            if (result_band_index < 0)
                                vj_start += abs(result_band_index);

                            typename Vector<DT2_>::ConstElementIterator j(vj->element_at(vj_start));
                            long res_end(2*result.size() - (vi.index() + 1));
                            for(typename Vector<DT1_>::ElementIterator r(result.band(result_band_index).element_at(res_start)),
                                    r_end(result.band(result_band_index).element_at(res_end)) ; r != r_end ; ++j, ++i, ++r)
                            {
                                *r += *i * *j;
                            }
                        }
                    }
            }

            return result;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT2_> value(const BandedMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
        {
            CONTEXT("When multiplying BandedMatrix with DenseMatrix:");

            if (a.columns() != b.rows())
                throw MatrixRowsDoNotMatch(b.rows(), a.columns());

            if (b.columns() != a.rows())
                throw MatrixRowsDoNotMatch(a.rows(), b.columns());

            DenseMatrix<DT2_> result(b.columns(), b.rows(), DT2_(0));
            unsigned long middle_index(a.size() -1);

            // Calculation for lower part
            for (typename BandedMatrix<DT1_>::ConstVectorIterator vi(a.begin_bands()),
                    vi_end(a.band_at(middle_index)) ; vi != vi_end ; ++vi)
            {
                if (! vi.exists())
                    continue;

                for (unsigned int s(0) ; s < b.columns() ; ++s)
                {
                    // Temporary container for efficient calculation of elementwise vector product.
                    DenseVector<DT2_> temp(b.rows(), DT2_(0));
                    unsigned long real_index(middle_index - vi.index());
                    typename Vector<DT2_>::ConstElementIterator c(b.column(s).begin_elements()),
                                d(vi->element_at(real_index));
                    for (typename Vector<DT2_>::ElementIterator x(temp.element_at(real_index)),
                                x_end(temp.end_elements()) ; x != x_end ; ++x)
                    {
                        *x = *c * *d;
                        ++c; ++d;
                    }

                    Sum<>::value(result.column(s), temp);
                }
            }

            // Calculation for diagonal
            typename BandedMatrix<DT1_>::ConstVectorIterator vi(a.band_at(middle_index));

            if (vi.exists())
            {
                for (unsigned int s(0) ; s < b.columns() ; ++s)
                {
                    DenseVector<DT2_> temp(b.column(s).copy());
                    Sum<>::value(result.column(s), ElementProduct<>::value(temp, *vi));
                }
            }

            // Calculation for upper part
            for (typename BandedMatrix<DT1_>::ConstVectorIterator vi(a.band_at(middle_index+1)),
                    vi_end(a.end_bands()) ; vi != vi_end ; ++vi)
            {
                if (! vi.exists())
                    continue;

                for (unsigned int s(0) ; s < b.columns() ; ++s)
                {
                    // Temporary container for efficient calculation of elementwise vector product.
                    DenseVector<DT2_> temp(b.rows(), DT2_(0));

                    unsigned long real_index(vi.index() - middle_index);
                    typename Vector<DT2_>::ConstElementIterator c(b.column(s).element_at(real_index)),
                            d(vi->begin_elements());
                    unsigned long end(temp.size() - real_index);

                    for (typename Vector<DT2_>::ElementIterator x(temp.begin_elements()),
                            x_end(temp.element_at(end)) ; x != x_end ; ++x)
                    {
                        *x = *c * *d;
                        ++c; ++d;
                    }
                        Sum<>::value(result.column(s), temp);
                }
            }

            return result;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> value(const BandedMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
        {
            CONTEXT("When multiplying BandedMatrix with SparseMatrix:");

            if (a.columns() != b.rows())
                throw MatrixRowsDoNotMatch(b.rows(), a.columns());

            if (b.columns() != a.rows())
                throw MatrixRowsDoNotMatch(a.rows(), b.columns());

            DenseMatrix<DT2_> result(a.columns(), a.rows(), DT2_(0));
            unsigned long middle_index(a.size() -1);

            // Calculation for lower part
            for (typename BandedMatrix<DT1_>::ConstVectorIterator vi(a.begin_bands()),
                    vi_end(a.band_at(middle_index)) ; vi != vi_end ; ++vi)
            {
                if (! vi.exists())
                    continue;

                unsigned long real_index(middle_index - vi.index());

                typename Vector<DT1_>::ConstElementIterator d(vi->element_at(real_index)),
                        d_end(vi->end_elements());
                for (unsigned int z(0) ; z < (b.rows()-real_index) ; ++z, ++d)
                {
                    const SparseVector<DT2_> row(b[z]);
                    typename Vector<DT2_>::ElementIterator x(result[d.index()].begin_elements());
                    for(typename Vector<DT2_>::ConstElementIterator c(row.begin_elements()),
                            c_end(row.end_elements()) ; c != c_end ; ++x, ++c)
                    {
                        *x += *d * *c;
                    }
                }
            }

            // Calculation for diagonal
            typename BandedMatrix<DT1_>::ConstVectorIterator vi(a.band_at(middle_index));

            if (vi.exists())
            {
                typename Vector<DT1_>::ConstElementIterator d(vi->begin_elements());
                for (unsigned int z(0) ; z < b.rows() ; ++z, ++d)
                {
                    const SparseVector<DT2_> row(b[z]);
                    typename Vector<DT2_>::ElementIterator x(result[z].begin_elements());
                    for(typename Vector<DT2_>::ConstElementIterator c(row.begin_elements()),
                            c_end(row.end_elements()) ; c != c_end ; ++x, ++c)
                    {
                        *x += *d * *c;
                    }
                }
            }

            // Calculation for upper part
            for (typename BandedMatrix<DT1_>::ConstVectorIterator vi(a.band_at(middle_index+1)),
                    vi_end(a.end_bands()) ; vi != vi_end ; ++vi)
            {
                if (! vi.exists())
                    continue;

                unsigned long real_index(vi.index() - middle_index);

                typename Vector<DT1_>::ConstElementIterator d(vi->begin_elements()),
                        d_end(vi->element_at(vi->size()-real_index));
                for (unsigned int z(real_index) ; z < b.rows() ; ++z, ++d)
                {
                    const SparseVector<DT2_> row(b[z]);
                    typename Vector<DT2_>::ElementIterator x(result[d.index()].begin_elements());
                    for(typename Vector<DT2_>::ConstElementIterator c(row.begin_elements()), c_end(row.end_elements()) ;
                            c != c_end ; ++x, ++c)
                    {
                        *x += *d * *c;
                    }
                }
            }

            return result;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> value(const DenseMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
        {
            CONTEXT("When multiplying DenseMatrix with BandedMatrix:");

            if (a.columns() != b.rows())
                throw MatrixRowsDoNotMatch(b.rows(), a.columns());

            DenseMatrix<DT1_> result(b.columns(), a.rows(), DT1_(0));

            unsigned long middle_index(b.size() -1);

            // Calculation for lower part
            for (typename BandedMatrix<DT2_>::ConstVectorIterator vi(b.begin_bands()),
                    vi_end(b.band_at(middle_index)) ; vi != vi_end ; ++vi)
            {
                if (! vi.exists())
                    continue;

                unsigned long real_index(vi.index() - middle_index);

                for (unsigned int z(0) ; z < a.rows() ; ++z)
                {
                    const DenseVector<DT1_> row(a[z]);
                    typename Vector<DT1_>::ConstElementIterator c(row.begin_elements());
                    typename Vector<DT1_>::ElementIterator x(result[z].element_at(real_index));
                    for(typename Vector<DT2_>::ConstElementIterator d(vi->begin_elements()),
                        d_end(vi->element_at(vi->size()-real_index)) ; d != d_end ; ++x, ++c, ++d)
                    {
                        *x += *d * *c;
                    }
                }
            }

            // Calculation for diagonal
            typename BandedMatrix<DT2_>::ConstVectorIterator vi(b.band_at(middle_index));

            if (vi.exists())
            {
                for (unsigned int z(0) ; z < a.rows() ; ++z)
                {
                    const DenseVector<DT1_> row(a[z]);
                    typename Vector<DT1_>::ConstElementIterator c(row.begin_elements());
                    typename Vector<DT1_>::ElementIterator x(result[z].begin_elements());
                    for(typename Vector<DT2_>::ConstElementIterator d(vi->begin_elements()),
                        d_end(vi->end_elements()) ; d != d_end ; ++x, ++c, ++d)
                    {
                        *x += *d * *c;
                    }
                }
            }

            // Calculation for upper part
            for (typename BandedMatrix<DT2_>::ConstVectorIterator vi(b.band_at(middle_index+1)),
                    vi_end(b.end_bands()) ; vi != vi_end ; ++vi)
            {
                if (! vi.exists())
                    continue;

                unsigned long real_index(middle_index - vi.index());

                for (unsigned int z(0) ; z < a.rows() ; ++z)
                {
                    const DenseVector<DT1_> row(a[z]);
                    typename Vector<DT1_>::ConstElementIterator c(row.element_at(real_index));
                    typename Vector<DT1_>::ElementIterator x(result[z].begin_elements());
                    for(typename Vector<DT2_>::ConstElementIterator d(vi->element_at(real_index)),
                        d_end(vi->end_elements()) ; d != d_end ; ++x, ++c, ++d)
                    {
                        *x += *d * *c;
                    }
                }
            }

            return result;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> value(const SparseMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
        {
            CONTEXT("When multiplying SparseMatrix with BandedMatrix:");

            if (a.columns() != b.rows())
                throw MatrixRowsDoNotMatch(b.rows(), a.columns());

            DenseMatrix<DT1_> result(b.columns(), a.rows(), DT1_(0));

            unsigned long middle_index(b.size() -1);

            // Calculation for lower part
            for (typename BandedMatrix<DT2_>::ConstVectorIterator vi(b.begin_bands()),
                    vi_end(b.band_at(middle_index)) ; vi != vi_end ; ++vi)
            {
                if (! vi.exists())
                    continue;

                unsigned long real_index(middle_index - vi.index());

                for (unsigned int z(0) ; z < a.rows() ; ++z)
                {
                    const SparseVector<DT1_> row(a[z]);
                    typename Vector<DT1_>::ConstElementIterator c(row.element_at(real_index));
                    typename Vector<DT1_>::ElementIterator x(result[z].begin_elements());
                    for(typename Vector<DT2_>::ConstElementIterator d(vi->element_at(real_index)),
                        d_end(vi->end_elements()) ; d != d_end ; ++x, ++c, ++d)
                    {
                        *x += *d * *c;
                    }
                }
            }
            // Calculation for diagonal part
            typename BandedMatrix<DT2_>::ConstVectorIterator vi(b.band_at(middle_index));

            if (vi.exists())
            {
                for (unsigned int z(0) ; z < a.rows() ; ++z)
                {
                    const SparseVector<DT1_> row(a[z]);
                    typename Vector<DT1_>::ConstElementIterator c(row.begin_elements());

                    typename Vector<DT1_>::ElementIterator x(result[z].begin_elements());
                    for(typename Vector<DT2_>::ConstElementIterator d(vi->begin_elements()),
                        d_end(vi->end_elements()) ; d != d_end ; ++x, ++c, ++d)
                    {
                        *x += *d * *c;
                    }
                }
            }

            // Calculation for upper part
            for (typename BandedMatrix<DT2_>::ConstVectorIterator vi(b.band_at(middle_index+1)),
                    vi_end(b.end_bands()) ; vi != vi_end ; ++vi)
            {
                if (! vi.exists())
                    continue;

                const unsigned long real_index(vi.index() - middle_index);

                for (unsigned int z(0) ; z < a.rows() ; ++z)
                {
                    const SparseVector<DT1_> row(a[z]);
                    typename Vector<DT1_>::ConstElementIterator c(row.begin_elements());
                    typename Vector<DT1_>::ElementIterator x(result[z].element_at(real_index));
                    for(typename Vector<DT2_>::ConstElementIterator d(vi->begin_elements()),
                        d_end(vi->element_at(vi->size()-real_index)) ; d != d_end ; ++x, ++c, ++d)
                    {
                        *x += *d * *c;
                    }
                }
            }

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

                for (unsigned int s(0) ; s < left.rows() ; ++s)
                {
                    const DenseVector<DT1_> left_row(left[s]);
                    for (unsigned int t=0; t < right.columns() ; ++t)
                    {
                        const DenseVector<DT2_> right_column(right.column(t));
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

   /// SSE implementiation.
    template <>
    struct Product<tags::CPU::SSE>
    {
        static DenseVector<float> value(const BandedMatrix<float> & a, const DenseVector<float> & b);
    };
}
#endif
