/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Sven Mallach <sven.mallach@honei.org>
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

#include <honei/libla/banded_matrix.hh>
#include <honei/libla/dense_vector.hh>
#include <honei/libla/dense_matrix.hh>
#include <honei/libla/dense_matrix_tile.hh>
#include <honei/libla/dot_product.hh>
#include <honei/libla/element_product.hh>
#include <honei/libla/matrix_error.hh>
#include <honei/libla/product-mc.hh>
#include <honei/libla/scale.hh>
#include <honei/libla/scaled_sum.hh>
#include <honei/libla/sparse_matrix.hh>
#include <honei/libla/sparse_vector.hh>
#include <honei/libla/sum.hh>
#include <honei/libutil/tags.hh>
#include <honei/libutil/benchmark_info.hh>

#include <cmath>

namespace honei
{
    /**
     * \brief Product of two entities.
     *
     * MatrixProduct is the class template for the product operations
     * \f[
     *     \texttt{Product}(a, b): \quad c \leftarrow a * b
     * \f]
     * which yields c, the product of entities a and b.
     *
     * The referenced containers are invariant under this operation.
     * In every case, a new object is created and returned.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <typename Tag_ = tags::CPU> struct Product
    {
        /**
         * \name Right-handed Matrix-Vector products
         * \{
         *
         * \brief Returns the right-handed Matrix-Vector product.
         *
         * \param a The matrix that is the left-hand factor of the operation.
         * \param b The vector that is the right-hand factor of the operation.
         *
         * \retval c Will create a new vector with Datatype of the first factor and return it.
         *
         * \exception VectorSizeDoesNotMatch is thrown if the vector's size does not match the matrix's number
         *            of columns.
         */

        template <typename DT1_, typename DT2_>
        static DenseVector<DT1_> value(const DenseMatrix<DT1_> & a, const DenseVectorBase<DT2_> & b)
        {
            CONTEXT("When multiplying DenseMatrix with DenseVector(Base):");

            if (b.size() != a.columns())
            {
                throw VectorSizeDoesNotMatch(b.size(), a.columns());
            }

            DenseVector<DT1_> result(a.rows());
            for (typename Vector<DT1_>::ElementIterator r(result.begin_elements()), r_end(result.end_elements()) ;
                    r != r_end ; ++r)
            {
                const DenseVectorRange<DT1_> dv(a[r.index()]);
                *r = DotProduct<Tag_>::value(dv, b);
            }

            return result;
        }

        template <typename DT1_, typename DT2_>
        static DenseVector<DT1_> value(const DenseMatrix<DT1_> & a, const SparseVector<DT2_> & b)
        {
            CONTEXT("When multiplying DenseMatrix with SparseVector:");

            if (b.size() != a.columns())
            {
                throw VectorSizeDoesNotMatch(b.size(), a.columns());
            }

            DenseVector<DT1_> result(a.rows());
            for (typename Vector<DT1_>::ElementIterator r(result.begin_elements()), r_end(result.end_elements()) ;
                    r != r_end ; ++r)
            {
                const DenseVectorRange<DT1_> dv(a[r.index()]);
                *r = DotProduct<Tag_>::value(dv, b);
            }

            return result;
        }

        template <typename DT1_, typename DT2_>
        static SparseVector<DT1_> value(const SparseMatrix<DT1_> & a, const DenseVectorBase<DT2_> & b)
        {
            CONTEXT("When multiplying SparseMatrix with DenseVectorBase:");

            if (b.size() != a.columns())
            {
                throw VectorSizeDoesNotMatch(b.size(), a.columns());
            }

            SparseVector<DT1_> result(a.rows(),1);
            for (typename Vector<DT1_>::ElementIterator r(result.begin_elements()), r_end(result.end_elements()) ;
                    r != r_end ; ++r)
            {
                const SparseVector<DT1_> sv(a[r.index()]);
                *r = DotProduct<Tag_>::value(sv, b);
            }

            return result;
        }

        template <typename DT1_, typename DT2_>
        static SparseVector<DT1_> value(const SparseMatrix<DT1_> & a, const SparseVector<DT2_> & b)
        {
            CONTEXT("When multiplying SparseMatrix with SparseVector:");

            if (b.size() != a.columns())
            {
                throw VectorSizeDoesNotMatch(b.size(), a.columns());
            }

            SparseVector<DT1_> result(a.rows(), a.rows());
            for (typename SparseVector<DT1_>::ElementIterator r(result.begin_elements()), r_end(result.end_elements()) ;
                    r != r_end ; ++r)
            {
                const SparseVector<DT1_> sv(a[r.index()]);
                *r = DotProduct<Tag_>::value(sv, b);
            }

            return result;
        }

        template <typename DT1_, typename DT2_>
        static DenseVector<DT1_> value(const BandedMatrix<DT1_> & a, const DenseVectorBase<DT2_> & b)
        {
            CONTEXT("When multiplying BandedMatrix with DenseVector(Base):");

            if (b.size() != a.columns())
            {
                throw VectorSizeDoesNotMatch(b.size(), a.columns());
            }

            DenseVector<DT1_> result(a.rows(), DT1_(0));
            int middle_index(a.rows() -1);
            for (typename BandedMatrix<DT1_>::ConstVectorIterator vi(a.begin_non_zero_bands()), vi_end(a.end_non_zero_bands()) ;
                    vi != vi_end ; ++vi)
            {
                // If we are above or on the diagonal band, we start at Element 0 and go on until Element band_size-band_index.
                if (middle_index < vi.index())
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
                else
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
            }
                return result;
        }

        template <typename DT1_, typename DT2_>
        static SparseVector<DT1_> value(const BandedMatrix<DT1_> & a, const SparseVector<DT2_> & b)
        {
            CONTEXT("When multiplying BandedMatrix with SparseVector:");

            if (b.size() != a.columns())
            {
                throw VectorSizeDoesNotMatch(b.size(), a.columns());
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
                    unsigned long ctr(0);
                    while (c.index() < (j.index() + move_index) && c != c_end) // Need a positive index here, so + is used!
                    {
                        ++c;
                        ctr++;
                    }

                    if (c != c_end)
                    {
                        result[move_index + ctr] += *c * *j;
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

        /// \}

        /**
         * \name Left-handed Matrix-Vector products
         * \{
         *
         * \brief Returns the left-handed Matrix-Vector product.
         *
         * \param b The vector that is the left-hand factor of the operation.
         * \param a The matrix that is the right-hand factor of the operation.
         *
         * \retval c Will create a new vector with Datatype of the first factor and return it.
         *
         * \exception MatrixRowsDoNotMatch is thrown if the vector's size does not match the matrix's number
         *            of rows.
         */

        template <typename DT1_, typename DT2_>
        static DenseVector<DT1_> value(const DenseVectorBase<DT1_> & a, const DenseMatrix<DT1_> & b)
        {
            CONTEXT("When multiplying DenseVector(Base) with DenseMatrix:");

            if (b.rows() != a.size())
            {
                throw MatrixRowsDoNotMatch(a.size(), b.rows());
            }

            DenseVector<DT1_> result(b[0].copy());
            Scale<tags::CPU>::value(result, a[0]);

            for (typename Vector<DT1_>::ConstElementIterator i(a.element_at(1)), i_end(a.end_elements()) ;
                    i != i_end ; ++i)
            {
                ScaledSum<tags::CPU>::value(result, b[i.index()], *i);
            }

            return result;
        }

        template <typename DT1_, typename DT2_>
        static DenseVector<DT1_> value(const DenseVectorBase<DT1_> & a, const DenseMatrixTile<DT2_> & b)
        {
            CONTEXT("When multiplying DenseVector(Base) with DenseMatrix:");

            if (b.rows() != a.size())
            {
                throw MatrixRowsDoNotMatch(a.size(), b.rows());
            }

            DenseVector<DT1_> result(b[0].copy());
            Scale<tags::CPU>::value(result, a[0]);

            for (typename Vector<DT1_>::ConstElementIterator i(a.element_at(1)), i_end(a.end_elements()) ;
                    i != i_end ; ++i)
            {
                ScaledSum<tags::CPU>::value(result, b[i.index()], *i);
            }

            return result;
        }

        /// \}

        /**
         * \name Matrix-Matrix products
         * \{
         *
         * \brief Returns the the product of two \link<Matrix>matrices\endlink.
         *
         * \param a The matrix that is the first factor of the operation.
         * \param b The vector that is the second factor of the operation.
         *
         * \retval c Will create a new matrix with Datatype of the first factor and return it.
         *
         * \exception MatrixSizeDoesNotMatch is thrown if two banded matrices do not have the same size.
         * \exception MatrixRowsDoNotMatch is thrown if the second matrix's number of rows does not match the
         *            first matrix's number of columns.
         */

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> value(const DenseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
        {
            CONTEXT("When multiplying DenseMatrix with DenseMatrix:");

            if (a.columns() != b.rows())
                throw MatrixRowsDoNotMatch(b.rows(), a.columns());

            DenseMatrix<DT1_> result(a.rows(), b.columns());
            typename MutableMatrix<DT1_>::ElementIterator i(result.begin_elements());

            for (unsigned int s(0) ; s < a.rows() ; ++s)
            {
                const DenseVectorRange<DT1_> a_row(a[s]);
                for (unsigned int t(0); t < b.columns() ; ++t)
                {
                    const DenseVectorSlice<DT2_> b_column(b.column(t));
                    *i = DotProduct<>::value(b_column, a_row);
                    ++i;
                }

            }
            return result;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrixTile<DT1_> & value(DenseMatrixTile<DT1_> & r, const DenseMatrixTile<DT1_> & a, const DenseMatrixTile<DT2_> & b)
        {
            CONTEXT("When multiplying DenseMatrixTile with DenseMatrixTile:");

            if (a.columns() != b.rows())
                throw MatrixRowsDoNotMatch(b.rows(), a.columns());

            typename MutableMatrix<DT1_>::ElementIterator i(r.begin_elements());

            for (unsigned long j(0) ; j < a.rows() ; ++j)
            {
                for (unsigned long k(0) ; k < b.columns() ; ++k)
                {
                    *i += DotProduct<>::value(a[j], b.column(k));
                    ++i;
                }
            }

            return r;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> value(const DenseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
        {
            CONTEXT("When multiplying DenseMatrix with SparseMatrix:");

            if (a.columns() != b.rows())
                throw MatrixRowsDoNotMatch(b.rows(), a.columns());

            DenseMatrix<DT1_> result(a.rows(), b.columns(), DT1_(0));
            for( typename Matrix<DT2_>::ConstElementIterator i(b.begin_non_zero_elements()), i_end(b.end_non_zero_elements()) ;
                    i < i_end ; ++i )
            {
                ScaledSum<>::value(result.column(i.column()), a.column(i.row()), *i);
            }

            return result;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> value(const SparseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
        {
            CONTEXT("When multiplying SparseMatrix with SparseMatrix:");

            if (a.columns() != b.rows())
                throw MatrixRowsDoNotMatch(b.rows(), a.columns());

            SparseMatrix<DT1_> result(a.rows(), b.columns());
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
            CONTEXT("When multiplying SparseMatrix with DenseMatrix:");

            if (a.columns() != b.rows())
                throw MatrixRowsDoNotMatch(b.rows(), a.columns());

            DenseMatrix<DT1_> result(a.rows(), b.columns(), DT1_(0));
            for( typename Matrix<DT1_>::ConstElementIterator i(a.begin_non_zero_elements()), i_end(a.end_non_zero_elements()) ;
                    i < i_end ; ++i )
            {
                ScaledSum<>::value(result[i.row()], b[i.column()], *i);
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

                            for(typename Vector<DT1_>::ElementIterator r(result.band(result_band_index).element_at(shift)),
                                    r_end(result.band(result_band_index).element_at(res_end)) ; r != r_end ; ++j, ++i, ++r)
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
                            for(typename Vector<DT1_>::ElementIterator r(result.band(result_band_index).element_at(abs(result_band_index))), 
                                    r_end(result.band(result_band_index).end_elements()) ; r != r_end ; ++j, ++i, ++r)
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
                        Sum<>::value(result.band(0), ElementProduct<>::value(temp, *vj));
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

            DenseMatrix<DT2_> result(b.rows(), b.columns(), DT2_(0));
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

            DenseMatrix<DT2_> result(a.rows(), b.columns(), DT2_(0));
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

            DenseMatrix<DT1_> result(a.rows(), b.columns(), DT1_(0));

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
                    /// \todo RowIterator
                    const DenseVectorRange<DT1_> row(a[z]);
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
                    const DenseVectorRange<DT1_> row(a[z]);
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
                    const DenseVectorRange<DT1_> row(a[z]);
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

            DenseMatrix<DT1_> result(a.rows(), b.columns(), DT1_(0));

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
                    typename Vector<DT1_>::ElementIterator x(result[z].begin_elements());
                    for(typename Vector<DT2_>::ConstElementIterator d(vi->element_at(real_index)),
                        d_end(vi->end_elements()) ; d != d_end ; ++x, ++d)
                    {
                        *x += *d * row[real_index + x.index()];
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

        /// \}

        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(BandedMatrix<DT1_> & a, DenseVectorBase<DT2_> & b)
        {
            BenchmarkInfo result;
            unsigned long temp(0);
            for (typename BandedMatrix<DT1_>::ConstVectorIterator vi(a.begin_non_zero_bands()), vi_end(a.end_non_zero_bands()) ;
                    vi != vi_end ; ++vi)
            {
                temp = vi.index();
                if (temp > (a.size() - 1))
                    temp = ((2 * a.size() - 2) - temp);
                ++temp;
                result.flops += temp * 2;
                result.load += temp * (sizeof(DT1_) + sizeof(DT2_) + sizeof(DT1_));
                result.store += temp * sizeof(DT1_);
            }
            result.size.push_back(a.size());
            result.size.push_back(b.size());
            return result;
        }

        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(DenseMatrix<DT1_> & a, DenseMatrix<DT2_> & b)
        {
            BenchmarkInfo result;
            DenseVector<DT1_> temp1(a.columns());
            DenseVector<DT2_> temp2(b.rows());
            result = DotProduct<>::get_benchmark_info(temp1, temp2) * (a.rows() * b.columns());
            result.size.clear();
            result.size.push_back(a.rows() * a.columns());
            result.size.push_back(b.rows() * b.columns());
            return result;
        }

        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(DenseMatrix<DT1_> & a, DenseVectorBase<DT2_> & b)
        {
            BenchmarkInfo result;
            DenseVector<DT1_> temp(a.columns());
            result = DotProduct<>::get_benchmark_info(temp, b) * a.rows();
            result.size.clear();
            result.size.push_back(a.rows() * a.columns());
            result.size.push_back(b.size());
            return result;
        }

        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(SparseMatrix<DT1_> & a, DenseMatrix<DT2_> & b)
        {
            BenchmarkInfo result;
            for (unsigned int s(0) ; s < a.rows() ; ++s)
            {
                SparseVector<DT1_> a_row(a[s]);
                for (unsigned int t(0); t < b.columns() ; ++t)
                {
                    DenseVectorSlice<DT2_> b_column(b.column(t));
                    result = result + DotProduct<>::get_benchmark_info(a_row, b_column);
                }

            }
            result.size.push_back(a.rows() * a.columns());
            result.size.push_back(b.rows() * b.columns());
            BenchmarkInfo result2;
            DenseVector<DT1_> temp1(a.columns());
            DenseVector<DT2_> temp2(b.rows());
            result2 = DotProduct<>::get_benchmark_info(temp1, temp2) * (a.rows() * a.columns());
            result.scale = ((double)result2.flops / result.flops);
            return result;
        }

        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(BandedMatrix<DT1_> & a, DenseMatrix<DT2_> & b)
        {
            BenchmarkInfo result;
            DenseMatrix<DT2_> presult(b.rows(), b.columns(), DT2_(0));
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
                        //*x = *c * *d;
                        result.flops += 1;
                        result.load += sizeof(DT2_) + sizeof(DT1_);
                        result.store += sizeof(DT2_);
                        ++c; ++d;
                    }

                    //Sum<>::value(presult.column(s), temp);
                    result = result + Sum<>::get_benchmark_info(presult.column(s), temp);
                }
            }

            // Calculation for diagonal
            typename BandedMatrix<DT1_>::ConstVectorIterator vi(a.band_at(middle_index));

            if (vi.exists())
            {
                for (unsigned int s(0) ; s < b.columns() ; ++s)
                {
                    DenseVector<DT2_> temp(b.column(s).copy());
                    //Sum<>::value(presult.column(s), ElementProduct<>::value(temp, *vi));
                    result = result + ElementProduct<>::get_benchmark_info(temp, *vi);
                    result = result + Sum<>::get_benchmark_info(presult.column(s), ElementProduct<>::value(temp, *vi));
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
                        //*x = *c * *d;
                        result.flops += 1;
                        result.load += sizeof(DT2_) + sizeof(DT1_);
                        result.store += sizeof(DT2_);
                        ++c; ++d;
                    }
                    //Sum<>::value(presult.column(s), temp);
                    result = result + Sum<>::get_benchmark_info(presult.column(s), temp);
                }
            }
            result.size.push_back(a.rows() * a.columns());
            result.size.push_back(b.rows() * b.columns());
            DenseMatrix<DT1_> t1(a.size(), a.size());
            BenchmarkInfo tinfo(get_benchmark_info(t1, t1));
            result.scale = (double(tinfo.flops) / result.flops);
            return result;
        }
    };

    /**
     * \brief Product of two entities.
     *
     * MatrixProduct is the class template for the product operation
     * \f[
     *     \texttt{Product}(a, b): \quad c \leftarrow a * b,
     * \f]
     * which yields c, the product of entities a and b.
     *
     * The referenced containers are invariant under this operation.
     * In every case, a new object is created and returned.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <> struct Product<tags::CPU::SSE>
    {
        /**
         * \name Products
         * \{
         *
         * \brief Returns the product of a Matrix and a Vector.
         *
         * \param a The matrix that is the first factor of the operation.
         * \param b The vector that is the second factor of the operation.
         *
         * \retval c Will create a new entity with Datatype of the first factor and return it.
         *
         * \exception MatrixSizeDoesNotMatch is thrown if two banded matrices do not have the same size.
         * \exception MatrixRowsDoNotMatch is thrown if two matrices do not have the same number of rows.
         * \exception MatrixColumnsDoNotMatch is thrown if two matrices do not have the same number of columns.
         * \exception MatrixIsNotSquare is thrown if a row access matrix's number of rows does not equal its number of columns.
         * \exception VectorSizeDoesNotMatch is thrown if two vectors do not have the same size.
         */

        static DenseVector<float> value(const BandedMatrix<float> & a, const DenseVectorContinuousBase<float> & b);

        static DenseVector<double> value(const BandedMatrix<double> & a, const DenseVectorContinuousBase<double> & b);

        static DenseVector<float> value(const DenseMatrix<float> & a, const DenseVectorContinuousBase<float> & b);

        static DenseVector<double> value(const DenseMatrix<double> & a, const DenseVectorContinuousBase<double> & b);

        static DenseMatrix<float> value(const DenseMatrix<float> &a, const DenseMatrix<float> & b);

        static DenseMatrix<double> value(const DenseMatrix<double> &a, const DenseMatrix<double> & b);

        static DenseMatrix<float> value(const SparseMatrix<float> &a, const DenseMatrix<float> & b);

        static DenseMatrix<double> value(const SparseMatrix<double> &a, const DenseMatrix<double> & b);

        static DenseMatrixTile<float> & value(DenseMatrixTile<float> & r, const DenseMatrixTile<float> & a, const DenseMatrixTile<float> & b);

        static DenseMatrixTile<double> & value(DenseMatrixTile<double> & r, const DenseMatrixTile<double> & a, const DenseMatrixTile<double> & b);

        /// \}
    };

    template <>
    struct Product<tags::Cell>
    {
        static DenseVector<float> value(const DenseMatrix<float> & a, const DenseVector<float> & b);
        static DenseMatrix<float> value(const DenseMatrix<float> & a, const DenseMatrix<float> & b);
        static DenseVector<float> value(const BandedMatrix<float> & a, const DenseVector<float> & b);
        static DenseVector<double> value(const BandedMatrix<double> & a, const DenseVector<double> & b);
        static DenseMatrix<float> value(const SparseMatrix<float> & a, const DenseMatrix<float> & b);
    };

    /**
     * \brief Product of two entities.
     *
     * MatrixProduct is the class template for the product operation
     * \f[
     *     \texttt{Product}(a, b): \quad c \leftarrow a * b,
     * \f]
     * which yields c, the product of entities a and b.
     *
     * The referenced containers are invariant under this operation.
     * In every case, a new object is created and returned.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <> struct Product <tags::CPU::MultiCore> : MCProduct <tags::CPU::MultiCore> {};
    template <> struct Product <tags::CPU::MultiCore::SSE> : MCProduct <tags::CPU::MultiCore::SSE> {};

}
#endif
