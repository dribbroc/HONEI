/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007, 2009 Sven Mallach <mallach@honei.org>
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2011 Markus Geveler <apryde@gmx.de>
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

#pragma once
#ifndef LIBLA_GUARD_PRODUCT_HH
#define LIBLA_GUARD_PRODUCT_HH 1

#include <honei/la/banded_matrix.hh>
#include <honei/la/algorithm.hh>
#include <honei/la/banded_matrix_qx.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/dense_matrix_tile.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/dot_product.hh>
#include <honei/la/element_product.hh>
#include <honei/la/matrix_error.hh>
#include <honei/la/scale.hh>
#include <honei/la/scaled_sum.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/la/sparse_matrix_ell.hh>
#include <honei/la/sparse_vector.hh>
#include <honei/la/sum.hh>
#include <honei/util/benchmark_info.hh>
#include <honei/util/tags.hh>
#include <honei/la/algorithm.hh>
#include <honei/mpi/operations.hh>
#include <honei/mpi/dense_vector_mpi-fwd.hh>
#include <honei/mpi/sparse_matrix_ell_mpi-fwd.hh>
#include <honei/mpi/sparse_matrix_csr_mpi-fwd.hh>

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
            for (typename DenseVector<DT1_>::ElementIterator r(result.begin_elements()), r_end(result.end_elements()) ;
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
            for (typename DenseVector<DT1_>::ElementIterator r(result.begin_elements()), r_end(result.end_elements()) ;
                    r != r_end ; ++r)
            {
                const DenseVectorRange<DT1_> dv(a[r.index()]);
                *r = DotProduct<Tag_>::value(dv, b);
            }

            return result;
        }

        template <typename DT1_, typename DT2_>
        static DenseVector<DT1_> value(const SparseMatrix<DT1_> & a, const DenseVectorBase<DT2_> & b)
        {
            CONTEXT("When multiplying SparseMatrix with DenseVectorBase:");

            if (b.size() != a.columns())
            {
                throw VectorSizeDoesNotMatch(b.size(), a.columns());
            }

            DenseVector<DT1_> result(a.rows(), DT1_(0));
            for (typename DenseVector<DT1_>::ElementIterator r(result.begin_elements()), r_end(result.end_elements()) ;
                    r != r_end ; ++r)
            {
                const SparseVector<DT1_> sv(a[r.index()]);
                *r = DotProduct<Tag_>::value(sv, b);
            }

            return result;
        }

        template <typename DT1_, typename DT2_>
        static DenseVector<DT1_> value(DenseVector<DT1_> & result, const SparseMatrixELL<DT1_> & a, const DenseVector<DT2_> & b,
                unsigned long row_start = 0, unsigned long row_end = 0)
        {
            CONTEXT("When multiplying SparseMatrixELL with DenseVector:");

            if (b.size() != a.columns())
            {
                throw VectorSizeDoesNotMatch(b.size(), a.columns());
            }
            if (a.rows() != result.size())
            {
                throw VectorSizeDoesNotMatch(a.rows(), result.size());
            }

            //DenseVector<DT1_> result(a.rows(), DT1_(0));
            //fill<tags::CPU>(result, DT1_(0));

            /*for(unsigned long n(0) ; n < a.num_cols_per_row() ; n++)
            {
                const unsigned long * Aj_n = a.Aj().elements() + n * a.stride();
                const DT1_ * Ax_n = a.Ax().elements() + n * a.stride();
                for(unsigned i(0) ; i < a.rows() ; i++)
                {
                    if(Ax_n[i] != DT1_(0))
                        result[i] += Ax_n[i] * b[Aj_n[i]];
                }
            }*/
            /*const unsigned long size(a.Ax().size());
            const unsigned long stride(a.stride());
            const DT1_ * aax(a.Ax().elements());
            const unsigned long * aaj(a.Aj().elements());
            for (unsigned long i(0) ; i < size ; ++i)
            {
                result.elements()[i % stride] += aax[i] * b.elements()[aaj[i]];
            }*/

            if (row_end == 0)
                row_end = a.rows();

            const unsigned long stride(a.stride());
            const DT1_ * bx(b.elements());
            const DT1_ * aax(a.Ax().elements());
            const unsigned long * aaj(a.Aj().elements());
            const unsigned long * aarl(a.Arl().elements());
            DT1_ sum(0);
            for (unsigned long row(row_start) ; row < row_end ; ++row)
            {
                sum = 0;
                for (unsigned long col(0), j(row * a.threads()) ; col < aarl[row] ; ++col, j+=stride)
                {
                    for (unsigned long thread(0) ; thread < a.threads() ; ++thread)
                        sum+= aax[j + thread] * bx[aaj[j + thread]];
                }
                result.elements()[row] = sum;
            }

            return result;
        }

        template <typename DT_>
        static DenseVector<DT_> & value(DenseVector<DT_> & rv, const SparseMatrixCSR<DT_> & a, const DenseVector<DT_> & bv,
                unsigned long row_start = 0, unsigned long row_end = 0)
        {
            CONTEXT("When multiplying SparseMatrixCSR with DenseVector:");

            if (bv.size() != a.columns())
            {
                throw VectorSizeDoesNotMatch(bv.size(), a.columns());
            }
            if (row_end == 0)
                row_end = a.rows();

            const unsigned long * const Ar(a.Ar().elements());
            const unsigned long * const Aj(a.Aj().elements());
            const DT_ * const Ax(a.Ax().elements());
            const DT_ * const b(bv.elements());
            DT_ * r(rv.elements());
            const unsigned long blocksize(a.blocksize());

            for (unsigned long row(row_start) ; row < row_end ; ++row)
            {
                DT_ sum(0);
                const unsigned long end(Ar[row+1]);
                for (unsigned long i(Ar[row]) ; i < end ; ++i)
                {
                    for (unsigned long blocki(0) ; blocki < blocksize ; ++blocki)
                    {
                        sum += Ax[(i*blocksize)+blocki] * b[Aj[i] + blocki];
                    }
                }
                r[row] = sum;
            }

            return rv;
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
            unsigned long middle_index(a.rows() -1);
            for (typename BandedMatrix<DT1_>::ConstBandIterator vi(a.begin_non_zero_bands()), vi_end(a.end_non_zero_bands()) ;
                    vi != vi_end ; ++vi)
            {
                // If we are above or on the diagonal band, we start at Element 0 and go on until Element band_size-band_index.
                if (middle_index < vi.index())
                {
                    typename DenseVector<DT2_>::ConstElementIterator j(b.begin_elements()), j_end(b.end_elements());
                    typename DenseVector<DT1_>::ElementIterator r(result.begin_elements()), r_end(result.end_elements());

                    for (unsigned int i(0) ; i < (vi.index()-middle_index) && j != j_end ; ++i)
                    {
                        ++j; // Get the right position in b.
                    }

                    //Calculation of the element-index to stop in iteration!
                    unsigned long end(vi->size() - (vi.index() - middle_index));
                    for(typename DenseVector<DT1_>::ConstElementIterator c(vi->begin_elements()),
                            c_end(vi->element_at(end)) ; c < c_end ; ++c)
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
                    typename DenseVector<DT2_>::ConstElementIterator j(b.begin_elements()), j_end(b.end_elements());
                    typename DenseVector<DT1_>::ElementIterator r(result.begin_elements()), r_end(result.end_elements());
                    unsigned long start(middle_index - vi.index()); //Calculation of the element-index to start in iteration!
                    for (unsigned int a(0); a < middle_index - vi.index() && r != r_end; ++a)
                    {
                        ++r; // Get the right position in result b.
                    }
                    for(typename DenseVector<DT1_>::ConstElementIterator c(vi->element_at(start)),
                            c_end(vi->end_elements()) ; c < c_end ; ++c)
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
        static DenseVector<DT1_> value(const BandedMatrixQx<Q1Type, DT1_> & a, const DenseVectorBase<DT2_> & b)
        {
            CONTEXT("When multiplying BandedMatrixQ1 with DenseVectorBase:");
            if (b.size() != a.columns())
            {
                throw VectorSizeDoesNotMatch(b.size(), a.columns());
            }

            DenseVector<DT1_> result(a.rows());
            long root(a.root());

            for (long index(0) ; index < (long)b.size() ; ++index)
            {
                result[index] = a.band(DD)[index] * b[index];
                if ((index - root - 1) >= 0)
                    result[index] += a.band(LL)[index] * b[index - root - 1];
                if ((index - root) >= 0)
                    result[index] += a.band(LD)[index] * b[index - root];
                if ((index - root + 1) >= 0)
                    result[index] += a.band(LU)[index] * b[index - root + 1];
                if ((index - 1) >= 0)
                    result[index] += a.band(DL)[index] * b[index - 1];
                if ((index + 1) < (long)b.size())
                    result[index] += a.band(DU)[index] * b[index + 1];
                if ((index + root - 1) < (long)b.size())
                    result[index] += a.band(UL)[index] * b[index + root - 1];
                if ((index + root) < (long)b.size())
                    result[index] += a.band(UD)[index] * b[index + root];
                if ((index + root + 1) < (long)b.size())
                    result[index] += a.band(UU)[index] * b[index + root + 1];
            }

            return result;
        }

        template <typename DT1_, typename DT2_>
        static DenseVector<DT1_> & value(DenseVector<DT1_> & result, const BandedMatrixQx<Q1Type, DT1_> & a, const DenseVectorBase<DT2_> & b)
        {
            CONTEXT("When multiplying BandedMatrixQ1 with DenseVectorBase:");
            if (b.size() != a.columns())
            {
                throw VectorSizeDoesNotMatch(b.size(), a.columns());
            }

            long root(a.root());

            for (long index(0) ; index < (long)b.size() ; ++index)
            {
                result[index] = a.band(DD)[index] * b[index];
                if ((index - root - 1) >= 0)
                    result[index] += a.band(LL)[index] * b[index - root - 1];
                if ((index - root) >= 0)
                    result[index] += a.band(LD)[index] * b[index - root];
                if ((index - root + 1) >= 0)
                    result[index] += a.band(LU)[index] * b[index - root + 1];
                if ((index - 1) >= 0)
                    result[index] += a.band(DL)[index] * b[index - 1];
                if ((index + 1) < (long)b.size())
                    result[index] += a.band(DU)[index] * b[index + 1];
                if ((index + root - 1) < (long)b.size())
                    result[index] += a.band(UL)[index] * b[index + root - 1];
                if ((index + root) < (long)b.size())
                    result[index] += a.band(UD)[index] * b[index + root];
                if ((index + root + 1) < (long)b.size())
                    result[index] += a.band(UU)[index] * b[index + root + 1];
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

            for (typename BandedMatrix<DT1_>::ConstBandIterator vi(a.begin_bands()),
                    vi_end(a.band_at(middle_index)) ; vi != vi_end; ++vi)
            {
                if (! vi.exists())
                    continue;

                // If we are below the diagonal band we correct the position in band and result b.
                unsigned long move_index(middle_index - vi.index()); // is also = the element index to start in iteration for b
                for (typename SparseVector<DT2_>::NonZeroConstElementIterator j(b.begin_non_zero_elements()),
                        j_end(b.end_non_zero_elements()) ; j != j_end ; ++j)
                {
                    typename DenseVector<DT1_>::ConstElementIterator c(vi->element_at(move_index)), c_end(vi->end_elements());
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

            for (typename BandedMatrix<DT1_>::ConstBandIterator vi(a.band_at(middle_index)),
                    vi_end(a.end_bands()) ;  vi != vi_end ; ++vi)
            {
                if (! vi.exists())
                    continue;

                // If we are above or on the diagonal band, we first search the first non-zero element, then
                // we make positions in band and result b meet it.
                unsigned long move_index(vi.index() - middle_index);
                unsigned long end(vi->size() - move_index); //Calculation of the element-index to stop in iteration!

                for (typename SparseVector<DT2_>::NonZeroConstElementIterator j(b.begin_non_zero_elements()),
                        j_end(b.end_non_zero_elements()) ; j != j_end ; ++j)
                {
                    if (j.index() < move_index)
                    {
                        continue;
                    }

                    typename DenseVector<DT1_>::ConstElementIterator c(vi->begin_elements()), c_end(vi->element_at(end));
                    typename SparseVector<DT1_>::ElementIterator r(result.begin_elements());

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

            for (typename DenseVector<DT1_>::ConstElementIterator i(a.element_at(1)), i_end(a.end_elements()) ;
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

            for (typename DenseVector<DT1_>::ConstElementIterator i(a.element_at(1)), i_end(a.end_elements()) ;
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
            typename DenseMatrix<DT1_>::ElementIterator i(result.begin_elements());

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

            typename DenseMatrixTile<DT1_>::ElementIterator i(r.begin_elements());

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
            for( typename SparseMatrix<DT2_>::NonZeroConstElementIterator i(b.begin_non_zero_elements()), i_end(b.end_non_zero_elements()) ;
                    i < i_end ; ++i )
            {
                typename DenseMatrix<DT1_>::Column column(result.column(i.column()));
                ScaledSum<>::value(column, a.column(i.row()), *i);
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
            for (typename SparseMatrix<DT1_>::NonZeroConstElementIterator i(a.begin_non_zero_elements()) ; i != a.end_non_zero_elements() ; ++i)
            {
                for (typename SparseVector<DT1_>::NonZeroConstElementIterator j(b[i.column()].begin_non_zero_elements()) ; j != b[i.column()].end_non_zero_elements() ; ++j)
                {
                    result(i.row(), j.index()) += *i * *j;
                }
            }
            return result;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrixELL<DT1_> value(const SparseMatrixELL<DT1_> & a, const SparseMatrixELL<DT2_> & b)
        {
            CONTEXT("When multiplying SparseMatrixELL with SparseMatrixELL:");

            if (a.columns() != b.rows())
                throw MatrixRowsDoNotMatch(b.rows(), a.columns());

            SparseMatrix<DT1_> result(a.rows(), b.columns(), b.num_cols_per_row());
            const DT1_ * aax = a.Ax().elements();
            const DT2_ * bax = b.Ax().elements();
            const unsigned long astride(a.stride());
            const unsigned long bstride(b.stride());
            const unsigned long * arl(a.Arl().elements());
            const unsigned long * brl(b.Arl().elements());
            const unsigned long * aaj(a.Aj().elements());
            const unsigned long * baj(b.Aj().elements());

            for (unsigned long row(0) ; row < a.rows() ; ++row)
            {
                for(unsigned long ac(0), i(row) ; ac < arl[row] ; i+=astride, ++ac)
                {
                    for (unsigned long bc(0), j(aaj[i]) ; bc < brl[aaj[i]] ; j+=bstride, ++bc)
                    {
                        result(row, baj[j]) += aax[i] * bax[j];
                    }
                }
            }

            SparseMatrixELL<DT1_> resultell(result);
            return resultell;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> value(const SparseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
        {
            CONTEXT("When multiplying SparseMatrix with DenseMatrix:");

            if (a.columns() != b.rows())
                throw MatrixRowsDoNotMatch(b.rows(), a.columns());

            DenseMatrix<DT1_> result(a.rows(), b.columns(), DT1_(0));
            for (typename SparseMatrix<DT1_>::NonZeroConstElementIterator i(a.begin_non_zero_elements()), i_end(a.end_non_zero_elements()) ;
                    i < i_end ; ++i)
            {
                typename DenseMatrix<DT1_>::Row row(result[i.row()]);
                ScaledSum<>::value(row, b[i.column()], *i);
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
            for (typename BandedMatrix<DT1_>::ConstBandIterator vi(a.begin_bands()),
                    vi_end(a.band_at(diag_index)) ; vi != vi_end ; ++vi)
            {
                if (! vi.exists())
                    continue;

                // We start at the zero_based band_index of b that is equal to abs(band_index(a) - diag_index)
                    unsigned long iteration_start(abs(vi.index() - diag_index));
                    for(typename BandedMatrix<DT1_>::ConstBandIterator vj(b.band_at(iteration_start)),
                            vj_end(b.end_bands()) ; vj != vj_end ; ++vj)
                    {
                        if (vj.index() == diag_index) // We are on diagonal of b
                        {
                            signed long result_band_index(vi.index() - diag_index); //index based on diag = 0
                            typename DenseVector<DT1_>::ConstElementIterator i(vi->element_at(abs(result_band_index)));
                            typename DenseVector<DT2_>::ConstElementIterator j(vj->begin_elements());
                            for (typename DenseVector<DT1_>::ElementIterator
                                    r(result.band(result_band_index).element_at(abs(result_band_index))),
                                    r_end(result.band(result_band_index).end_elements()) ; r < r_end ;
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
                            typename DenseVector<DT1_>::ConstElementIterator i(vi->element_at(shift));
                            typename DenseVector<DT2_>::ConstElementIterator j(vj->begin_elements());
                            long res_end(result.size());
                            if (result_band_index > 0)
                                res_end -= result_band_index;

                            for(typename DenseVector<DT1_>::ElementIterator r(result.band(result_band_index).element_at(shift)),
                                    r_end(result.band(result_band_index).element_at(res_end)) ; r < r_end ; ++j, ++i, ++r)
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
                            typename DenseVector<DT1_>::ConstElementIterator i(vi->element_at(abs(result_band_index)));
                            typename DenseVector<DT2_>::ConstElementIterator j(vj->element_at(shift));
                            for(typename DenseVector<DT1_>::ElementIterator r(result.band(result_band_index).element_at(abs(result_band_index))), 
                                    r_end(result.band(result_band_index).end_elements()) ; r < r_end ; ++j, ++i, ++r)
                            {
                                *r += *i * *j;
                            }
                        }
                    }
                }

            // Diagonal of a
            typename BandedMatrix<DT1_>::ConstBandIterator vi(a.band_at(diag_index));

            if (vi.exists())
            {
                // Lower part of b
                for (typename BandedMatrix<DT2_>::ConstBandIterator vj(b.begin_bands()),
                        vj_end(b.band_at(diag_index)) ; vj != vj_end ; ++vj)
                {
                    if (! vj.exists())
                        continue;

                    signed long result_band_index(vj.index() - diag_index);  //index based on diag = 0
                    typename DenseVector<DT1_>::ElementIterator r(result.band(result_band_index).element_at(abs(result_band_index)));
                    typename DenseVector<DT1_>::ConstElementIterator i(vi->element_at(abs(result_band_index)));
                    for(typename DenseVector<DT2_>::ConstElementIterator j(vj->element_at(abs(result_band_index))),
                            j_end(vj->end_elements()) ; j < j_end ; ++j, ++i, ++r)
                    {
                        *r += *i * *j;
                    }
                }

                // Diagonal of b
                for(typename BandedMatrix<DT2_>::ConstBandIterator vj(b.band_at(diag_index)),
                        vj_end(b.band_at(diag_index + 1)) ; vj != vj_end ; ++vj)
                {
                        if (! vj.exists())
                            continue;

                        DenseVector<DT1_> temp(vi->copy());
                        Sum<>::value(result.band(0), ElementProduct<>::value(temp, *vj));
                }

                // Upper part of b
                for(typename BandedMatrix<DT2_>::ConstBandIterator vj(b.band_at(diag_index+1)),
                        vj_end(b.end_bands()) ; vj != vj_end ; ++vj)
                {
                    signed long result_band_index(vj.index() - diag_index); //index based on diag = 0
                    typename DenseVector<DT1_>::ElementIterator r(result.band(result_band_index).begin_elements());
                    typename DenseVector<DT1_>::ConstElementIterator i(vi->begin_elements());
                    for(typename DenseVector<DT2_>::ConstElementIterator j(vj->begin_elements()),
                            j_end(vj->element_at(vj->size() - result_band_index)) ; j < j_end ; ++j, ++i, ++r)
                    {
                        *r += *i * *j;
                    }
                }
            }

            // Upper part of a
            for(typename BandedMatrix<DT1_>::ConstBandIterator vi(a.band_at(diag_index+1)),
                    vi_end(a.end_bands()) ; vi != vi_end ; ++vi)
            {
                if (! vi.exists())
                    continue;

                    // We only go on until band_index of b is equal to [numbers of bands] - actual band_index of a (zero_based).
                    unsigned long iteration_end((2*b.size()-1) - (vi.index() - diag_index));
                    for(typename BandedMatrix<DT1_>::ConstBandIterator vj(b.begin_bands()), vj_end(b.band_at(iteration_end)) ; vj != vj_end ; ++vj)
                    {
                        if (vj.index() == diag_index) // We are on diagonal of b
                        {
                            signed long result_band_index(vi.index() - diag_index); //index based on diag = 0
                            typename DenseVector<DT1_>::ConstElementIterator i(vi->begin_elements());
                            typename DenseVector<DT2_>::ConstElementIterator j(vj->element_at(result_band_index));
                            for(typename DenseVector<DT1_>::ElementIterator r(result.band(result_band_index).begin_elements()),
                                    r_end(result.band(result_band_index).element_at(result.size() - result_band_index)) ; r < r_end ; ++j, ++i, ++r)
                            {
                                *r += *i * *j;
                            }
                        }
                        else if (vj.index() > diag_index) // We are above diagonal of b
                        {
                            signed long diag_based_index_a(vi.index() - diag_index);
                            signed long diag_based_index_b(vj.index() - diag_index);
                            signed long result_band_index(diag_based_index_a + diag_based_index_b);
                            typename DenseVector<DT1_>::ConstElementIterator i(vi->begin_elements());
                            unsigned long shift(vi.index() - diag_index);
                            typename DenseVector<DT2_>::ConstElementIterator j(vj->element_at(shift));
                            for(typename DenseVector<DT1_>::ElementIterator r(result.band(result_band_index).begin_elements()),
                                    r_end(result.band(result_band_index).element_at(result.size() - result_band_index)) ; r < r_end ; ++j, ++i, ++r)
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

                            typename DenseVector<DT1_>::ConstElementIterator i(vi->element_at(res_start));
                            long vj_start(vi.index() - diag_index);
                            if (result_band_index < 0)
                                vj_start += abs(result_band_index);

                            typename DenseVector<DT2_>::ConstElementIterator j(vj->element_at(vj_start));
                            long res_end(2*result.size() - (vi.index() + 1));
                            for(typename DenseVector<DT1_>::ElementIterator r(result.band(result_band_index).element_at(res_start)),
                                    r_end(result.band(result_band_index).element_at(res_end)) ; r < r_end ; ++j, ++i, ++r)
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
            for (typename BandedMatrix<DT1_>::ConstBandIterator vi(a.begin_bands()),
                    vi_end(a.band_at(middle_index)) ; vi != vi_end ; ++vi)
            {
                if (! vi.exists())
                    continue;

                for (unsigned int s(0) ; s < b.columns() ; ++s)
                {
                    // Temporary container for efficient calculation of elementwise vector product.
                    DenseVector<DT2_> temp(b.rows(), DT2_(0));
                    unsigned long real_index(middle_index - vi.index());
                    typename DenseVector<DT2_>::ConstElementIterator c(b.column(s).begin_elements()),
                                d(vi->element_at(real_index));
                    for (typename DenseVector<DT2_>::ElementIterator x(temp.element_at(real_index)),
                                x_end(temp.end_elements()) ; x != x_end ; ++x)
                    {
                        *x = *c * *d;
                        ++c; ++d;
                    }

                    typename DenseMatrix<DT2_>::Column col(result.column(s));
                    Sum<>::value(col, temp);
                }
            }

            // Calculation for diagonal
            typename BandedMatrix<DT1_>::ConstBandIterator vi(a.band_at(middle_index));

            if (vi.exists())
            {
                for (unsigned int s(0) ; s < b.columns() ; ++s)
                {
                    DenseVector<DT2_> temp(b.column(s).copy());
                    typename DenseMatrix<DT2_>::Column col(result.column(s));
                    Sum<>::value(col, ElementProduct<>::value(temp, *vi));
                }
            }

            // Calculation for upper part
            for (typename BandedMatrix<DT1_>::ConstBandIterator vi(a.band_at(middle_index+1)),
                    vi_end(a.end_bands()) ; vi != vi_end ; ++vi)
            {
                if (! vi.exists())
                    continue;

                for (unsigned int s(0) ; s < b.columns() ; ++s)
                {
                    // Temporary container for efficient calculation of elementwise vector product.
                    DenseVector<DT2_> temp(b.rows(), DT2_(0));

                    unsigned long real_index(vi.index() - middle_index);
                    typename DenseVector<DT2_>::ConstElementIterator c(b.column(s).element_at(real_index)),
                            d(vi->begin_elements());
                    unsigned long end(temp.size() - real_index);

                    for (typename DenseVector<DT2_>::ElementIterator x(temp.begin_elements()),
                            x_end(temp.element_at(end)) ; x != x_end ; ++x)
                    {
                        *x = *c * *d;
                        ++c; ++d;
                    }

                    typename DenseMatrix<DT2_>::Column col(result.column(s));
                    Sum<>::value(col, temp);
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
            for (typename BandedMatrix<DT1_>::ConstBandIterator vi(a.begin_bands()),
                    vi_end(a.band_at(middle_index)) ; vi != vi_end ; ++vi)
            {
                if (! vi.exists())
                    continue;

                unsigned long real_index(middle_index - vi.index());

                typename DenseVector<DT1_>::ConstElementIterator d(vi->element_at(real_index)),
                        d_end(vi->end_elements());
                for (unsigned int z(0) ; z < (b.rows()-real_index) ; ++z, ++d)
                {
                    const SparseVector<DT2_> row(b[z]);
                    typename DenseVector<DT2_>::ElementIterator x(result[d.index()].begin_elements());
                    for(typename SparseVector<DT2_>::ConstElementIterator c(row.begin_elements()),
                            c_end(row.end_elements()) ; c != c_end ; ++x, ++c)
                    {
                        *x += *d * *c;
                    }
                }
            }

            // Calculation for diagonal
            typename BandedMatrix<DT1_>::ConstBandIterator vi(a.band_at(middle_index));

            if (vi.exists())
            {
                typename DenseVector<DT1_>::ConstElementIterator d(vi->begin_elements());
                for (unsigned int z(0) ; z < b.rows() ; ++z, ++d)
                {
                    const SparseVector<DT2_> row(b[z]);
                    typename DenseVector<DT2_>::ElementIterator x(result[z].begin_elements());
                    for(typename SparseVector<DT2_>::ConstElementIterator c(row.begin_elements()),
                            c_end(row.end_elements()) ; c != c_end ; ++x, ++c)
                    {
                        *x += *d * *c;
                    }
                }
            }

            // Calculation for upper part
            for (typename BandedMatrix<DT1_>::ConstBandIterator vi(a.band_at(middle_index+1)),
                    vi_end(a.end_bands()) ; vi != vi_end ; ++vi)
            {
                if (! vi.exists())
                    continue;

                unsigned long real_index(vi.index() - middle_index);

                typename DenseVector<DT1_>::ConstElementIterator d(vi->begin_elements()),
                        d_end(vi->element_at(vi->size()-real_index));
                for (unsigned int z(real_index) ; z < b.rows() ; ++z, ++d)
                {
                    const SparseVector<DT2_> row(b[z]);
                    typename DenseVector<DT2_>::ElementIterator x(result[d.index()].begin_elements());
                    for(typename SparseVector<DT2_>::ConstElementIterator c(row.begin_elements()), c_end(row.end_elements()) ;
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
            for (typename BandedMatrix<DT2_>::ConstBandIterator vi(b.begin_bands()),
                    vi_end(b.band_at(middle_index)) ; vi != vi_end ; ++vi)
            {
                if (! vi.exists())
                    continue;

                unsigned long real_index(vi.index() - middle_index);

                for (unsigned int z(0) ; z < a.rows() ; ++z)
                {
                    /// \todo RowIterator
                    const DenseVectorRange<DT1_> row(a[z]);
                    typename DenseVector<DT1_>::ConstElementIterator c(row.begin_elements());
                    typename DenseVector<DT1_>::ElementIterator x(result[z].element_at(real_index));
                    for(typename DenseVector<DT2_>::ConstElementIterator d(vi->begin_elements()),
                        d_end(vi->element_at(vi->size()-real_index)) ; d < d_end ; ++x, ++c, ++d)
                    {
                        *x += *d * *c;
                    }
                }
            }

            // Calculation for diagonal
            typename BandedMatrix<DT2_>::ConstBandIterator vi(b.band_at(middle_index));

            if (vi.exists())
            {
                for (unsigned int z(0) ; z < a.rows() ; ++z)
                {
                    const DenseVectorRange<DT1_> row(a[z]);
                    typename DenseVector<DT1_>::ConstElementIterator c(row.begin_elements());
                    typename DenseVector<DT1_>::ElementIterator x(result[z].begin_elements());
                    for(typename DenseVector<DT2_>::ConstElementIterator d(vi->begin_elements()),
                        d_end(vi->end_elements()) ; d < d_end ; ++x, ++c, ++d)
                    {
                        *x += *d * *c;
                    }
                }
            }

            // Calculation for upper part
            for (typename BandedMatrix<DT2_>::ConstBandIterator vi(b.band_at(middle_index+1)),
                    vi_end(b.end_bands()) ; vi != vi_end ; ++vi)
            {
                if (! vi.exists())
                    continue;

                unsigned long real_index(middle_index - vi.index());

                for (unsigned int z(0) ; z < a.rows() ; ++z)
                {
                    const DenseVectorRange<DT1_> row(a[z]);
                    typename DenseVector<DT1_>::ConstElementIterator c(row.element_at(real_index));
                    typename DenseVector<DT1_>::ElementIterator x(result[z].begin_elements());
                    for(typename DenseVector<DT2_>::ConstElementIterator d(vi->element_at(real_index)),
                        d_end(vi->end_elements()) ; d < d_end ; ++x, ++c, ++d)
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
            for (typename BandedMatrix<DT2_>::ConstBandIterator vi(b.begin_bands()),
                    vi_end(b.band_at(middle_index)) ; vi != vi_end ; ++vi)
            {
                if (! vi.exists())
                    continue;

                unsigned long real_index(middle_index - vi.index());

                for (unsigned int z(0) ; z < a.rows() ; ++z)
                {
                    const SparseVector<DT1_> row(a[z]);
                    typename DenseVectorRange<DT1_>::ElementIterator x(result[z].begin_elements());
                    for(typename DenseVector<DT2_>::ConstElementIterator d(vi->element_at(real_index)),
                        d_end(vi->end_elements()) ; d < d_end ; ++x, ++d)
                    {
                        *x += *d * row[real_index + x.index()];
                    }
                }
            }
            // Calculation for diagonal part
            typename BandedMatrix<DT2_>::ConstBandIterator vi(b.band_at(middle_index));

            if (vi.exists())
            {
                for (unsigned int z(0) ; z < a.rows() ; ++z)
                {
                    const SparseVector<DT1_> row(a[z]);
                    typename SparseVector<DT1_>::ConstElementIterator c(row.begin_elements());
                    typename DenseVectorRange<DT1_>::ElementIterator x(result[z].begin_elements());
                    for(typename DenseVector<DT2_>::ConstElementIterator d(vi->begin_elements()),
                        d_end(vi->end_elements()) ; d < d_end ; ++x, ++c, ++d)
                    {
                        *x += *d * *c;
                    }
                }
            }

            // Calculation for upper part
            for (typename BandedMatrix<DT2_>::ConstBandIterator vi(b.band_at(middle_index+1)),
                    vi_end(b.end_bands()) ; vi != vi_end ; ++vi)
            {
                if (! vi.exists())
                    continue;

                const unsigned long real_index(vi.index() - middle_index);

                for (unsigned int z(0) ; z < a.rows() ; ++z)
                {
                    const SparseVector<DT1_> row(a[z]);
                    typename SparseVector<DT1_>::ConstElementIterator c(row.begin_elements());
                    typename DenseVectorRange<DT1_>::ElementIterator x(result[z].element_at(real_index));
                    for(typename DenseVector<DT2_>::ConstElementIterator d(vi->begin_elements()),
                        d_end(vi->element_at(vi->size()-real_index)) ; d < d_end ; ++x, ++c, ++d)
                    {
                        *x += *d * *c;
                    }
                }
            }

            return result;
        }

        template<typename DT1_, typename DT2_>
        static DenseVectorContinuousBase<DT1_> & value(DenseVectorContinuousBase<DT1_> & y, const DenseVectorContinuousBase<DT1_> & a, const DenseVectorContinuousBase<DT2_> & b)
        {
            ElementProduct<Tag_>::value(y, a, b);
            return y;
        }


        template <typename DT_>
        static inline DenseVectorMPI<DT_> & value(DenseVectorMPI<DT_> & r, const SparseMatrixELLMPI<DT_> & a, const DenseVectorMPI<DT_> & b)
        {
            MPIOps<tags::CPU>::product(r, a, b);
            return r;
        }

        template <typename DT_>
        static inline DenseVectorMPI<DT_> & value(DenseVectorMPI<DT_> & r, const SparseMatrixCSRMPI<DT_> & a, const DenseVectorMPI<DT_> & b)
        {
            MPIOps<tags::CPU>::product(r, a, b);
            return r;
        }

        template<typename DT_>
        static inline DenseVectorMPI<DT_> & value(DenseVectorMPI<DT_> & r, const DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y)
        {
            MPIOps<tags::CPU>::element_product(r, x, y);
            return r;
        }

        /// \}

        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(BandedMatrix<DT1_> & a, DenseVectorBase<DT2_> & b)
        {
            BenchmarkInfo result;
            unsigned long temp(0);
            for (typename BandedMatrix<DT1_>::ConstBandIterator vi(a.begin_non_zero_bands()), vi_end(a.end_non_zero_bands()) ;
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
            for (typename BandedMatrix<DT1_>::ConstBandIterator vi(a.begin_bands()),
                    vi_end(a.band_at(middle_index)) ; vi != vi_end ; ++vi)
            {
                if (! vi.exists())
                    continue;

                for (unsigned int s(0) ; s < b.columns() ; ++s)
                {
                    // Temporary container for efficient calculation of elementwise vector product.
                    DenseVector<DT2_> temp(b.rows(), DT2_(0));
                    unsigned long real_index(middle_index - vi.index());
                    typename DenseVectorSlice<DT2_>::ConstElementIterator c(b.column(s).begin_elements()),
                                d(vi->element_at(real_index));
                    for (typename DenseVectorRange<DT2_>::ElementIterator x(temp.element_at(real_index)),
                                x_end(temp.end_elements()) ; x != x_end ; ++x)
                    {
                        //*x = *c * *d;
                        result.flops += 1;
                        result.load += sizeof(DT2_) + sizeof(DT1_);
                        result.store += sizeof(DT2_);
                        ++c; ++d;
                    }

                    //Sum<>::value(presult.column(s), temp);
                    result = result + Sum<>::get_benchmark_info(temp, presult.column(s));
                }
            }

            // Calculation for diagonal
            typename BandedMatrix<DT1_>::ConstBandIterator vi(a.band_at(middle_index));

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
            for (typename BandedMatrix<DT1_>::ConstBandIterator vi(a.band_at(middle_index+1)),
                    vi_end(a.end_bands()) ; vi != vi_end ; ++vi)
            {
                if (! vi.exists())
                    continue;

                for (unsigned int s(0) ; s < b.columns() ; ++s)
                {
                    // Temporary container for efficient calculation of elementwise vector product.
                    DenseVector<DT2_> temp(b.rows(), DT2_(0));

                    unsigned long real_index(vi.index() - middle_index);
                    typename DenseVectorSlice<DT2_>::ConstElementIterator c(b.column(s).element_at(real_index)),
                            d(vi->begin_elements());
                    unsigned long end(temp.size() - real_index);

                    for (typename DenseVector<DT2_>::ElementIterator x(temp.begin_elements()),
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

        template <typename DT_>
        static inline BenchmarkInfo get_benchmark_info(const DenseVectorContinuousBase<DT_> & r, const SparseMatrixELL<DT_> & a, const DenseVectorContinuousBase<DT_> & b)
        {
            BenchmarkInfo result;
            result.flops = a.used_elements() * 2;
            result.load = (a.used_elements() + b.size() )* sizeof(DT_);
            result.store = r.size() * sizeof(DT_);
            return result;
        }
    };

    template <> struct Product<tags::CPU::Generic>
    {
        template <typename DT_>
        static DenseVector<DT_> & value(DenseVector<DT_> & r, const SparseMatrix<DT_> & a, const DenseVector<DT_> & b,
                unsigned long row_start = 0, unsigned long row_end = 0)
        {
            if (b.size() != a.columns())
            {
                throw VectorSizeDoesNotMatch(b.size(), a.columns());
            }
            if (row_end == 0)
                row_end = a.rows();

            BENCHADD(Product<tags::CPU>::get_benchmark_info(r, a, b));

            const DT_ * const b_e(b.elements());
            for (unsigned long row(row_start) ; row < row_end ; ++row)
            {
                const DT_ * const row_e(a[row].elements());
                const unsigned long * const row_i(a[row].indices());
                const unsigned long ue(a[row].used_elements());
                DT_ sum(0);
                for (unsigned long i(0) ; i < ue ; ++i)
                {
                    const unsigned long idx(row_i[i]);
                    sum += row_e[i] * b_e[idx];
                }
                r[row] = sum;
            }

            return r;
        }

        template <typename DT_>
        static DenseVector<DT_> & value(DenseVector<DT_> & r, const SparseMatrixELL<DT_> & a, const DenseVector<DT_> & bv,
                unsigned long row_start = 0, unsigned long row_end = 0)
        {
            if (bv.size() != a.columns())
            {
                throw VectorSizeDoesNotMatch(bv.size(), a.columns());
            }
            if (row_end == 0)
                row_end = a.rows();

            BENCHADD(Product<tags::CPU>::get_benchmark_info(r, a, bv));

            DT_ * result(r.elements());
            const unsigned long * Aj(a.Aj().elements());
            const DT_ * Ax(a.Ax().elements());
            const unsigned long * Arl(a.Arl().elements());
            const DT_ * b(bv.elements());
            const unsigned long stride(a.stride());
            const unsigned long threads(a.threads());

            for (unsigned long row(row_start) ; row < row_end ; ++row)
            {
                const unsigned long * tAj(Aj);
                const DT_ * tAx(Ax);
                DT_ sum(0);
                tAj += row * threads;
                tAx += row * threads;

                const unsigned long max(Arl[row]);
                for(unsigned long n = 0; n < max ; n++)
                {
                    for (unsigned long thread(0) ; thread < threads ; ++thread)
                    {
                        const DT_ A_ij = *(tAx + thread);

                        //if (A_ij != 0)
                        {
                            const unsigned long col = *(tAj + thread);
                            sum += A_ij * b[col];
                        }
                    }

                    tAj += stride;
                    tAx += stride;
                }
                result[row] = sum;
            }

            return r;
        }

        template <typename DT_>
        static DenseVector<DT_> & value(DenseVector<DT_> & rv, const SparseMatrixCSR<DT_> & a, const DenseVector<DT_> & bv,
                unsigned long row_start = 0, unsigned long row_end = 0)
        {
            CONTEXT("When multiplying SparseMatrixCSR with DenseVector:");
            if (bv.size() != a.columns())
            {
                throw VectorSizeDoesNotMatch(bv.size(), a.columns());
            }
            if (row_end == 0)
                row_end = a.rows();

            BENCHADD(Product<tags::CPU>::get_benchmark_info(rv, a, bv));

            const unsigned long * const Ar(a.Ar().elements());
            const unsigned long * const Aj(a.Aj().elements());
            const DT_ * const Ax(a.Ax().elements());
            const DT_ * const b(bv.elements());
            DT_ * r(rv.elements());
            const unsigned long blocksize(a.blocksize());

            for (unsigned long row(row_start) ; row < row_end ; ++row)
            {
                DT_ sum(0);
                const unsigned long end(Ar[row+1]);
                for (unsigned long i(Ar[row]) ; i < end ; ++i)
                {
                    for (unsigned long blocki(0) ; blocki < blocksize ; ++blocki)
                    {
                        sum += Ax[(i*blocksize)+blocki] * b[Aj[i] + blocki];
                    }
                }
                r[row] = sum;
            }

            return rv;
        }

        template<typename DT_>
        static inline DenseVectorContinuousBase<DT_> & value(DenseVectorContinuousBase<DT_> & r, const DenseVectorContinuousBase<DT_> & x, const DenseVectorContinuousBase<DT_> & y)
        {
            ElementProduct<tags::CPU::Generic>::value(r, x, y);
            return r;
        }

        template <typename DT_>
        static inline DenseVectorMPI<DT_> & value(DenseVectorMPI<DT_> & r, const SparseMatrixELLMPI<DT_> & a, const DenseVectorMPI<DT_> & b)
        {
            MPIOps<tags::CPU::Generic>::product(r, a, b);
            return r;
        }

        template <typename DT_>
        static inline DenseVectorMPI<DT_> & value(DenseVectorMPI<DT_> & r, const SparseMatrixCSRMPI<DT_> & a, const DenseVectorMPI<DT_> & b)
        {
            MPIOps<tags::CPU::Generic>::product(r, a, b);
            return r;
        }

        template<typename DT_>
        static inline DenseVectorMPI<DT_> & value(DenseVectorMPI<DT_> & r, const DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y)
        {
            MPIOps<tags::CPU::Generic>::element_product(r, x, y);
            return r;
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
    template <> struct Product<tags::OpenCL::CPU>
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

        template <typename DT_>
        static DenseVectorContinuousBase<DT_> & value(DenseVectorContinuousBase<DT_> & result, const SparseMatrixELL<DT_> & a, const DenseVectorContinuousBase<DT_> & b);

        template <typename DT_>
        static DenseVector<DT_> value(const BandedMatrixQx<Q1Type, DT_> & a, const DenseVectorContinuousBase<DT_> & b);

        template <typename DT_>
        static DenseVectorContinuousBase<DT_> & value(DenseVectorContinuousBase<DT_> & result, const BandedMatrixQx<Q1Type, DT_> & a, const DenseVectorContinuousBase<DT_> & b);

        template<typename DT1_, typename DT2_>
        static DenseVectorContinuousBase<DT1_> & value(DenseVectorContinuousBase<DT1_> & y, const DenseVectorContinuousBase<DT1_> & a, const DenseVectorContinuousBase<DT2_> & b)
        {
            ElementProduct<tags::OpenCL::CPU>::value(y, a, b);
            return y;
        }

        /// \}
    };

    template <> struct Product<tags::OpenCL::GPU>
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

        template <typename DT_>
        static DenseVectorContinuousBase<DT_> & value(DenseVectorContinuousBase<DT_> & result, const SparseMatrixELL<DT_> & a, const DenseVectorContinuousBase<DT_> & b);

        template <typename DT_>
        static DenseVector<DT_> value(const BandedMatrixQx<Q1Type, DT_> & a, const DenseVectorContinuousBase<DT_> & b);

        template <typename DT_>
        static DenseVectorContinuousBase<DT_> & value(DenseVectorContinuousBase<DT_> & result, const BandedMatrixQx<Q1Type, DT_> & a, const DenseVectorContinuousBase<DT_> & b);

        template<typename DT1_, typename DT2_>
        static DenseVectorContinuousBase<DT1_> & value(DenseVectorContinuousBase<DT1_> & y, const DenseVectorContinuousBase<DT1_> & a, const DenseVectorContinuousBase<DT2_> & b)
        {
            ElementProduct<tags::OpenCL::GPU>::value(y, a, b);
            return y;
        }

        /// \}
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
    template <> struct Product<tags::GPU::CUDA>
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

        static DenseVector<float> value(const BandedMatrixQx<Q1Type, float> & a, const DenseVectorContinuousBase<float> & b);

        static DenseVector<double> value(const BandedMatrixQx<Q1Type, double> & a, const DenseVectorContinuousBase<double> & b);

        static DenseVector<float> & value(DenseVector<float> & result, const BandedMatrixQx<Q1Type, float> & a, const DenseVectorContinuousBase<float> & b);

        static DenseVector<double> & value(DenseVector<double> & result, const BandedMatrixQx<Q1Type, double> & a, const DenseVectorContinuousBase<double> & b);

        static DenseVector<float> & value(DenseVector<float> & result, const SparseMatrixELL<float> & a, const DenseVector<float> & b);

        static DenseVector<double> & value(DenseVector<double> & result, const SparseMatrixELL<double> & a, const DenseVector<double> & b);

        static DenseVector<float> & value(DenseVector<float> & result, const SparseMatrixCSR<float> & a, const DenseVector<float> & b);

        static DenseVector<double> & value(DenseVector<double> & result, const SparseMatrixCSR<double> & a, const DenseVector<double> & b);

        template<typename DT1_, typename DT2_>
        static DenseVectorContinuousBase<DT1_> & value(DenseVectorContinuousBase<DT1_> & y, const DenseVectorContinuousBase<DT1_> & a, const DenseVectorContinuousBase<DT2_> & b)
        {
            ElementProduct<tags::GPU::CUDA>::value(y, a, b);
            return y;
        }

        template <typename DT_>
        static inline DenseVectorMPI<DT_> & value(DenseVectorMPI<DT_> & r, const SparseMatrixELLMPI<DT_> & a, const DenseVectorMPI<DT_> & b)
        {
            MPIOps<tags::GPU::CUDA>::product(r, a, b);
            return r;
        }

        template <typename DT_>
        static inline DenseVectorMPI<DT_> & value(DenseVectorMPI<DT_> & r, const SparseMatrixCSRMPI<DT_> & a, const DenseVectorMPI<DT_> & b)
        {
            MPIOps<tags::GPU::CUDA>::product(r, a, b);
            return r;
        }

        template<typename DT_>
        static inline DenseVectorMPI<DT_> & value(DenseVectorMPI<DT_> & r, const DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y)
        {
            MPIOps<tags::GPU::CUDA>::element_product(r, x, y);
            return r;
        }
        /// \}
    };

    template <> struct Product<tags::GPU::MultiCore::CUDA>
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

        static DenseVector<float> & value(DenseVector<float> & result, const SparseMatrixELL<float> & a, const DenseVector<float> & b);

        static DenseVector<double> & value(DenseVector<double> & result, const SparseMatrixELL<double> & a, const DenseVector<double> & b);

        template<typename DT1_, typename DT2_>
        static DenseVectorContinuousBase<DT1_> & value(DenseVectorContinuousBase<DT1_> & y, const DenseVectorContinuousBase<DT1_> & a, const DenseVectorContinuousBase<DT2_> & b)
        {
            ElementProduct<tags::GPU::MultiCore::CUDA>::value(y, a, b);
            return y;
        }

        /// \}
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

        static DenseVector<float> value(const BandedMatrixQx<Q1Type, float> & a, const DenseVectorContinuousBase<float> & b);

        static DenseVector<double> value(const BandedMatrixQx<Q1Type, double> & a, const DenseVectorContinuousBase<double> & b);

        static DenseVector<float> & value(DenseVector<float> & result, const BandedMatrixQx<Q1Type, float> & a, const DenseVectorContinuousBase<float> & b);

        static DenseVector<double> & value(DenseVector<double> & result, const BandedMatrixQx<Q1Type, double> & a, const DenseVectorContinuousBase<double> & b);

        static DenseVector<float> value(const DenseMatrix<float> & a, const DenseVectorContinuousBase<float> & b);

        static DenseVector<double> value(const DenseMatrix<double> & a, const DenseVectorContinuousBase<double> & b);

        static DenseMatrix<float> value(const DenseMatrix<float> &a, const DenseMatrix<float> & b);

        static DenseMatrix<double> value(const DenseMatrix<double> &a, const DenseMatrix<double> & b);

        static DenseMatrix<float> value(const SparseMatrix<float> &a, const DenseMatrix<float> & b);

        static DenseMatrix<double> value(const SparseMatrix<double> &a, const DenseMatrix<double> & b);

        static DenseMatrixTile<float> & value(DenseMatrixTile<float> & r, const DenseMatrixTile<float> & a, const DenseMatrixTile<float> & b);

        static DenseMatrixTile<double> & value(DenseMatrixTile<double> & r, const DenseMatrixTile<double> & a, const DenseMatrixTile<double> & b);

        static DenseVector<float> & value(DenseVector<float> & result, const SparseMatrixELL<float> & a, const DenseVector<float> & b,
                unsigned long row_start = 0, unsigned long row_end = 0);

        static DenseVector<double> & value(DenseVector<double> & result, const SparseMatrixELL<double> & a, const DenseVector<double> & b,
                unsigned long row_start = 0, unsigned long row_end = 0);

        static DenseVector<float> & value(DenseVector<float> & result, const SparseMatrixCSR<float> & a, const DenseVector<float> & b,
                unsigned long row_start = 0, unsigned long row_end = 0);

        static DenseVector<double> & value(DenseVector<double> & result, const SparseMatrixCSR<double> & a, const DenseVector<double> & b,
                unsigned long row_start = 0, unsigned long row_end = 0);

        template<typename DT1_, typename DT2_>
        static DenseVectorContinuousBase<DT1_> & value(DenseVectorContinuousBase<DT1_> & y, const DenseVectorContinuousBase<DT1_> & a, const DenseVectorContinuousBase<DT2_> & b)
        {
            ElementProduct<tags::CPU::SSE>::value(y, a, b);
            return y;
        }

        template <typename DT_>
        static inline DenseVectorMPI<DT_> & value(DenseVectorMPI<DT_> & r, const SparseMatrixELLMPI<DT_> & a, const DenseVectorMPI<DT_> & b)
        {
            MPIOps<tags::CPU::SSE>::product(r, a, b);
            return r;
        }

        template <typename DT_>
        static inline DenseVectorMPI<DT_> & value(DenseVectorMPI<DT_> & r, const SparseMatrixCSRMPI<DT_> & a, const DenseVectorMPI<DT_> & b)
        {
            MPIOps<tags::CPU::SSE>::product(r, a, b);
            return r;
        }

        template<typename DT_>
        static inline DenseVectorMPI<DT_> & value(DenseVectorMPI<DT_> & r, const DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y)
        {
            MPIOps<tags::CPU::SSE>::element_product(r, x, y);
            return r;
        }
        /// \}
    };

    template <>
    struct Product<tags::Cell>
    {
        static DenseVector<float> value(const DenseMatrix<float> & a, const DenseVectorContinuousBase<float> & b);
        static DenseVector<double> value(const DenseMatrix<double> & a, const DenseVectorContinuousBase<double> & b);
        static DenseMatrix<float> value(const DenseMatrix<float> & a, const DenseMatrix<float> & b);
        static DenseMatrix<double> value(const DenseMatrix<double> & a, const DenseMatrix<double> & b);

        static DenseVector<float> value(const BandedMatrix<float> & a, const DenseVectorContinuousBase<float> & b);
        static DenseVector<double> value(const BandedMatrix<double> & a, const DenseVectorContinuousBase<double> & b);
        static DenseMatrix<float> value(const SparseMatrix<float> & a, const DenseMatrix<float> & b);

        template<typename DT1_, typename DT2_>
        static DenseVectorContinuousBase<DT1_> & value(DenseVectorContinuousBase<DT1_> & y, const DenseVectorContinuousBase<DT1_> & a, const DenseVectorContinuousBase<DT2_> & b)
        {
            //TODO
            DenseVector<DT1_> temp(a.copy());
            ElementProduct<tags::Cell>::value(temp, b);
            copy<tags::Cell>(temp, y);
            return y;
        }
    };

    namespace mc
    {
        template <typename Tag_> struct Product
        {
            template<typename DT1_, typename DT2_>
            static DenseVectorContinuousBase<DT1_> & value(DenseVectorContinuousBase<DT1_> & y, const DenseVectorContinuousBase<DT1_> & a, const DenseVectorContinuousBase<DT2_> & b)
            {
                ElementProduct<Tag_>::value(y, a, b);
                return y;
            }

            template <typename DT1_, typename DT2_>
                static DenseVector<DT1_> value(const BandedMatrixQx<Q1Type, DT1_> & a, const DenseVectorContinuousBase<DT2_> & b)
            {
                CONTEXT("When multiplying BandedMatrixQ1 with DenseVectorContinuousBase using backend : " + Tag_::name);
                if (b.size() != a.columns())
                {
                    throw VectorSizeDoesNotMatch(b.size(), a.columns());
                }

                signed long root(a.root());
                const unsigned long size(b.size());
                DenseVector<DT1_> result(size, DT1_(0));

                const DenseVectorRange<DT1_> ll_band(a.band_range(LL));
                const DenseVectorRange<DT1_> ld_band(a.band_range(LD));
                const DenseVectorRange<DT1_> lu_band(a.band_range(LU));
                const DenseVectorRange<DT1_> dl_band(a.band_range(DL));
                const DenseVectorRange<DT1_> du_band(a.band_range(DU));
                const DenseVectorRange<DT1_> ud_band(a.band_range(UD));
                const DenseVectorRange<DT1_> uu_band(a.band_range(UU));
                const DenseVectorRange<DT1_> ul_band(a.band_range(UL));

                DenseVectorRange<DT1_> res_ll(result.range(size - root - 1, root + 1));
                DenseVectorRange<DT1_> res_ld(result.range(size - root, root));
                DenseVectorRange<DT1_> res_lu(result.range(size - root + 1, root - 1));
                DenseVectorRange<DT1_> res_dl(result.range(size - 1, 1));
                DenseVectorRange<DT1_> res_du(result.range(size - 1, 0));
                DenseVectorRange<DT1_> res_ud(result.range(size - root, 0));
                DenseVectorRange<DT1_> res_uu(result.range(size - root - 1, 0));
                DenseVectorRange<DT1_> res_ul(result.range(size - root + 1, 0));

                const DenseVectorRange<DT1_> b_ll(b.range(size - root - 1, 0));
                const DenseVectorRange<DT1_> b_ld(b.range(size - root, 0));
                const DenseVectorRange<DT1_> b_lu(b.range(size - root + 1, 0));
                const DenseVectorRange<DT1_> b_dl(b.range(size - 1, 0));
                const DenseVectorRange<DT1_> b_du(b.range(size - 1, 1));
                const DenseVectorRange<DT1_> b_ud(b.range(size - root, root));
                const DenseVectorRange<DT1_> b_uu(b.range(size - root - 1, root + 1));
                const DenseVectorRange<DT1_> b_ul(b.range(size - root + 1, root - 1));

                honei::ScaledSum<Tag_>::value(res_ll, ll_band, b_ll);
                honei::ScaledSum<Tag_>::value(res_ld, ld_band, b_ld);
                honei::ScaledSum<Tag_>::value(res_lu, lu_band, b_lu);
                honei::ScaledSum<Tag_>::value(res_dl, dl_band, b_dl);
                honei::ScaledSum<Tag_>::value(result, const_cast<const DenseVector<DT1_> &>(a.band(DD)), b);
                honei::ScaledSum<Tag_>::value(res_du, du_band, b_du);
                honei::ScaledSum<Tag_>::value(res_ud, ud_band, b_ud);
                honei::ScaledSum<Tag_>::value(res_uu, uu_band, b_uu);
                honei::ScaledSum<Tag_>::value(res_ul, ul_band, b_ul);

                return result;
            }

            template <typename DT1_, typename DT2_>
            static DenseVector<DT1_> & value(DenseVector<DT1_> & result, const BandedMatrixQx<Q1Type, DT1_> & a, const DenseVectorContinuousBase<DT2_> & b)
            {
                CONTEXT("When multiplying BandedMatrixQ1 with DenseVectorContinuousBase using backend : " + Tag_::name);
                if (b.size() != a.columns())
                {
                    throw VectorSizeDoesNotMatch(b.size(), a.columns());
                }

                signed long root(a.root());
                const unsigned long size(b.size());
                fill<Tag_>(result, DT1_(0));

                const DenseVectorRange<DT1_> ll_band(a.band_range(LL));
                const DenseVectorRange<DT1_> ld_band(a.band_range(LD));
                const DenseVectorRange<DT1_> lu_band(a.band_range(LU));
                const DenseVectorRange<DT1_> dl_band(a.band_range(DL));
                const DenseVectorRange<DT1_> du_band(a.band_range(DU));
                const DenseVectorRange<DT1_> ud_band(a.band_range(UD));
                const DenseVectorRange<DT1_> uu_band(a.band_range(UU));
                const DenseVectorRange<DT1_> ul_band(a.band_range(UL));

                DenseVectorRange<DT1_> res_ll(result.range(size - root - 1, root + 1));
                DenseVectorRange<DT1_> res_ld(result.range(size - root, root));
                DenseVectorRange<DT1_> res_lu(result.range(size - root + 1, root - 1));
                DenseVectorRange<DT1_> res_dl(result.range(size - 1, 1));
                DenseVectorRange<DT1_> res_du(result.range(size - 1, 0));
                DenseVectorRange<DT1_> res_ud(result.range(size - root, 0));
                DenseVectorRange<DT1_> res_uu(result.range(size - root - 1, 0));
                DenseVectorRange<DT1_> res_ul(result.range(size - root + 1, 0));

                const DenseVectorRange<DT1_> b_ll(b.range(size - root - 1, 0));
                const DenseVectorRange<DT1_> b_ld(b.range(size - root, 0));
                const DenseVectorRange<DT1_> b_lu(b.range(size - root + 1, 0));
                const DenseVectorRange<DT1_> b_dl(b.range(size - 1, 0));
                const DenseVectorRange<DT1_> b_du(b.range(size - 1, 1));
                const DenseVectorRange<DT1_> b_ud(b.range(size - root, root));
                const DenseVectorRange<DT1_> b_uu(b.range(size - root - 1, root + 1));
                const DenseVectorRange<DT1_> b_ul(b.range(size - root + 1, root - 1));

                honei::ScaledSum<Tag_>::value(res_ll, ll_band, b_ll);
                honei::ScaledSum<Tag_>::value(res_ld, ld_band, b_ld);
                honei::ScaledSum<Tag_>::value(res_lu, lu_band, b_lu);
                honei::ScaledSum<Tag_>::value(res_dl, dl_band, b_dl);
                honei::ScaledSum<Tag_>::value(result, const_cast<const DenseVector<DT1_> &>(a.band(DD)), b);
                honei::ScaledSum<Tag_>::value(res_du, du_band, b_du);
                honei::ScaledSum<Tag_>::value(res_ud, ud_band, b_ud);
                honei::ScaledSum<Tag_>::value(res_uu, uu_band, b_uu);
                honei::ScaledSum<Tag_>::value(res_ul, ul_band, b_ul);

                return result;
            }

            template <typename DT_>
            static DenseVector<DT_> & value(DenseVector<DT_> & result, const SparseMatrix<DT_> & a, const DenseVector<DT_> & b)
            {
                if (b.size() != a.columns())
                {
                    throw VectorSizeDoesNotMatch(b.size(), a.columns());
                }
                if (a.rows() != result.size())
                {
                    throw VectorSizeDoesNotMatch(a.rows(), result.size());
                }

                unsigned long max_count(Configuration::instance()->get_value("mc::Product(DV,SMELL,DV)::max_count",
                            mc::ThreadPool::instance()->num_threads()));

                TicketVector tickets;

                unsigned long limits[max_count + 1];
                limits[0] = 0;
                for (unsigned long i(1) ; i < max_count; ++i)
                {
                    limits[i] = limits[i-1] + a.rows() / max_count;
                }
                limits[max_count] = a.rows();

                for (unsigned long i(0) ; i < max_count ; ++i)
                {
                    OperationWrapper<honei::Product<typename Tag_::DelegateTo>, DenseVector<DT_>,
                        DenseVector<DT_>, SparseMatrix<DT_>, DenseVector<DT_>, unsigned long, unsigned long > wrapper(result);
                    tickets.push_back(mc::ThreadPool::instance()->enqueue(bind(wrapper, result, a, b, limits[i], limits[i+1])));
                }

                tickets.wait();

                return result;
            }

            template <typename DT_>
            static DenseVector<DT_> & value(DenseVector<DT_> & result, const SparseMatrixELL<DT_> & a, const DenseVector<DT_> & b)
            {
                if (b.size() != a.columns())
                {
                    throw VectorSizeDoesNotMatch(b.size(), a.columns());
                }
                if (a.rows() != result.size())
                {
                    throw VectorSizeDoesNotMatch(a.rows(), result.size());
                }

                //fill<typename Tag_::DelegateTo>(result, DT_(0));

                unsigned long max_count(Configuration::instance()->get_value("mc::Product(DV,SMELL,DV)::max_count",
                            mc::ThreadPool::instance()->num_threads()));

                TicketVector tickets;

                unsigned long limits[max_count + 1];
                limits[0] = 0;
                for (unsigned long i(1) ; i < max_count; ++i)
                {
                    limits[i] = limits[i-1] + a.rows() / max_count;
                }
                limits[max_count] = a.rows();

                for (unsigned long i(0) ; i < max_count ; ++i)
                {
                    OperationWrapper<honei::Product<typename Tag_::DelegateTo>, DenseVector<DT_>,
                        DenseVector<DT_>, SparseMatrixELL<DT_>, DenseVector<DT_>, unsigned long, unsigned long > wrapper(result);
                    tickets.push_back(mc::ThreadPool::instance()->enqueue(bind(wrapper, result, a, b, limits[i], limits[i+1])));
                }

                tickets.wait();

                return result;
            }

            template <typename DT_>
            static DenseVector<DT_> & value(DenseVector<DT_> & result, const SparseMatrixCSR<DT_> & a, const DenseVector<DT_> & b)
            {
                CONTEXT("When multiplying SparseMatrixCSR with DenseVector (MC):");
                if (b.size() != a.columns())
                {
                    throw VectorSizeDoesNotMatch(b.size(), a.columns());
                }
                if (a.rows() != result.size())
                {
                    throw VectorSizeDoesNotMatch(a.rows(), result.size());
                }

                //fill<typename Tag_::DelegateTo>(result, DT_(0));

                unsigned long max_count(Configuration::instance()->get_value("mc::Product(DV,SMELL,DV)::max_count",
                            mc::ThreadPool::instance()->num_threads()));

                TicketVector tickets;

                unsigned long limits[max_count + 1];
                limits[0] = 0;
                for (unsigned long i(1) ; i < max_count; ++i)
                {
                    limits[i] = limits[i-1] + a.rows() / max_count;
                }
                limits[max_count] = a.rows();

                for (unsigned long i(0) ; i < max_count ; ++i)
                {
                    OperationWrapper<honei::Product<typename Tag_::DelegateTo>, DenseVector<DT_>,
                        DenseVector<DT_>, SparseMatrixCSR<DT_>, DenseVector<DT_>, unsigned long, unsigned long > wrapper(result);
                    tickets.push_back(mc::ThreadPool::instance()->enqueue(bind(wrapper, result, a, b, limits[i], limits[i+1])));
                }

                tickets.wait();

                return result;
            }

            template <typename DT_>
            static inline DenseVectorMPI<DT_> & value(DenseVectorMPI<DT_> & r, const SparseMatrixELLMPI<DT_> & a, const DenseVectorMPI<DT_> & b)
            {
                MPIOps<Tag_>::product(r, a, b);
                return r;
            }

            template <typename DT_>
            static inline DenseVectorMPI<DT_> & value(DenseVectorMPI<DT_> & r, const SparseMatrixCSRMPI<DT_> & a, const DenseVectorMPI<DT_> & b)
            {
                MPIOps<Tag_>::product(r, a, b);
                return r;
            }

            template<typename DT_>
            static inline DenseVectorMPI<DT_> & value(DenseVectorMPI<DT_> & r, const DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y)
            {
                MPIOps<Tag_>::element_product(r, x, y);
                return r;
            }

            // Dummy
            template <typename DT1_, typename DT2_>
                static DenseVector<DT1_> value(const DenseMatrix<DT1_> & a, const DenseVectorBase<DT2_> & b)
                {
                    CONTEXT("When multiplying DenseMatrix with DenseVector(Base) using backend : " + Tag_::name);

                    if (b.size() != a.columns())
                    {
                        throw VectorSizeDoesNotMatch(b.size(), a.columns());
                    }

                    return honei::Product<tags::CPU>::value(a, b);
                }

            // Dummy
            template <typename DT1_, typename DT2_>
                static DenseVector<DT1_> value(const BandedMatrix<DT1_> & a, const DenseVectorBase<DT2_> & b)
                {
                    CONTEXT("When multiplying BandedMatrix with DenseVector(Base) using backend : " + Tag_::name);

                    if (b.size() != a.columns())
                    {
                        throw VectorSizeDoesNotMatch(b.size(), a.columns());
                    }

                    return honei::Product<tags::CPU>::value(a, b);
                }

            // Dummy
            template <typename DT1_, typename DT2_>
                static DenseMatrix<DT1_> value(const DenseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
                {
                    CONTEXT("When multiplying DenseMatrix with DenseMatrix using backend : " + Tag_::name);

                    if (a.columns() != b.rows())
                        throw MatrixRowsDoNotMatch(b.rows(), a.columns());

                    return honei::Product<tags::CPU>::value(a, b);
                }

            // Dummy
            template <typename DT1_, typename DT2_>
                static DenseMatrix<DT1_> value(const SparseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
                {
                    CONTEXT("When multiplying SparseMatrix with DenseMatrix using backend : " + Tag_::name);

                    if (a.columns() != b.rows())
                        throw MatrixRowsDoNotMatch(b.rows(), a.columns());

                    return honei::Product<tags::CPU>::value(a, b);
                }

        };

    }

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

    template <> struct Product<tags::CPU::MultiCore> :
        public mc::Product<tags::CPU::MultiCore>
        {
        };

    template <> struct Product<tags::CPU::MultiCore::Generic> :
        public mc::Product<tags::CPU::MultiCore::Generic>
        {
        };

    template <> struct Product<tags::CPU::MultiCore::SSE> :
        public mc::Product<tags::CPU::MultiCore::SSE>
        {
        };
}
#endif
