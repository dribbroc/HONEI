/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007, 2008, 2009 Sven Mallach <mallach@honei.org>
 * Copyright (c) 2008 Danny van Dyk <dyk@honei.org>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef LIBLA_GUARD_DIFFERENCE_HH
#define LIBLA_GUARD_DIFFERENCE_HH 1

#include <honei/backends/multicore/operation.hh>
#include <honei/backends/multicore/thread_pool.hh>
#include <honei/la/banded_matrix.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/matrix_error.hh>
#include <honei/la/scale.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/la/dense_vector.hh>
#include <honei/util/benchmark_info.hh>
#include <honei/util/configuration.hh>
#include <honei/util/tags.hh>

namespace honei
{
    template <typename Tag_ = tags::CPU> struct Difference;

    /**
     * \brief Difference of two entities.
     *
     * Difference is the class template for the subtraction operation
     * \f[
     *     \texttt{Difference}(a, b): \quad r \leftarrow a - b,
     * \f]
     * which yields r, the difference of entities a and b.
     *
     * Usually, the return value is the minuend a after modification. However,
     * there are signatures for which b is returned. For these cases a short
     * notice is added.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <> struct Difference<tags::CPU>
    {
        /**
         * \name Matrix differences that return minuend
         * \{
         *
         * \brief Returns the the difference of two given matrices
         *
         * \param a The entity that is the minuend of the operation.
         * \param b The entity that is the subtrahend of the operation.
         *
         * \retval r Will modify the minuend a and return it.
         *
         * \exception MatrixSizeDoesNotMatch is thrown if two banded matrices do not have the same size.
         * \exception MatrixRowsDoNotMatch is thrown if two matrices do not have the same number of rows.
         * \exception MatrixColumnsDoNotMatch is thrown if two matrices do not have the same number of columns.
         * \exception MatrixIsNotSquare is thrown if a RowAccessMatrix's number of rows does not equal its number of columns.
         */

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
        {
            CONTEXT("When subtracting DenseMatrix from DenseMatrix:");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            typename DenseMatrix<DT2_>::ConstElementIterator r(b.begin_elements());
            for (typename DenseMatrix<DT1_>::ElementIterator l(a.begin_elements()),
                    l_end(a.end_elements()) ; l != l_end ; ++l, ++r)
            {
                *l -= *r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
        {
            CONTEXT("When subtracting SparseMatrix from DenseMatrix:");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            typename SparseMatrix<DT2_>::NonZeroConstElementIterator r(b.begin_non_zero_elements()),
                     r_end(b.end_non_zero_elements());
            for ( ; r != r_end ; ++r )
            {
                a(r.row(), r.column()) -= *r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
        {
            CONTEXT("When subtracting SparseMatrix from SparseMatrix:");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            typename SparseMatrix<DT1_>::NonZeroElementIterator l(a.begin_non_zero_elements());
            for (typename SparseMatrix<DT2_>::NonZeroConstElementIterator r(b.begin_non_zero_elements()),
                    r_end(b.end_non_zero_elements()) ; r != r_end ; )
            {
               if (r.index() < l.index())
                {
                    a[r.row()][r.column()] = -(*r);
                    ++r;
                }
                else if (l.index() < r.index())
                {
                    ++l;
                }
                else
                {
                    *l -= *r;
                    ++l; ++r;
                }
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static BandedMatrix<DT1_> & value(BandedMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
        {
            CONTEXT("When subtracting BandedMatrix from BandedMatrix:");

            if (a.size() != b.size())
            {
                throw MatrixSizeDoesNotMatch(b.size(), a.size());
            }

            typename BandedMatrix<DT1_>::BandIterator l(a.begin_bands()), l_end(a.end_bands());
            typename BandedMatrix<DT2_>::ConstBandIterator r(b.begin_bands()), r_end(b.end_bands());
            for ( ; ((l != l_end) && (r != r_end)) ; ++l, ++r)
            {
                if (! r.exists())
                    continue;

                if (l.exists())
                {
                    DenseVector<DT1_> band(*l);
                    Difference<>::value(band, *r);
                }
                else
                {
                    DenseVector<DT2_> band(r->copy());
                    Scale<>::value(band, DT1_(-1));
                    a.band(r.index() - a.size() + 1) = band;
                }
            }

            return a;
        }

        /// \}

        /**
         * \name Matrix differences that return subtrahend
         * \{
         *
         * \brief Returns the the difference of two given matrices.
         *
         * \param a The matrix that is the minuend of the operation.
         * \param b The matrix that is the subtrahend of the operation.
         *
         * \retval b Will modify the subtrahend b and return it.
         *
         * \exception MatrixSizeDoesNotMatch is thrown if two banded matrices do not have the same size.
         * \exception MatrixRowsDoNotMatch is thrown if two matrices do not have the same number of rows.
         * \exception MatrixColumnsDoNotMatch is thrown if two matrices do not have the same number of columns.
         * \exception MatrixIsNotSquare is thrown if a row access matrix's number of rows does not equal its number of columns.
         */

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT2_> & value(const BandedMatrix<DT1_> & a, DenseMatrix<DT2_> & b)
        {
            CONTEXT("When subtracting DenseMatrix from BandedMatrix:");

            if (b.columns() != b.rows())
            {
                throw MatrixIsNotSquare(b.rows(), b.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

#if defined (HONEI_SSE)
            Scale<tags::CPU::SSE>::value(b, -1);
#else
            Scale<tags::CPU>::value(b, -1);
#endif

            unsigned long middle_index(a.rows() -1);
            for (typename BandedMatrix<DT1_>::ConstBandIterator vi(a.begin_non_zero_bands()), vi_end(a.end_non_zero_bands()) ;
                    vi != vi_end ; ++vi)
            {
                // If we are below the diagonal band, we start at Element index and go on until the last element.
                if (vi.index() < middle_index)
                {
                    unsigned long start(middle_index - vi.index()); //Calculation of the element-index to start in iteration!
                    unsigned long i(0);
                    for(typename DenseVector<DT1_>::ConstElementIterator c(vi->element_at(start)),
                            c_end(vi->end_elements()) ; c < c_end ; ++c)
                    {
                        b(start, i) += *c;
                        ++start, ++i;
                    }
                }

                // If we are above or on the diagonal band, we start at Element 0 and go on until Element band_size-band_index.
                else
                {

                    //Calculation of the element-index to stop in iteration!
                    unsigned long offset(vi.index() - middle_index);
                    unsigned long end(vi->size() - offset);
                    unsigned long i(0);
                    for(typename DenseVector<DT1_>::ConstElementIterator c(vi->begin_elements()),
                            c_end(vi->element_at(end)) ; c < c_end ; ++c)
                    {
                        b(i, offset) +=  *c;
                        ++offset, ++i;
                    }
                }
            }

            return b;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT2_> & value(const BandedMatrix<DT1_> & a, SparseMatrix<DT2_> & b)
        {
            CONTEXT("When subtracting SparseMatrix from BandedMatrix:");

            if (b.columns() != b.rows())
            {
                throw MatrixIsNotSquare(b.rows(), b.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            Scale<tags::CPU>::value(b, -1);

            unsigned long  middle_index(a.rows() -1);
            for (typename BandedMatrix<DT2_>::ConstBandIterator vi(a.begin_non_zero_bands()), vi_end(a.end_non_zero_bands()) ;
                    vi != vi_end ; ++vi)
            {
                // If we are below the diagonal band, we start at Element index and go on until the last element.
                if (vi.index() < middle_index)
                {
                    unsigned long start(middle_index - vi.index()); //Calculation of the element-index to start in iteration!
                    unsigned long i(0);
                    for (typename DenseVector<DT1_>::ConstElementIterator c(vi->element_at(start)),
                            c_end(vi->end_elements()) ; c < c_end ; ++c)
                    {
                        b[start][i] += *c;
                        ++start, ++i;
                    }
                }

                // If we are above or on the diagonal band, we start at Element 0 and go on until Element band_size-band_index.
                else
                {
                    //Calculation of the element-index to stop in iteration!
                    unsigned long offset(vi.index() - middle_index);
                    unsigned long end(vi->size() - offset);
                    unsigned long i(0);
                    for(typename DenseVector<DT1_>::ConstElementIterator c(vi->begin_elements()),
                            c_end(vi->element_at(end)) ; c < c_end ; ++c)
                    {
                        b[i][offset] +=  *c;
                        ++offset, ++i;
                    }
                }
            }

            return b;
        }

        /// \}

        /**
         * \name Vector differences
         * \{
         *
         * \brief Returns the the difference of two given vectors.
         *
         * \param a The vector that is the left-hand minuend of the operation.
         * \param b The vector that is the right-hand subtrahend of the operation.
         *
         * \retval a Will normally modify the minuend a and return it. Only for
         *           Difference(SparseVector, DenseVector) the subtrahend b is modified and returned.
         *
         * \exception VectorSizeDoesNotMatch is thrown if the two vectors don't have the same size.
         */

        template <typename DT1_, typename DT2_>
        static DenseVectorBase<DT1_> & value(DenseVectorBase<DT1_> & a, const DenseVectorBase<DT2_> & b)
        {
            CONTEXT("When subtracting DenseVectorBase from DenseVectorBase:");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());

            typename DenseVectorBase<DT2_>::ConstElementIterator r(b.begin_elements());
            for (typename DenseVectorBase<DT1_>::ElementIterator l(a.begin_elements()),
                    l_end(a.end_elements()) ; l != l_end ; ++l)
            {
                *l -= *r;
                ++r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static inline DenseVectorContinuousBase<DT1_> & value(DenseVectorContinuousBase<DT1_> & a, const DenseVectorBase<DT2_> & b)
        {
            CONTEXT("When subtracting DenseVectorBase from DenseVectorContinuousBase:");

            DenseVectorBase<DT1_> & temp = a;
            Difference<>::value(temp, b);
            return a;
        }

        template <typename DT1_, typename DT2_>
        static SparseVector<DT1_> & value(SparseVector<DT1_> & a, const SparseVector<DT2_> & b)
        {
            CONTEXT("When subtracting SparseVector from SparseVector:");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());

            typename SparseVector<DT1_>::NonZeroElementIterator l(a.begin_non_zero_elements());
            for (typename SparseVector<DT2_>::NonZeroConstElementIterator r(b.begin_non_zero_elements()),
                    r_end(b.end_non_zero_elements()) ; r != r_end ; )
            {
                if (r.index() < l.index())
                {
                    a[r.index()] = -(*r);
                    ++r;
                }
                else if (l.index() < r.index())
                {
                    ++l;
                }
                else
                {
                    *l -= *r;
                    ++l; ++r;
                }
            }
            return a;
            ///\todo: perhaps sparsify - i.e. substraction of 7 and 7 possible.
        }

        template <typename DT1_, typename DT2_>
        static DenseVectorBase<DT1_> & value(DenseVectorBase<DT1_> & a, const SparseVector<DT2_> & b)
        {
            CONTEXT("When subtracting SparseVector from DenseVectorBase:");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());

            for (typename SparseVector<DT2_>::NonZeroConstElementIterator r(b.begin_non_zero_elements()),
                    r_end(b.end_non_zero_elements()) ; r != r_end ; ++r)
            {
                a[r.index()] -= *r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static inline DenseVector<DT1_> & value(DenseVector<DT1_> & a, const SparseVector<DT2_> & b)
        {
            DenseVectorBase<DT1_> & temp = a;
            Difference<>::value(temp, b);
            return a;
        }

        template <typename DT1_, typename DT2_>
        static inline DenseVectorRange<DT1_> & value(DenseVectorRange<DT1_> & a, const SparseVector<DT2_> & b)
        {
            DenseVectorBase<DT1_> & temp = a;
            Difference<>::value(temp, b);
            return a;
        }

        template <typename DT1_, typename DT2_>
        static inline DenseVectorSlice<DT1_> & value(DenseVectorSlice<DT1_> & a, const SparseVector<DT2_> & b)
        {
            DenseVectorBase<DT1_> & temp = a;
            Difference<>::value(temp, b);
            return a;
        }

        template <typename DT1_, typename DT2_>
        static DenseVectorBase<DT2_> & value(const SparseVector<DT1_> & a, DenseVectorBase<DT2_> & b)
        {
            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());

            typename SparseVector<DT1_>::ConstElementIterator l(a.begin_elements());
            for (typename DenseVectorBase<DT2_>::ElementIterator r(b.begin_elements()),
                    r_end(b.end_elements()) ; r != r_end ; ++r , ++l)
            {
                *r = *l - *r;
            }

            return b;
        }

        template <typename DT1_, typename DT2_>
        static inline DenseVector<DT2_> & value(const SparseVector<DT1_> & a, DenseVector<DT2_> & b)
        {
            DenseVectorBase<DT2_> & temp = b;
            Difference<>::value(a, temp);
            return b;
        }

        template <typename DT1_, typename DT2_>
        static inline DenseVectorRange<DT2_> & value(const SparseVector<DT1_> & a, DenseVectorRange<DT2_> & b)
        {
            DenseVectorBase<DT2_> & temp = b;
            Difference<>::value(a, temp);
            return b;
        }

        template <typename DT1_, typename DT2_>
        static inline DenseVectorSlice<DT2_> & value(const SparseVector<DT1_> & a, DenseVectorSlice<DT2_> & b)
        {
            DenseVectorBase<DT2_> & temp = b;
            Difference<>::value(a, temp);
            return b;
        }
        /// \}

        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(BandedMatrix<DT1_> & a, SparseMatrix<DT2_> & b)
        {
            BenchmarkInfo result;
            result.flops = a.rows() * a.columns();
            result.load = a.rows() * a.columns() * (sizeof(DT1_) + sizeof(DT2_));
            result.store = a.rows() * a.columns() * sizeof(DT1_);
            result.size.push_back(a.rows() * a.columns());
            result.size.push_back(b.rows() * b.columns());
            return result; 
        }
        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(BandedMatrix<DT1_> & a, DenseMatrix<DT2_> & b)
        {
            BenchmarkInfo result;
/*#if defined (HONEI_SSE)
            Scale<tags::CPU::SSE>::value(b, -1);
#else
            Scale<tags::CPU>::value(b, -1);
#endif*/

            result = result + Scale<>::get_benchmark_info(b, DT2_(-1));
            unsigned long middle_index(a.rows() -1);
            for (typename BandedMatrix<DT1_>::ConstBandIterator vi(a.begin_non_zero_bands()), vi_end(a.end_non_zero_bands()) ;
                    vi != vi_end ; ++vi)
            {
                // If we are below the diagonal band, we start at Element index and go on until the last element.
                if (vi.index() < middle_index)
                {
                    unsigned long start(middle_index - vi.index()); //Calculation of the element-index to start in iteration!
                    unsigned long i(0);
                    for (typename DenseVector<DT1_>::ConstElementIterator c(vi->element_at(start)),
                            c_end(vi->end_elements()) ; c < c_end ; ++c)
                    {
                        //b(start, i) += *c;
                        result.flops += 1;
                        result.load += sizeof(DT1_) + sizeof(DT2_);
                        result.store += sizeof(DT2_);
                        ++start, ++i;
                    }
                }

                // If we are above or on the diagonal band, we start at Element 0 and go on until Element band_size-band_index.
                else
                {

                    //Calculation of the element-index to stop in iteration!
                    unsigned long offset(vi.index() - middle_index);
                    unsigned long end(vi->size() - offset);
                    unsigned long i(0);
                    for (typename DenseVector<DT1_>::ConstElementIterator c(vi->begin_elements()),
                            c_end(vi->element_at(end)) ; c < c_end ; ++c)
                    {
                        //b(i, offset) +=  *c;
                        result.flops += 1;
                        result.load += sizeof(DT1_) + sizeof(DT2_);
                        result.store += sizeof(DT2_);
                        ++offset, ++i;
                    }
                }
            }
            result.size.push_back(a.size() * a.size());
            result.size.push_back(b.rows() * b.columns());
            return result; 
        }

        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(DenseVector<DT1_> & a, DenseVector<DT2_> & b)
        {
            BenchmarkInfo result;
            result.flops = a.size();
            result.load = a.size() * (sizeof(DT1_) + sizeof(DT2_));
            result.store = a.size() * sizeof(DT1_);
            result.size.push_back(a.size());
            result.size.push_back(b.size());
            return result; 
        }
    };

    /**
     * \brief Difference of two entities.
     *
     * Difference is the class template for the subtraction operation
     * \f[
     *     \texttt{Difference}(a, b): \quad r \leftarrow a - b,
     * \f]
     * which yields r, the difference of entities a and b.
     *
     * Usually, the return value is the minuend a after modification. However,
     * there are signatures for which b is returned. For these cases a short
     * notice is added.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <>
    struct Difference<tags::GPU::CUDA>
    {
        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b);

        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b);

        static DenseMatrix<float> & value(DenseMatrix<float> & a, const DenseMatrix<float> & b);
    };

    /**
     * \brief Difference of two entities.
     *
     * Difference is the class template for the subtraction operation
     * \f[
     *     \texttt{Difference}(a, b): \quad r \leftarrow a - b,
     * \f]
     * which yields r, the difference of entities a and b.
     *
     * Usually, the return value is the minuend a after modification. However,
     * there are signatures for which b is returned. For these cases a short
     * notice is added.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <>
    struct Difference<tags::CPU::SSE>
    {
        /**
         * \name Vector differences
         * \{
         *
         * \brief Returns the the difference of two given vectors.
         *
         * \param a The vector that is the left-hand minuend of the operation.
         * \param b The vector that is the right-hand subtrahend of the operation.
         *
         * \retval a Will normally modify the minuend a and return it. Only for
         *           Difference(SparseVector, DenseVector) the subtrahend b is modified and returned.
         *
         * \exception VectorSizeDoesNotMatch is thrown if the two vectors don't have the same size.
         */

        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b);

        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b);

        static DenseMatrix<float> & value(DenseMatrix<float> & a, const DenseMatrix<float> & b);

        static DenseMatrix<double> & value(DenseMatrix<double> & a, const DenseMatrix<double> & b);

        static SparseMatrix<float> value(const BandedMatrix<float> &a, SparseMatrix<float> & b)
        {
            CONTEXT("When subtracting SparseMatrix from BandedMatrix (SSE forwarding to CPU):");
            return Difference<tags::CPU>::value(a, b);
        }

        static SparseMatrix<double> value(const BandedMatrix<double> &a, SparseMatrix<double> & b)
        {
            CONTEXT("When subtracting SparseMatrix from BandedMatrix (SSE forwarding to CPU):");
            return Difference<tags::CPU>::value(a, b);
        }

        static DenseMatrix<float> value(const BandedMatrix<float> &a, DenseMatrix<float> & b)
        {
            CONTEXT("When subtracting DenseMatrix from BandedMatrix (SSE forwarding to CPU):");
            return Difference<tags::CPU>::value(a, b);
        }

        static DenseMatrix<double> value(const BandedMatrix<double> &a, DenseMatrix<double> & b)
        {
            CONTEXT("When subtracting DenseMatrix from BandedMatrix (SSE forwarding to CPU):");
            return Difference<tags::CPU>::value(a, b);
        }

        /// \}
    };

   /**
     * \brief Difference of two entities.
     *
     * Difference is the class template for the subtraction operation
     * \f[
     *     \texttt{Difference}(a, b): \quad r \leftarrow a - b,
     * \f]
     * which yields r, the difference of entities a and b.
     *
     * Usually, the return value is the minuend a after modification. However,
     * there are signatures for which b is returned. For these cases a short
     * notice is added.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <>
    struct Difference<tags::Cell>
    {
        /**
         * \name Vector differences
         * \{
         *
         * \brief Returns the the difference of two given vectors.
         *
         * \param a The vector that is the left-hand minuend of the operation.
         * \param b The vector that is the right-hand subtrahend of the operation.
         *
         * \retval a Will normally modify the minuend a and return it. Only for
         *           Difference(SparseVector, DenseVector) the subtrahend b is modified and returned.
         *
         * \exception VectorSizeDoesNotMatch is thrown if the two vectors don't have the same size.
         */

        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & a,
                const DenseVectorContinuousBase<float> & b);

        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & a,
                const DenseVectorContinuousBase<double> & b);

        static DenseMatrix<float> & value(DenseMatrix<float> & a, const DenseMatrix<float> & b);

        static DenseVector<float> & value(const SparseVector<float> & a, DenseVector<float> & b);

        // Dummy
        static DenseMatrix <float> & value(const BandedMatrix<float> & a, DenseMatrix<float> & b)
        {
            Difference<tags::CPU>::value(a, b);

            return b;
        }
        /// \}
    };

    /**
     * \brief Difference of two entities.
     *
     * Difference is the class template for the subtraction operation
     * \f[
     *     \texttt{Difference}(a, b): \quad r \leftarrow a - b,
     * \f]
     * which yields r, the difference of entities a and b.
     *
     * Usually, the return value is the minuend a after modification. However,
     * there are signatures for which b is returned. For these cases a short
     * notice is added.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */

    namespace mc
    {
        template <typename Tag_> struct Difference
        {
            template <typename DT1_, typename DT2_>
            static DenseVectorBase<DT1_> & value(DenseVectorBase<DT1_> & x, const DenseVectorBase<DT2_> & y)
            {
                CONTEXT("When calculating Difference (DenseVectorBase, DenseVectorBase) using backend : " + Tag_::name);
                unsigned long min_part_size(Configuration::instance()->get_value("mc::Difference(DVB,DVB)::min_part_size", 128));
                unsigned long max_count(Configuration::instance()->get_value("mc::Difference(DVB,DVB)::max_count",
                            mc::ThreadPool::instance()->num_threads()));

                Operation<honei::Difference<typename Tag_::DelegateTo> >::op(x, y, min_part_size, max_count);

                return x;
            }

            template <typename DT1_, typename DT2_>
            static DenseVectorContinuousBase<DT1_> & value(DenseVectorContinuousBase<DT1_> & x, const DenseVectorContinuousBase<DT2_> & y)
            {
                CONTEXT("When calculating Difference (DenseVectorContinuousBase, DenseVectorContinuousBase) using backend : " + Tag_::name);

                unsigned long min_part_size(Configuration::instance()->get_value("mc::Difference(DVCB,DVCB)::min_part_size", 128));
                unsigned long max_count(Configuration::instance()->get_value("mc::Difference(DVCB,DVCB)::max_count",
                            mc::ThreadPool::instance()->num_threads()));

                Operation<honei::Difference<typename Tag_::DelegateTo> >::op(x, y, min_part_size, max_count);

                return x;
            }

            // Dummy
            template <typename DT1_, typename DT2_>
            static DenseVectorBase<DT1_> & value(DenseVectorBase<DT1_> & a, const SparseVector<DT2_> & b)
            {
                CONTEXT("When subtracting SparseVector from DenseVectorBase using backend : " + Tag_::name);

                if (a.size() != b.size())
                    throw VectorSizeDoesNotMatch(b.size(), a.size());

                return honei::Difference<tags::CPU>::value(a, b);
            }

            // Dummy
            template <typename DT1_, typename DT2_>
            static SparseVector<DT1_> & value(SparseVector<DT1_> & a, const DenseVectorBase<DT2_> & b)
            {
                CONTEXT("When calculating Difference (SparseVector, DenseVectorBase) using backend : " + Tag_::name);

                if (a.size() != b.size())
                    throw VectorSizeDoesNotMatch(b.size(), a.size());

                return honei::Difference<tags::CPU>::value(a, b);
             }

            // Dummy
            template <typename DT1_, typename DT2_>
            static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
            {
                CONTEXT("When subtracting SparseMatrix from DenseMatrix using backend : " + Tag_::name);

                if (a.columns() != b.columns())
                {
                    throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
                }

                if (a.rows() != b.rows())
                {
                    throw MatrixRowsDoNotMatch(b.rows(), a.rows());
                }
                return honei::Difference<tags::CPU>::value(a, b);
             }

            // Dummy
            template <typename DT1_, typename DT2_>
            static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
            {
                CONTEXT("When calculating Difference (SparseMatrix, SparseMatrix) using backend : " + Tag_::name);

                if (a.columns() != b.columns())
                {
                    throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
                }

                if (a.rows() != b.rows())
                {
                    throw MatrixRowsDoNotMatch(b.rows(), a.rows());
                }

                return honei::Difference<tags::CPU>::value(a, b);
             }

            // Dummy
            template <typename DT1_, typename DT2_>
            static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
            {
                CONTEXT("When calculating Difference (DenseMatrix, DenseMatrix) using backend : " + Tag_::name);

                if (a.columns() != b.columns())
                {
                    throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
                }

                if (a.rows() != b.rows())
                {
                    throw MatrixRowsDoNotMatch(b.rows(), a.rows());
                }

                return honei::Difference<tags::CPU>::value(a, b);
            }

            // Dummy
            template <typename DT1_, typename DT2_>
            static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
            {
                CONTEXT("When calculating Difference (SparseMatrix, DenseMatrix) using backend : " + Tag_::name);

                if (a.columns() != b.columns())
                {
                    throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
                }

                if (a.rows() != b.rows())
                {
                    throw MatrixRowsDoNotMatch(b.rows(), a.rows());
                }

                return honei::Difference<tags::CPU>::value(a, b);
            }

            // Dummy
            template <typename DT1_, typename DT2_>
            static BandedMatrix<DT1_> & value(BandedMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
            {
                CONTEXT("When calculating Difference (BandedMatrix, BandedMatrix) using backend : " + Tag_::name);

                if (a.rows() != b.rows())
                {
                    throw MatrixSizeDoesNotMatch(b.rows(), a.rows());
                }

                return honei::Difference<tags::CPU>::value(a, b);
            }

            // Dummy
            template <typename DT1_, typename DT2_>
            static SparseMatrix<DT2_> & value(const BandedMatrix<DT1_> & a, SparseMatrix<DT2_> & b)
            {
                CONTEXT("When subtracting SparseMatrix from BandedMatrix using backend : " + Tag_::name);

                if (b.columns() != b.rows())
                {
                    throw MatrixIsNotSquare(b.rows(), b.columns());
                }

                if (a.rows() != b.rows())
                {
                    throw MatrixRowsDoNotMatch(b.rows(), a.rows());
                }

                return honei::Difference<tags::CPU>::value(a, b);
             }

            // Dummy
            template <typename DT1_, typename DT2_>
            static DenseMatrix<DT2_> & value(const BandedMatrix<DT1_> & a, DenseMatrix<DT2_> & b)
            {
                CONTEXT("When subtracting DenseMatrix from BandedMatrix using backend : " + Tag_::name);

                if (b.columns() != b.rows())
                {
                    throw MatrixIsNotSquare(b.rows(), b.columns());
                }

                if (a.rows() != b.rows())
                {
                    throw MatrixRowsDoNotMatch(b.rows(), a.rows());
                }

                return honei::Difference<tags::CPU>::value(a, b);
             }

        };
    }

    template <> struct Difference<tags::CPU::MultiCore> :
        public mc::Difference<tags::CPU::MultiCore>
    {
    };

    template <> struct Difference<tags::CPU::MultiCore::SSE> :
        public mc::Difference<tags::CPU::MultiCore::SSE>
    {
    };
}
#endif
