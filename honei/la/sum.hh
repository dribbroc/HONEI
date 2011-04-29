/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007, 2008, 2009 Sven Mallach <mallach@honei.org>
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

#pragma once
#ifndef LIBLA_GUARD_SUM_HH
#define LIBLA_GUARD_SUM_HH 1

#include <honei/backends/multicore/operation.hh>
#include <honei/backends/multicore/thread_pool.hh>
#include <honei/la/banded_matrix.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/matrix_error.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/util/benchmark_info.hh>
#include <honei/util/configuration.hh>
#include <honei/util/operation_wrapper.hh>
#include <honei/util/partitioner.hh>
#include <honei/util/tags.hh>

#include <algorithm>

namespace honei
{
    template <typename Tag_ = tags::CPU> struct Sum;

    /**
     * \brief Sum of two entities
     *
     * Sum is the class template for the addition operation
     * \f[
     *     \texttt{Sum}(a, b): \quad r \leftarrow a + b,
     * \f]
     * which yields r, the sum of entities a and b.
     *
     * The return value is the summand a after modification.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <> struct Sum<tags::CPU>
    {
        /**
         * \name Sums of two matrices
         * \{
         *
         * \brief Returns the the sum of two given entities.
         *
         * \param a The entity that is the left-hand summand of the operation.
         * \param b The entity that is the right-hand summand of the operation.
         *
         * \retval a Will modify the summand a and return it.
         *
         * \exception MatrixSizeDoesNotMatch is thrown if two banded matrices do not have the same size.
         * \exception MatrixRowsDoNotMatch is thrown if two matrices do not have the same number of rows.
         * \exception MatrixColumnsDoNotMatch is thrown if two matrices do not have the same number of columns.
         * \exception MatrixIsNotSquare is thrown if a row access matrix's number of rows does not equal its number of columns.
         * \exception VectorSizeDoesNotMatch is thrown if the two vectors don't have the same size.
         */

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
        {
            CONTEXT("When adding DenseMatrix to DenseMatrix:");

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
                *l += *r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
        {
            CONTEXT("When adding SparseMatrix to DenseMatrix:");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            for (typename SparseMatrix<DT2_>::NonZeroConstElementIterator r(b.begin_non_zero_elements()),
                    r_end(b.end_non_zero_elements()) ; r != r_end ; ++r)
            {
                a(r.row(), r.column()) += *r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
        {
            CONTEXT("When adding SparseMatrix to SparseMatrix:");

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
                    a[r.row()][r.column()] = (*r);
                    ++r;
                }
                else if (l.index() < r.index())
                {
                    ++l;
                }
                else
                {
                    *l += *r;
                    ++l; ++r;
                }
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static BandedMatrix<DT1_> & value(BandedMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
        {
            CONTEXT("When adding BandedMatrix to BandedMatrix:");

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
                    Sum<>::value(band, *r);
                }
                else
                    a.band(r.index()) = r->copy();
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
        {
            CONTEXT("When adding BandedMatrix to DenseMatrix:");

            if (a.columns() != a.rows())
            {
                throw MatrixIsNotSquare(a.rows(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            for (typename BandedMatrix<DT2_>::ConstBandIterator r(b.begin_non_zero_bands()), r_end(b.end_non_zero_bands()) ;
                    r != r_end ; ++r)
            {
                unsigned long size(b.size());
                unsigned long row_index(std::max(long(-(r.index() - size + 1)), long(0)));
                unsigned long col_index(std::max(long(r.index() - size + 1), long(0)));

                typename DenseVector<DT2_>::ConstElementIterator c(r->begin_elements()), c_end(r->end_elements());

                if (r.index() < size - 1)
                {
                        c += ((size-1) - r.index());
                }

                for ( ; c != c_end ; ++c)
                {
                    if (row_index >= size)
                        break;

                    if (col_index >= size)
                        break;

                    a(row_index, col_index) += *c;
                    ++row_index;
                    ++col_index;
                }
            }
            return a;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
        {
            CONTEXT("When adding BandedMatrix to SparseMatrix:");

            if (a.columns() != a.rows())
            {
                throw MatrixIsNotSquare(a.rows(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            typename SparseMatrix<DT1_>::NonZeroElementIterator l(a.begin_non_zero_elements());
            for (typename BandedMatrix<DT2_>::ConstElementIterator r(b.begin_elements()), r_end(b.end_elements()) ;
                    r != r_end ; ++r)
            {
                if (r.index() < l.index())
                {
                    a[r.row()][r.column()] = *r;
                }
                else // l.index() < r.index() not possible
                {
                *l += *r;
                ++l;
                }
            }

            return a;
        }

        /// \}

        /**
         * \name Sums of scalar and matrix
         * \{
         *
         * Returns the matrix after adding a scalar to every element.
         *
         * \param a The DenseMatrix to be used.
         * \param b The scalar to be added.
         *
         * \retval a The referenced matrix is changed by adding the given scalar to each of its elements.
         */

        template <typename DT_>
        static DenseMatrix<DT_> & value(DenseMatrix<DT_> & x, const DT_ a)
        {
            CONTEXT("When adding scalar to DenseMatrix:");

            for (typename DenseMatrix<DT_>::ElementIterator l(x.begin_elements()),
                    l_end(x.end_elements()) ; l != l_end ; ++l)
            {
                *l += a;
            }

            return x;
        }

        /// \}

        /**
         * \name Sums of two vectors
         * \{
         *
         * \brief Returns the the sum of two given vectors.
         *
         * \param a The vector that is the left-hand summand of the operation.
         * \param b The vector that is the right-hand summand of the operation.
         *
         * \retval Will modify the summand a and return it.
         *
         * \exception VectorSizeDoesNotMatch is thrown if the two vectors don't have the same size.
         */

        template <typename DT1_, typename DT2_>
        static DenseVectorBase<DT1_> & value(DenseVectorBase<DT1_> & a, const DenseVectorBase<DT2_> & b)
        {
            CONTEXT("When adding DenseVectorBase to DenseVectorBase:");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());

            typename DenseVectorBase<DT2_>::ConstElementIterator r(b.begin_elements());
            for (typename DenseVectorBase<DT1_>::ElementIterator l(a.begin_elements()),
                    l_end(a.end_elements()) ; l != l_end ; ++l)
            {
                *l += *r;
                ++r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static inline DenseVectorContinuousBase<DT1_> & value(DenseVectorContinuousBase<DT1_> & a, const DenseVectorBase<DT2_> & b)
        {
            DenseVectorBase<DT1_> & temp = a;
            Sum<tags::CPU>::value(temp, b);
            return a;
        }

        template <typename DT1_, typename DT2_>
        static SparseVector<DT1_> & value(SparseVector<DT1_> & a, const SparseVector<DT2_> & b)
        {
            CONTEXT("When adding SparseVector to SparseVector:");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());

            typename SparseVector<DT1_>::NonZeroElementIterator l(a.begin_non_zero_elements());
            for (typename SparseVector<DT2_>::NonZeroConstElementIterator r(b.begin_non_zero_elements()),
                    r_end(b.end_non_zero_elements()) ; r != r_end ; )
            {
                if (r.index() < l.index())
                {
                    a[r.index()] = *r;
                    ++r;
                }
                else if (l.index() < r.index())
                {
                    ++l;
                }
                else
                {
                    *l += *r;
                    ++l; ++r;
                }
            }
            ///\todo: perhaps sparsify - i.e. addition of -7 and 7 possible
            return a;
        }

        template <typename DT1_, typename DT2_>
        static DenseVectorBase<DT1_> & value(DenseVectorBase<DT1_> & a, const SparseVector<DT2_> & b)
        {
            CONTEXT("When adding DenseVectorBase to SparseVector:");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());

            for (typename SparseVector<DT2_>::NonZeroConstElementIterator r(b.begin_non_zero_elements()),
                    r_end(b.end_non_zero_elements()) ; r != r_end ; ++r )
            {
                a[r.index()] += *r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static inline DenseVectorContinuousBase<DT1_> & value(DenseVectorContinuousBase<DT1_> & a, const SparseVector<DT2_> & b)
        {
            DenseVectorBase<DT1_> & temp = a;
            Sum<>::value(temp, b);
            return a;
        }

        /// \}

        /**
         * \name Sums of scalar and vector
         * \{
         *
         * Returns the vector after adding a scalar to every element.
         *
         * \param x The Vector that shall be added.
         * \param a The scalar that shall be added.
         *
         * \retval x Will return x after modification.
         */

        template <typename DT_>
        static DenseVectorBase<DT_> & value(DenseVectorBase<DT_> & x, const DT_ a)
        {
            CONTEXT("When adding scalar to DenseVectorBase:");

            for (typename DenseVectorBase<DT_>::ElementIterator l(x.begin_elements()),
                    l_end(x.end_elements()) ; l != l_end ; ++l)
            {
                *l += a;
            }

            return x;
        }

        template <typename DT_>
        static inline DenseVectorContinuousBase<DT_> & value(DenseVectorContinuousBase<DT_> & x, const DT_ a)
        {
            DenseVectorBase<DT_> & temp = x;
            Sum<>::value(temp, a);
            return x;
        }

        /// \}

        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(const DenseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
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
        static inline BenchmarkInfo get_benchmark_info(HONEI_UNUSED const DenseMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
        {
            BenchmarkInfo result;
            for (typename BandedMatrix<DT2_>::ConstBandIterator r(b.begin_non_zero_bands()), r_end(b.end_non_zero_bands()) ;
                    r != r_end ; ++r)
            {
                unsigned long size(b.size());
                unsigned long row_index(std::max(long(-(r.index() - size + 1)), long(0)));
                unsigned long col_index(std::max(long(r.index() - size + 1), long(0)));

                typename DenseVector<DT2_>::ConstElementIterator c(r->begin_elements()), c_end(r->end_elements());

                if (r.index() < size - 1)
                {
                        c += ((size-1) - r.index());
                }

                for ( ; c != c_end ; ++c)
                {
                    if (row_index >= size)
                        break;

                    if (col_index >= size)
                        break;

                    //a(row_index, col_index) += *c;
                    result.flops += 1;
                    result.load += sizeof(DT1_) + sizeof(DT2_);
                    result.store += sizeof(DT1_);
                    ++row_index;
                    ++col_index;
                }
            }
            return result;
        }

        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(const DenseVectorBase<DT1_> & a, const DenseVectorBase<DT2_> & b)
        {
            BenchmarkInfo result;
            result.flops = a.size();
            result.load = a.size() * (sizeof(DT1_) + sizeof(DT2_));
            result.store = a.size() * sizeof(DT1_);
            result.size.push_back(a.size());
            result.size.push_back(b.size());
            return result;
        }

        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(const DenseMatrix<DT1_> & a, HONEI_UNUSED DT2_ b)
        {
            BenchmarkInfo result;
            result.flops = a.rows() * a.columns();
            result.load = a.rows() * a.columns() * sizeof(DT1_) + sizeof(DT2_);
            result.store = a.rows() * a.columns() * sizeof(DT1_);
            result.size.push_back(a.rows() * a.columns());
            return result;
        }

        template <typename DT1_>
        static inline BenchmarkInfo get_benchmark_info(const DenseVectorBase<DT1_> & a, HONEI_UNUSED float b)
        {
            BenchmarkInfo result;
            result.flops = a.size();
            result.load = a.size() * sizeof(DT1_) + sizeof(float);
            result.store = a.size() * sizeof(DT1_);
            result.size.push_back(a.size());
            return result;
        }

        template <typename DT1_>
        static inline BenchmarkInfo get_benchmark_info(const DenseVectorBase<DT1_> & a, HONEI_UNUSED double b)
        {
            BenchmarkInfo result;
            result.flops = a.size();
            result.load = a.size() * sizeof(DT1_) + sizeof(double);
            result.store = a.size() * sizeof(DT1_);
            result.size.push_back(a.size());
            return result;
        }
    };

    /**
     * \brief Sum of two entities
     *
     * Sum is the class template for the addition operation
     * \f[
     *     \texttt{Sum}(a, b): \quad r \leftarrow a + b,
     * \f]
     * which yields r, the sum of entities a and b.
     *
     * The return value is the summand a after modification.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <> struct Sum<tags::Cell>
    {
        /**
         * \name Sums of two matrices.
         * \{
         *
         * Returns the first matrix after adding each element of the second one.
         *
         * \param a The DenseMatrix to be used.
         * \param b The DenseMatrix to be added.
         *
         * \retval a The referenced matrix is changed by adding each corresponding element to each of its elements.
         */

        static DenseMatrix<float> & value(DenseMatrix<float> & a, const DenseMatrix<float> & b);

        static DenseMatrix<double> & value(DenseMatrix<double> & a, const DenseMatrix<double> & b);

        // Dummy
        static DenseMatrix<float> & value(DenseMatrix<float> & a, const BandedMatrix<float> & b)
        {
            Sum<tags::CPU>::value(a, b);
            return a;
        }

        /// \}


        /**
         * \name Sums of scalar and matrix
         * \{
         *
         * Returns the matrix after adding a scalar to every element.
         *
         * \param a The DenseMatrix to be used.
         * \param b The scalar to be added.
         *
         * \retval a The referenced matrix is changed by adding the given scalar to each of its elements.
         */

        static DenseMatrix<float> & value(DenseMatrix<float> & a, const float & b);

        static DenseMatrix<double> & value(DenseMatrix<double> & a, const double & b);

        /// \}

        /**
         * \name Sums of two vectors
         * \{
         *
         * \brief Returns the the sum of two given vectors.
         *
         * \param a The vector that is the left-hand summand of the operation.
         * \param b The vector that is the right-hand summand of the operation.
         *
         * \retval Will modify the summand a and return it.
         *
         * \exception VectorSizeDoesNotMatch is thrown if the two vectors don't have the same size.
         */

        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b);

        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b);

        static DenseVector<float> & value(DenseVector<float> & a, const SparseVector<float> & b);

        /// \}

        /**
         * \name Sums of scalar and vector
         * \{
         *
         * Returns the vector after adding a scalar to every element.
         *
         * \param a The Vector that shall be added.
         * \param b The scalar that shall be added.
         *
         * \retval Will modify the summand a and return it.
         */

        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & a, const float & b);

        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & a, const double & b);

        /// \}
    };

    /**
     * \brief Sum of two entities
     *
     * Sum is the class template for the addition operation
     * \f[
     *     \texttt{Sum}(a, b): \quad r \leftarrow a + b,
     * \f]
     * which yields r, the sum of entities a and b.
     *
     * The return value is the summand a after modification.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <> struct Sum<tags::GPU::CUDA>
    {
        /**
         * \name Sums of two matrices.
         * \{
         *
         * Returns the first matrix after adding each element of the second one.
         *
         * \param a The DenseMatrix to be used.
         * \param b The DenseMatrix to be added.
         *
         * \retval a The referenced matrix is changed by adding each corresponding element to each of its elements.
         */

        static DenseMatrix<float> & value(DenseMatrix<float> & a, const DenseMatrix<float> & b);

        //static DenseMatrix<double> & value(DenseMatrix<double> & a, const DenseMatrix<double> & b);

        /**
         * \name Sums of two vectors
         * \{
         *
         * \brief Returns the the sum of two given vectors.
         *
         * \param a The vector that is the left-hand summand of the operation.
         * \param b The vector that is the right-hand summand of the operation.
         *
         * \retval Will modify the summand a and return it.
         *
         * \exception VectorSizeDoesNotMatch is thrown if the two vectors don't have the same size.
         */

        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b);

        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b);
        /// \}

    };

    /**
     * \brief Sum of two entities
     *
     * Sum is the class template for the addition operation
     * \f[
     *     \texttt{Sum}(a, b): \quad r \leftarrow a + b,
     * \f]
     * which yields r, the sum of entities a and b.
     *
     * The return value is the summand a after modification.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <> struct Sum<tags::GPU::MultiCore::CUDA>
    {
        /**
         * \name Sums of two matrices.
         * \{
         *
         * Returns the first matrix after adding each element of the second one.
         *
         * \param a The DenseMatrix to be used.
         * \param b The DenseMatrix to be added.
         *
         * \retval a The referenced matrix is changed by adding each corresponding element to each of its elements.
         */

        //static DenseMatrix<float> & value(DenseMatrix<float> & a, const DenseMatrix<float> & b);

        //static DenseMatrix<double> & value(DenseMatrix<double> & a, const DenseMatrix<double> & b);

        /**
         * \name Sums of two vectors
         * \{
         *
         * \brief Returns the the sum of two given vectors.
         *
         * \param a The vector that is the left-hand summand of the operation.
         * \param b The vector that is the right-hand summand of the operation.
         *
         * \retval Will modify the summand a and return it.
         *
         * \exception VectorSizeDoesNotMatch is thrown if the two vectors don't have the same size.
         */

        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b);

        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b);
        /// \}

    };

    /**
     * \brief Sum of two entities
     *
     * Sum is the class template for the addition operation
     * \f[
     *     \texttt{Sum}(a, b): \quad r \leftarrow a + b,
     * \f]
     * which yields r, the sum of entities a and b.
     *
     * The return value is the summand a after modification.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <> struct Sum<tags::CPU::SSE>
    {
        /**
         * \name Sums of two vectors
         * \{
         *
         * \brief Returns the the sum of two given vectors.
         *
         * \param a The vector that is the left-hand summand of the operation.
         * \param b The vector that is the right-hand summand of the operation.
         *
         * \retval Will modify the summand a and return it.
         *
         * \exception VectorSizeDoesNotMatch is thrown if the two vectors don't have the same size.
         */

        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b);

        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b);

        static DenseMatrix<float> & value(DenseMatrix<float> & a, const DenseMatrix<float> & b);

        static DenseMatrix<double> & value(DenseMatrix<double> & a, const DenseMatrix<double> & b);

        static DenseMatrix<float> value(DenseMatrix<float> &a, const BandedMatrix<float> & b)
        {
            CONTEXT("When adding DenseMatrix and BandedMatrix (SSE forwarding to CPU):");
            return Sum<tags::CPU>::value(a, b);
        }

        static DenseMatrix<double> value(DenseMatrix<double> &a, const BandedMatrix<double> & b)
        {
            CONTEXT("When adding DenseMatrix and BandedMatrix (SSE forwarding to CPU):");
            return Sum<tags::CPU>::value(a, b);
        }

        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & x, const float a);

        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & x, const double a);

        static DenseMatrix<float> & value(DenseMatrix<float> & x, const float a);

        static DenseMatrix<double> & value(DenseMatrix<double> & x, const double a);
        /// \}
    };

    template <> struct Sum<tags::OpenCL::CPU>
    {
        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b);

        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b);
    };

    template <> struct Sum<tags::OpenCL::GPU>
    {
        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b);

        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b);
    };

    /**
     * \brief Sum of two entities
     *
     * Sum is the class template for the addition operation
     * \f[
     *     \texttt{Sum}(a, b): \quad r \leftarrow a + b,
     * \f]
     * which yields r, the sum of entities a and b.
     *
     * The return value is the summand a after modification.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */

    namespace mc
    {
        template <typename Tag_> struct Sum
        {
            template <typename DT1_, typename DT2_>
            static DenseVectorBase<DT1_> & value(DenseVectorBase<DT1_> & x, const DenseVectorBase<DT2_> & y)
            {
                CONTEXT("When calculating Sum (DenseVectorBase, DenseVectorBase) using backend : " + Tag_::name);

                if (x.size() != y.size())
                    throw VectorSizeDoesNotMatch(y.size(), x.size());

                unsigned long min_part_size(Configuration::instance()->get_value("mc::Sum(DVB,DVB)::min_part_size", 128));
                unsigned long max_count(Configuration::instance()->get_value("mc::Sum(DVB,DVB)::max_count",
                            mc::ThreadPool::instance()->num_threads()));

                Operation<honei::Sum<typename Tag_::DelegateTo> >::op(x, y, min_part_size, max_count);

                return x;
            }

            template <typename DT1_, typename DT2_>
            static DenseVectorContinuousBase<DT1_> & value(DenseVectorContinuousBase<DT1_> & x, const DenseVectorContinuousBase<DT2_> & y)
            {
                CONTEXT("When calculating Sum (DenseVectorContinuousBase, DenseVectorContinuousBase) using backend : " + Tag_::name);

                if (x.size() != y.size())
                    throw VectorSizeDoesNotMatch(y.size(), x.size());

                unsigned long min_part_size(Configuration::instance()->get_value("mc::Sum(DVCB,DVCB)::min_part_size", 128));
                unsigned long max_count(Configuration::instance()->get_value("mc::Sum(DVCB,DVCB)::max_count",
                            mc::ThreadPool::instance()->num_threads()));

                Operation<honei::Sum<typename Tag_::DelegateTo> >::op(x, y, min_part_size, max_count);

                return x;
            }
/*
            template <typename DT1_, typename DT2_>
            static DenseVectorBase<DT1_> & value(DenseVectorBase<DT1_> & x, const DT2_ & a)
            {
                CONTEXT("When calculating Sum (DenseVectorBase, DT) using backend : " + Tag_::name);

                unsigned long min_part_size(Configuration::instance()->get_value("mc::Sum(DVB,DT)::min_part_size", 128));
                unsigned long max_count(Configuration::instance()->get_value("mc::Sum(DVB,DT)::max_count",
                            mc::ThreadPool::instance()->num_threads()));
                Operation<honei::Sum<typename Tag_::DelegateTo> >::op(x, a, min_part_size, max_count);

                return x;
            }

            template <typename DT1_, typename DT2_>
            static DenseVectorContinuousBase<DT1_> & value(DenseVectorContinuousBase<DT1_> & x, const DT2_ & a)
            {
                CONTEXT("When calculating Sum (DenseVectorContinuousBase, DT) using backend : " + Tag_::name);

                unsigned long min_part_size(Configuration::instance()->get_value("mc::Sum(DVCB,DT)::min_part_size", 128));
                unsigned long max_count(Configuration::instance()->get_value("mc::Sum(DVCB,DT)::max_count",
                            mc::ThreadPool::instance()->num_threads()));

                Operation<honei::Sum<typename Tag_::DelegateTo> >::op(x, a, min_part_size, max_count);

                return x;
            }
*/

            // Dummy
            template <typename DT1_, typename DT2_>
            static DenseVectorBase<DT1_> & value(DenseVectorBase<DT1_> & a, const SparseVector<DT2_> & b)
            {
                CONTEXT("When calculating Sum (DenseVectorBase, SparseVector) using backend : " + Tag_::name);

                if (a.size() != b.size())
                    throw VectorSizeDoesNotMatch(b.size(), a.size());

                return honei::Sum<tags::CPU>::value(a, b);
             }

            // Dummy
            template <typename DT1_, typename DT2_>
            static DenseVectorContinuousBase<DT1_> & value(DenseVectorContinuousBase<DT1_> & a, const SparseVector<DT2_> & b)
            {
                CONTEXT("When calculating Sum (DenseVectorContinuousBase, SparseVector) using backend : " + Tag_::name);

                if (a.size() != b.size())
                    throw VectorSizeDoesNotMatch(b.size(), a.size());

                return honei::Sum<tags::CPU>::value(a, b);
             }

            // Dummy
            template <typename DT1_, typename DT2_>
            static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & x, const DT2_ & a)
            {
                CONTEXT("When calculating Sum (DenseMatrix, DT) using backend : " + Tag_::name);

                return honei::Sum<tags::CPU>::value(x, a);
            }

            // Dummy
            template <typename DT1_, typename DT2_>
            static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
            {
                CONTEXT("When calculating Sum (SparseMatrix, SparseMatrix) using backend : " + Tag_::name);

                if (a.columns() != b.columns())
                {
                    throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
                }

                if (a.rows() != b.rows())
                {
                    throw MatrixRowsDoNotMatch(b.rows(), a.rows());
                }

                return honei::Sum<tags::CPU>::value(a, b);
             }

            // Dummy
            template <typename DT1_, typename DT2_>
            static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
            {
                CONTEXT("When calculating Sum (DenseMatrix, DenseMatrix) using backend : " + Tag_::name);

                if (a.columns() != b.columns())
                {
                    throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
                }

                if (a.rows() != b.rows())
                {
                    throw MatrixRowsDoNotMatch(b.rows(), a.rows());
                }

                return honei::Sum<tags::CPU>::value(a, b);
            }

            // Dummy
            template <typename DT1_, typename DT2_>
            static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
            {
                CONTEXT("When calculating Sum (DenseMatrix, SparseMatrix) using backend : " + Tag_::name);

                if (a.columns() != b.columns())
                {
                    throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
                }

                if (a.rows() != b.rows())
                {
                    throw MatrixRowsDoNotMatch(b.rows(), a.rows());
                }

                return honei::Sum<tags::CPU>::value(a, b);
            }

            // Dummy
            template <typename DT1_, typename DT2_>
            static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
            {
                CONTEXT("When calculating Sum (SparseMatrix, BandedMatrix) using backend : " + Tag_::name);

                return honei::Sum<tags::CPU>::value(a, b);
            }

            // Dummy
            template <typename DT1_, typename DT2_>
            static BandedMatrix<DT1_> & value(BandedMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
            {
                CONTEXT("When calculating Sum (BandedMatrix, BandedMatrix) using backend : " + Tag_::name);

                if (a.rows() != b.rows())
                {
                    throw MatrixSizeDoesNotMatch(b.rows(), a.rows());
                }

                return honei::Sum<tags::CPU>::value(a, b);
            }

            // Dummy
            template <typename DT1_, typename DT2_>
            static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
            {
                CONTEXT("When calculating Sum (DenseMatrix, BandedMatrix) using backend : " + Tag_::name);

                if (a.columns() != a.rows())
                {
                    throw MatrixIsNotSquare(a.rows(), a.columns());
                }

                if (a.rows() != b.rows())
                {
                    throw MatrixRowsDoNotMatch(b.rows(), a.rows());
                }

                return honei::Sum<tags::CPU>::value(a, b);
            }
        };
    }

    template <> struct Sum<tags::CPU::MultiCore> :
        public mc::Sum<tags::CPU::MultiCore>
    {
    };

    template <> struct Sum<tags::CPU::MultiCore::SSE> :
        public mc::Sum<tags::CPU::MultiCore::SSE>
    {
    };
}

#endif
