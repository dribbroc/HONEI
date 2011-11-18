/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007, 2009 Sven Mallach <mallach@honei.org>
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
#ifndef LIBLA_GUARD_ELEMENT_PRODUCT_HH
#define LIBLA_GUARD_ELEMENT_PRODUCT_HH 1

#include <honei/backends/multicore/operation.hh>
#include <honei/backends/multicore/thread_pool.hh>
#include <honei/la/banded_matrix.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/dense_vector_range.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/la/sparse_vector.hh>
#include <honei/la/vector_error.hh>
#include <honei/util/benchmark_info.hh>
#include <honei/util/configuration.hh>
#include <honei/util/tags.hh>
#include <honei/mpi/operations.hh>
#include <honei/mpi/dense_vector_mpi-fwd.hh>

namespace honei
{
    template <typename Tag_ = tags::CPU> struct ElementProduct;
    /**
     * \brief Multiplication of the elements of two given entities.
     *
     * ElementProduct is the template for the multiplication of the elements 
     * \f[
     *     \texttt{ElementProduct}(a,b): \quad a \leftarrow a[i] \cdot b[i],
     * \f]
     *
     * of a given entity.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template<> struct ElementProduct<tags::CPU>
    {
        /**
         * \name Element products
         * \{
         *
         * \brief Returns the the result of elementwise multiplication of two entities.
         *
         * \param a Entity that is a factor of the operation.
         * \param b idem
         *
         * \retval a Will modify the factor a and return it.
         *
         * \exception MatrixSizeDoesNotMatch is thrown if two banded matrices do not have the same size.
         * \exception MatrixRowsDoNotMatch is thrown if two matrices do not have the same number of rows.
         * \exception MatrixColumnsDoNotMatch is thrown if two matrices do not have the same number of columns.
         * \exception MatrixIsNotSquare is thrown if a row access matrix's number of rows does not equal its number of columns.
         */

        template <typename DT1_, typename DT2_>
        static DenseVectorBase<DT1_> & value(DenseVectorBase<DT1_> & a, const DenseVectorBase<DT2_> & b)
        {
            CONTEXT("When calculating the product of DenseVectorBases elements:");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());

            typename DenseVectorBase<DT2_>::ConstElementIterator r(b.begin_elements());
            for (typename DenseVectorBase<DT1_>::ElementIterator l(a.begin_elements()),
                    l_end(a.end_elements()) ; l != l_end ; ++l)
            {
                *l *= *r;
                ++r;
            }

            return a;
        }
        template <typename DT1_, typename DT2_>
        static DenseVectorBase<DT1_> & value(DenseVectorBase<DT1_> & result, DenseVectorBase<DT1_> & a, const DenseVectorBase<DT2_> & b)
        {
            CONTEXT("When calculating the product of DenseVectorBases elements:");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());

            if (a.size() != result.size())
                throw VectorSizeDoesNotMatch(result.size(), a.size());

            typename DenseVectorBase<DT1_>::ConstElementIterator a_i(a.begin_elements());
            typename DenseVectorBase<DT2_>::ConstElementIterator b_i(b.begin_elements());
            for (typename DenseVectorBase<DT1_>::ElementIterator res_i(result.begin_elements()),
                    a_i_end(a.end_elements()) ; a_i != a_i_end ; ++a_i)
            {
                *res_i = *a_i * *b_i;
                ++b_i;
                ++res_i;
            }

            return result;
        }

        template <typename DT1_, typename DT2_>
        static DenseVectorContinuousBase<DT1_> & value(DenseVectorContinuousBase<DT1_> & a, const DenseVectorBase<DT2_> & b)
        {
            CONTEXT("When calculating the product of DenseVectorContinuousBases elements:");

            DenseVectorBase<DT1_> & temp(a);
            ElementProduct<tags::CPU>::value(temp, b);

            return a;
        }

        template <typename DT1_, typename DT2_>
        static DenseVectorContinuousBase<DT1_> & value(DenseVectorContinuousBase<DT1_> & y, const DenseVectorContinuousBase<DT1_> & a, const DenseVectorContinuousBase<DT2_> & b)
        {
            CONTEXT("When calculating the product of DenseVectorContinuousBases elements:");

            typename DenseVectorBase<DT1_>::ConstElementIterator a_i(a.begin_elements());
            typename DenseVectorBase<DT2_>::ConstElementIterator b_i(b.begin_elements());
            typename DenseVectorBase<DT2_>::ElementIterator res_i(y.begin_elements());
            for (typename DenseVectorBase<DT1_>::ConstElementIterator a_i_end(a.end_elements()) ; a_i != a_i_end ; ++a_i)
            {
                *res_i = *a_i * *b_i;
                ++b_i;
                ++res_i;
            }

            return y;
        }

        template <typename DT1_, typename DT2_>
        static SparseVector<DT1_> value(SparseVector<DT1_> & a, const SparseVector<DT2_> & b)
        {
            CONTEXT("When calculating the product of SparseVector elements:");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());

            typename SparseVector<DT2_>::NonZeroConstElementIterator r(b.begin_non_zero_elements());
            for (typename SparseVector<DT1_>::NonZeroElementIterator l(a.begin_non_zero_elements()),
                    l_end(a.end_non_zero_elements()) ; l != l_end ;)
            {
                if (r.index() < l.index())
                {
                    ++r;
                }
                else if (l.index() < r.index())
                {
                    *l = DT1_(0);
                    ++l;
                }
                else
                {
                    *l *= *r;
                    ++l; ++r;
                }
            }
            return a;
            ///\todo: perhaps sparsify - in case l.index < r.index Write of 0 possible.
        }

        template <typename DT1_, typename DT2_>
        static SparseVector<DT1_> & value(SparseVector<DT1_> & a, const DenseVectorBase<DT2_> & b)
        {
            CONTEXT("When calculating the product of SparseVector and DenseVectorBase elements:");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());

            for (typename SparseVector<DT1_>::NonZeroElementIterator l(a.begin_non_zero_elements()),
                    l_end(a.end_non_zero_elements()) ; l != l_end ; ++l )
            {
                *l *= b[l.index()];
            }
            return a;
            ///\todo: perhaps sparsify - if *b[l.index()] == 0 -> write of zero.
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
        {
            CONTEXT("When calculating the product of DenseMatrix elements:");

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
                    l_end(a.end_elements()) ; l != l_end ; ++l)
            {
                *l *= *r;
                ++r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
        {
            CONTEXT("When calculating the product of SparseMatrix and DenseMatrix elements:");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            for (typename SparseMatrix<DT1_>::NonZeroElementIterator l(a.begin_non_zero_elements()),
                    l_end(a.end_non_zero_elements()); l != l_end ; ++l)
            {
                *l *= b(l.row(), l.column());
            }

            return a; ///\todo: perhaps sparsify, dense_matrix[row][col] may be zero.
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
        {
            CONTEXT("When calculating the product of SparseMatrix elements:");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            typename SparseMatrix<DT2_>::NonZeroConstElementIterator r(b.begin_non_zero_elements());
            for (typename SparseMatrix<DT1_>::NonZeroElementIterator l(a.begin_non_zero_elements()),
                    l_end(a.end_non_zero_elements()) ; l != l_end ; )
            {
                if (l.index() < r.index())
                {
                    *l = DT1_(0);
                    ++l;
                }

                else if (r.index() < l.index())
                {
                    ++r;
                }

                else
                {
                    *l *= *r;
                    ++l; ++r;
                }
            }
            ///\todo: perhaps sparsify - in case l.index < r.index set to zero possible.
            return a;
        }

        template <typename DT1_, typename DT2_>
        static BandedMatrix<DT1_> & value(BandedMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
        {
            CONTEXT("When calculating the product of BandedMatrix elements:");

            if (a.rows() != b.rows())
            {
                throw MatrixSizeDoesNotMatch(b.rows(), a.rows());
            }

            typename BandedMatrix<DT2_>::ConstBandIterator r(b.begin_bands());
            for (typename BandedMatrix<DT1_>::BandIterator l(a.begin_bands()),
                    l_end(a.end_bands()) ; l != l_end ; ++l)
            {
                if (! r.exists()) //won't check l.exists() here, cause it is always created by Iterator.
                {
                    ++r;
                    continue;
                }

                DenseVector<DT1_> band(*l);
                ElementProduct<>::value(band, *r);
                ++r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static BandedMatrix<DT1_> & value(BandedMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
        {
            CONTEXT("When calculating the product of BandedMatrix and a DenseMatrix elements:");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            typename DenseMatrix<DT2_>::ConstElementIterator r(b.begin_elements());
            for (typename BandedMatrix<DT1_>::ElementIterator l(a.begin_elements()),
                    l_end(a.end_elements()) ; l != l_end ; ++l)
            {
                *l *= *r;
                ++r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static BandedMatrix<DT1_> & value(BandedMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
        {
            CONTEXT("When calculating the product of BandedMatrix and a SparseMatrix elements:");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            typename BandedMatrix<DT1_>::ElementIterator l(a.begin_elements());
            for (typename SparseMatrix<DT2_>::NonZeroConstElementIterator r(b.begin_non_zero_elements()),
                    r_end(b.end_non_zero_elements()); r != r_end ; )
            {
                while (l.index() < r.index())
                {
                    *l = DT1_(0);
                    ++l;
                }

                *l *= *r;
                ++l; ++r;
            }

            for (typename BandedMatrix<DT1_>::ElementIterator l_end(a.end_elements()) ;
                l != l_end ; ++l)
            {
                *l = DT1_(0);
            }
            /// \todo Left is complete dense at this point.
            return a;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
        {
            CONTEXT("When calculating the product of DenseMatrix and BandedMatrix elements:");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            typename BandedMatrix<DT2_>::ConstElementIterator r(b.begin_elements());
            for (typename DenseMatrix<DT1_>::ElementIterator l(a.begin_elements()),
                    l_end(a.end_elements()) ; l != l_end ; ++l)
            {
                *l *= *r;
                ++r;
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
        {
            CONTEXT("When calculating the product of SparseMatrix and BandedMatrix elements:");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            typename BandedMatrix<DT2_>::ConstElementIterator r(b.begin_elements());
            for (typename SparseMatrix<DT1_>::NonZeroElementIterator l(a.begin_non_zero_elements()),
                    l_end(a.end_non_zero_elements()) ; l != l_end ; )
            {
                while (r.index() < l.index())
                    ++r;

                *l *= *r;
                ++r; ++l;
            }

            return a;
        }

        template <typename DT_>
        static inline DenseVectorMPI<DT_> & value(DenseVectorMPI<DT_> & r, const DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y)
        {
            MPIOps<tags::CPU>::element_product(r, x, y);
            return r;
        }

        /// \}

        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(DenseMatrix<DT1_> & a, DenseMatrix<DT2_> & b)
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
        static inline BenchmarkInfo get_benchmark_info(SparseMatrix<DT1_> & a, DenseMatrix<DT2_> & b)
        {
            BenchmarkInfo result;
            for (typename SparseMatrix<DT1_>::NonZeroElementIterator l(a.begin_non_zero_elements()),
                    l_end(a.end_non_zero_elements()); l != l_end ; ++l)
            {
                result.flops += 1;
                result.load += sizeof(DT1_) + sizeof(DT2_);
                result.store += sizeof(DT1_);
            }
            result.size.push_back(a.rows() * a.columns());
            result.size.push_back(b.rows() * b.columns());
            result.scale = (double(a.rows() * a.columns()) / result.flops);
            return result;
        }

        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(BandedMatrix<DT1_> & a, DenseMatrix<DT2_> & b)
        {
            BenchmarkInfo result;
            result.flops = a.size() * a.size();
            result.load = a.size() * a.size() * (sizeof(DT1_) + sizeof(DT2_));
            result.store = a.size() * a.size() * sizeof(DT1_);
            result.size.push_back(a.rows() * a.columns());
            result.size.push_back(b.rows() * b.columns());
            return result;
        }

        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(DenseVector<DT1_> & a, const DenseVector<DT2_> & b)
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
     * \brief Multiplication of the elements of two given entities.
     *
     * ElementProduct is the template for the multiplication of the elements 
     * \f[
     *     \texttt{ElementProduct}(a,b): \quad a \leftarrow a[i] \cdot b[i],
     * \f]
     *
     * of a given entity.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <> struct ElementProduct<tags::Cell>
    {
        /**
         * \name Element products
         * \{
         *
         * \brief Returns the the result of elementwise multiplication of two entities.
         *
         * \param a Entity that is a factor of the operation.
         * \param b idem
         *
         * \retval a Will modify the factor a and return it.
         *
         * \exception MatrixSizeDoesNotMatch is thrown if two banded matrices do not have the same size.
         * \exception MatrixRowsDoNotMatch is thrown if two matrices do not have the same number of rows.
         * \exception MatrixColumnsDoNotMatch is thrown if two matrices do not have the same number of columns.
         * \exception MatrixIsNotSquare is thrown if a row access matrix's number of rows does not equal its number of columns.
         */

        static DenseMatrix<float> & value(DenseMatrix<float> & a, const DenseMatrix<float> & b);
        static DenseMatrix<double> & value(DenseMatrix<double> & a, const DenseMatrix<double> & b);
        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b);
        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b);

        /// \}
    };

    /**
     * \brief Multiplication of the elements of two given entities.
     *
     * ElementProduct is the template for the multiplication of the elements 
     * \f[
     *     \texttt{ElementProduct}(a,b): \quad a \leftarrow a[i] \cdot b[i],
     * \f]
     *
     * of a given entity.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <> struct ElementProduct<tags::GPU::CUDA>
    {
        /**
         * \name Element products
         * \{
         *
         * \brief Returns the the result of elementwise multiplication of two entities.
         *
         * \param a Entity that is a factor of the operation.
         * \param b idem
         *
         * \retval a Will modify the factor a and return it.
         *
         * \exception MatrixSizeDoesNotMatch is thrown if two banded matrices do not have the same size.
         * \exception MatrixRowsDoNotMatch is thrown if two matrices do not have the same number of rows.
         * \exception MatrixColumnsDoNotMatch is thrown if two matrices do not have the same number of columns.
         * \exception MatrixIsNotSquare is thrown if a row access matrix's number of rows does not equal its number of columns.
         */

        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b);

        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b);

        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & result, const DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b);

        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & result, const DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b);

        static DenseMatrix<float> & value(DenseMatrix<float> & a, const DenseMatrix<float> & b);
        /// \}
    };

    template <> struct ElementProduct<tags::GPU::MultiCore::CUDA>
    {
        /**
         * \name Element products
         * \{
         *
         * \brief Returns the the result of elementwise multiplication of two entities.
         *
         * \param a Entity that is a factor of the operation.
         * \param b idem
         *
         * \retval a Will modify the factor a and return it.
         *
         * \exception MatrixSizeDoesNotMatch is thrown if two banded matrices do not have the same size.
         * \exception MatrixRowsDoNotMatch is thrown if two matrices do not have the same number of rows.
         * \exception MatrixColumnsDoNotMatch is thrown if two matrices do not have the same number of columns.
         * \exception MatrixIsNotSquare is thrown if a row access matrix's number of rows does not equal its number of columns.
         */

        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b);

        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b);

        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & result, const DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b);

        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & result, const DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b);

        /// \}
    };

    /**
     * \brief Multiplication of the elements of two given entities.
     *
     * ElementProduct is the template for the multiplication of the elements
     * \f[
     *     \texttt{ElementProduct}(a,b): \quad a \leftarrow a[i] \cdot b[i],
     * \f]
     *
     * of a given entity.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <> struct ElementProduct<tags::CPU::SSE>
    {
        /**
         * \name Element products
         * \{
         *
         * \brief Returns the the result of elementwise multiplication of two entities.
         *
         * \param a Entity that is a factor of the operation.
         * \param b idem
         *
         * \retval a Will modify the factor a and return it.
         *
         * \exception MatrixSizeDoesNotMatch is thrown if two banded matrices do not have the same size.
         * \exception MatrixRowsDoNotMatch is thrown if two matrices do not have the same number of rows.
         * \exception MatrixColumnsDoNotMatch is thrown if two matrices do not have the same number of columns.
         * \exception MatrixIsNotSquare is thrown if a row access matrix's number of rows does not equal its number of columns.
         */

        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b);

        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b);

        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b, const DenseVectorContinuousBase<float> & c);

        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b, const DenseVectorContinuousBase<double> & c);

        static DenseMatrix<float> & value(DenseMatrix<float> & a, const DenseMatrix<float> & b);

        static DenseMatrix<double> & value(DenseMatrix<double> & a, const DenseMatrix<double> & b);

        static SparseMatrix<float> & value(SparseMatrix<float> & a, const DenseMatrix<float> & b)
        {
            CONTEXT("When multiplying SparseMatrix with DenseMatrix elementwise (SSE forwarding to CPU):");
            return ElementProduct<tags::CPU>::value(a, b);
        }

        static SparseMatrix<double> & value(SparseMatrix<double> & a, const DenseMatrix<double> & b)
        {
            CONTEXT("When multiplying SparseMatrix with DenseMatrix elementwise (SSE forwarding to CPU):");
            return ElementProduct<tags::CPU>::value(a, b);
        }

        template <typename DT_>
        static inline DenseVectorMPI<DT_> & value(DenseVectorMPI<DT_> & r, const DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y)
        {
            MPIOps<tags::CPU::SSE>::element_product(r, x, y);
            return r;
        }

        /// \}
    };

    template <> struct ElementProduct<tags::OpenCL::CPU>
    {
        template <typename DT_>
        static DenseVectorContinuousBase<DT_> & value(DenseVectorContinuousBase<DT_> & a, const DenseVectorContinuousBase<DT_> & b);

        template <typename DT_>
        static DenseVectorContinuousBase<DT_> & value(DenseVectorContinuousBase<DT_> & a, const DenseVectorContinuousBase<DT_> & b, const DenseVectorContinuousBase<DT_> & c);
    };

    template <> struct ElementProduct<tags::OpenCL::GPU>
    {
        template <typename DT_>
        static DenseVectorContinuousBase<DT_> & value(DenseVectorContinuousBase<DT_> & a, const DenseVectorContinuousBase<DT_> & b);

        template <typename DT_>
        static DenseVectorContinuousBase<DT_> & value(DenseVectorContinuousBase<DT_> & a, const DenseVectorContinuousBase<DT_> & b, const DenseVectorContinuousBase<DT_> & c);
    };

    /**
     * \brief Multiplication of the elements of two given entities.
     *
     * ElementProduct is the template for the multiplication of the elements 
     * \f[
     *     \texttt{ElementProduct}(a,b): \quad a \leftarrow a[i] \cdot b[i],
     * \f]
     *
     * of a given entity.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */

    namespace mc
    {
        template <typename Tag_> struct ElementProduct
        {
            template <typename DT1_, typename DT2_>
            static DenseVectorBase<DT1_> & value(DenseVectorBase<DT1_> & x, const DenseVectorBase<DT2_> & y)
            {
                CONTEXT("When calculating ElementProduct (DenseVectorBase, DenseVectorBase) using backend : " + Tag_::name);
                unsigned long min_part_size(Configuration::instance()->get_value("mc::ElementProduct(DVB,DVB)::min_part_size", 128));
                unsigned long max_count(Configuration::instance()->get_value("mc::ElementProduct(DVB,DVB)::max_count",
                            mc::ThreadPool::instance()->num_threads()));

                Operation<honei::ElementProduct<typename Tag_::DelegateTo> >::op(x, y, min_part_size, max_count);

                return x;
            }

            template <typename DT1_, typename DT2_>
            static DenseVectorContinuousBase<DT1_> & value(DenseVectorContinuousBase<DT1_> & x, const DenseVectorContinuousBase<DT2_> & y)
            {
                CONTEXT("When calculating ElementProduct (DenseVectorContinuousBase, DenseVectorContinuousBase) using backend : " + Tag_::name);

                unsigned long min_part_size(Configuration::instance()->get_value("mc::ElementProduct(DVCB,DVCB)::min_part_size", 128));
                unsigned long max_count(Configuration::instance()->get_value("mc::ElementProduct(DVCB,DVCB)::max_count",
                            mc::ThreadPool::instance()->num_threads()));

                Operation<honei::ElementProduct<typename Tag_::DelegateTo> >::op(x, y, min_part_size, max_count);

                return x;
            }

            template <typename DT1_, typename DT2_>
            static DenseVectorContinuousBase<DT1_> & value(DenseVectorContinuousBase<DT1_> & result, const DenseVectorContinuousBase<DT1_> & x, const DenseVectorContinuousBase<DT2_> & y)
            {
                CONTEXT("When calculating ElementProduct (DenseVectorContinuousBase, DenseVectorContinuousBase) using backend : " + Tag_::name);

                unsigned long min_part_size(Configuration::instance()->get_value("mc::ElementProduct(DVCB,DVCB)::min_part_size", 128));
                unsigned long max_count(Configuration::instance()->get_value("mc::ElementProduct(DVCB,DVCB)::max_count",
                            mc::ThreadPool::instance()->num_threads()));

                Operation<honei::ElementProduct<typename Tag_::DelegateTo> >::op(result, x, y, min_part_size, max_count);

                return result;
            }

            template <typename DT_>
            static inline DenseVectorMPI<DT_> & value(DenseVectorMPI<DT_> & r, const DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y)
            {
                MPIOps<Tag_>::element_product(r, x, y);
                return r;
            }

            // Dummy
            template <typename DT1_, typename DT2_>
            static SparseVector<DT1_> & value(SparseVector<DT1_> & a, const DenseVectorBase<DT2_> & b)
            {
                CONTEXT("When calculating ElementProduct (SparseVector, DenseVectorBase) using backend : " + Tag_::name);

                if (a.size() != b.size())
                    throw VectorSizeDoesNotMatch(b.size(), a.size());

                return honei::ElementProduct<tags::CPU>::value(a, b);
             }

            // Dummy
            template <typename DT1_, typename DT2_>
            static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
            {
                CONTEXT("When calculating ElementProduct (SparseMatrix, SparseMatrix) using backend : " + Tag_::name);

                if (a.columns() != b.columns())
                {
                    throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
                }

                if (a.rows() != b.rows())
                {
                    throw MatrixRowsDoNotMatch(b.rows(), a.rows());
                }

                return honei::ElementProduct<tags::CPU>::value(a, b);
             }

            // Dummy
            template <typename DT1_, typename DT2_>
            static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
            {
                CONTEXT("When calculating ElementProduct (DenseMatrix, DenseMatrix) using backend : " + Tag_::name);

                if (a.columns() != b.columns())
                {
                    throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
                }

                if (a.rows() != b.rows())
                {
                    throw MatrixRowsDoNotMatch(b.rows(), a.rows());
                }

                return honei::ElementProduct<tags::CPU>::value(a, b);
            }

            // Dummy
            template <typename DT1_, typename DT2_>
            static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
            {
                CONTEXT("When calculating ElementProduct (SparseMatrix, DenseMatrix) using backend : " + Tag_::name);

                if (a.columns() != b.columns())
                {
                    throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
                }

                if (a.rows() != b.rows())
                {
                    throw MatrixRowsDoNotMatch(b.rows(), a.rows());
                }

                return honei::ElementProduct<tags::CPU>::value(a, b);
            }

            // Dummy
            template <typename DT1_, typename DT2_>
            static BandedMatrix<DT1_> & value(BandedMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
            {
                CONTEXT("When calculating ElementProduct (BandedMatrix, BandedMatrix) using backend : " + Tag_::name);

                if (a.rows() != b.rows())
                {
                    throw MatrixSizeDoesNotMatch(b.rows(), a.rows());
                }

                return honei::ElementProduct<tags::CPU>::value(a, b);
            }

        };
    }

    template <> struct ElementProduct<tags::CPU::MultiCore> :
        public mc::ElementProduct<tags::CPU::MultiCore>
    {
    };

    template <> struct ElementProduct<tags::CPU::MultiCore::SSE> :
        public mc::ElementProduct<tags::CPU::MultiCore::SSE>
    {
    };
}
#endif
