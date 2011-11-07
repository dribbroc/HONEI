/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Sven Mallach <sven.mallach@honei.org>
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
#ifndef LIBLA_GUARD_DOT_PRODUCT_HH
#define LIBLA_GUARD_DOT_PRODUCT_HH 1

#include <honei/la/dense_vector.hh>
#include <honei/la/dense_vector_range.hh>
#include <honei/la/sparse_vector.hh>
#include <honei/la/vector_error.hh>
#include <honei/util/benchmark_info.hh>
#include <honei/util/configuration.hh>
#include <honei/util/operation_wrapper.hh>
#include <honei/util/partitioner.hh>
#include <honei/backends/multicore/operation.hh>
#include <honei/backends/multicore/thread_pool.hh>
#include <honei/util/tags.hh>
#include <honei/mpi/operations.hh>
#include <honei/mpi/dense_vector_mpi-fwd.hh>

namespace honei
{
    template <typename Tag_ = tags::CPU> struct DotProduct;

    /**
     * \brief DotProduct of two vectors.
     *
     * DotProduct is the class template for the operation
     * \f[
     *     \texttt{DotProduct}(x, y): \quad r \leftarrow x \cdot y,
     * \f]
     * which yields the dor or inner product of the given vectors x and y.
     *
     * \ingroup grplaoperations
     * \ingroup grplavectoroperations
     */
    template <> struct DotProduct<tags::CPU>
    {
        /**
         * \name Dot products
         * \{
         *
         * Returns the dot product of two given vectors.
         *
         * \param x One of the vectors of which the scalar product shall be computed.
         * \param y idem
         *
         * \retval r Will return an instance of the used data type containing the scalar product.
         *
         * \exception VectorSizeDoesNotMatch is thrown if the two vectors don't have the same size.
         */

        template <typename DT1_, typename DT2_>
        static DT1_ value(const DenseVectorBase<DT1_> & x, const DenseVectorBase<DT2_> & y)
        {
            CONTEXT("When calculating DenseVectorBase-DenseVectorBase dot product:");

            if (x.size() != y.size())
                throw VectorSizeDoesNotMatch(y.size(), x.size());

            DT1_ result(0);

            for (typename DenseVector<DT1_>::ConstElementIterator l(x.begin_elements()),
                    l_end(x.end_elements()) ; l != l_end ; ++l )
            {
                result += (*l) * y[l.index()];
            }

            return result;
        }

        template <typename DT1_, typename DT2_>
        static DT1_ value(const SparseVector<DT1_> & x, const DenseVectorBase<DT2_> & y)
        {
            CONTEXT("When calculating SparseVector-Vector dot product:");

            if (x.size() != y.size())
                throw VectorSizeDoesNotMatch(y.size(), x.size());

            DT1_ result(0);

            for (typename SparseVector<DT1_>::NonZeroConstElementIterator l(x.begin_non_zero_elements()),
                    l_end(x.end_non_zero_elements()) ; l != l_end ; ++l )
            {
                result += (*l) * y[l.index()];
            }

            return result;
        }

        template <typename DT1_, typename DT2_>
        static DT1_ value(const DenseVectorBase<DT1_> & x, const SparseVector<DT2_> & y)
        {
            CONTEXT("When calculating Vector-SparseVector dot product:");

            if (x.size() != y.size())
                throw VectorSizeDoesNotMatch(y.size(), x.size());

            DT1_ result(0);

            for (typename SparseVector<DT2_>::NonZeroConstElementIterator l(y.begin_non_zero_elements()),
                    l_end(y.end_non_zero_elements()) ; l != l_end ; ++l )
            {
                result += (*l) * x[l.index()];
            }

            return result;
        }


        template <typename DT1_, typename DT2_>
        static DT1_ value(const SparseVector<DT1_> & x, const SparseVector<DT2_> & y)
        {
            CONTEXT("When calculating SparseVector-SparseVector dot product:");

            if (x.size() != y.size())
                throw VectorSizeDoesNotMatch(y.size(), x.size());

            DT1_ result(0);

            /*typename SparseVector<DT2_>::NonZeroConstElementIterator r(y.begin_non_zero_elements());
            for (typename SparseVector<DT1_>::NonZeroConstElementIterator l(x.begin_non_zero_elements()),
                    l_end(x.end_non_zero_elements()) ; l != l_end ; )
            {
                if (l.index() == r.index())
                {
                    result += (*l) * (*r);
                    ++l; ++r;
                }
                else if (l.index() < r.index())
                {
                    ++l;
                }
                else
                {
                    ++r;
                }
            }*/
            unsigned long xi(0), yi(0);
            const unsigned long * xidx(x.indices());
            const unsigned long * yidx(y.indices());
            const DT1_ * xval(x.elements());
            const DT2_ * yval(y.elements());
            while (xidx[xi] < x.size() && yidx[yi] < y.size())
            {
                if(xidx[xi] == yidx[yi])
                {
                    result += xval[xi] * yval[yi];
                    ++xi, ++yi;
                }
                else if (xidx[xi] < yidx[yi])
                    ++xi;
                else
                    ++yi;
            }

            return result;
        }

        template <typename IT1_, typename IT2_>
        static typename IT1_::value_type value(IT1_ & x, const IT1_ & x_end,
                IT2_ & y, const IT2_ & y_end)
        {
            CONTEXT("When calculating iterator-based dot product:");

            typename IT1_::value_type result(0);

            for ( ; (x != x_end) && (y != y_end) ; ++x, ++y)
            {
                result += (*x) * (*y);
            }

            return result;
        }

        template <typename DT_>
        static inline DT_ value(const DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y)
        {
            return MPIOps<tags::CPU>::dot_product(x, y);
        }

        /// \}

        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(DenseVectorBase<DT1_> & a, DenseVectorBase<DT2_> & b)
        {
            BenchmarkInfo result;
            result.flops = 2 * a.size();
            result.load = a.size() * (sizeof(DT1_) + sizeof(DT2_));
            result.store = sizeof(DT1_);
            result.size.push_back(a.size());
            result.size.push_back(b.size());
            return result; 
        }

        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(SparseVector<DT1_> & a, DenseVectorBase<DT2_> & b)
        {
            BenchmarkInfo result;
            result.flops = 2 * a.used_elements();
            result.load = a.used_elements() * (sizeof(DT1_) + sizeof(DT2_));
            result.store = sizeof(DT1_);
            result.size.push_back(a.size());
            result.size.push_back(b.size());
            result.scale = (double(a.size()) / a.used_elements());
            return result;
        }

        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(SparseVector<DT1_> & a, SparseVector<DT2_> & b)
        {
            BenchmarkInfo result;
            result.store = sizeof(DT1_);
            typename SparseVector<DT2_>::NonZeroConstElementIterator r(b.begin_non_zero_elements());
            for (typename SparseVector<DT1_>::NonZeroConstElementIterator l(a.begin_non_zero_elements()), l_end(a.end_non_zero_elements()) ; l != l_end ; )
            {
                if (l.index() == r.index())
                {
                    result.flops += 2;
                    result.load += (sizeof(DT1_) + sizeof(DT2_));
                    ++l; ++r;
                }
                else if (l.index() < r.index())
                {
                    ++l;
                }
                else
                {
                    ++r;
                }
            }
            result.size.push_back(a.size());
            result.size.push_back(b.size());
            result.scale = (double(2 * a.size) / result.flops);
            return result;
        }
    };

    /**
     * \brief DotProduct of two vectors.
     *
     * DotProduct is the class template for the operation
     * \f[
     *     \texttt{DotProduct}(x, y): \quad r \leftarrow x \cdot y,
     * \f]
     * which yields the dot or inner product of the given vectors x and y.
     *
     * \ingroup grplaoperations
     * \ingroup grplavectoroperations
     */
    template <> struct DotProduct<tags::Cell>
    {
        /**
         * \name Dot products
         * \{
         *
         * Returns the dot product of two given vectors.
         *
         * \param x One of the vectors of which the dot product shall be computed.
         * \param y idem
         *
         * \retval r Will return an instance of the used data type containing the scalar product.
         *
         * \exception VectorSizeDoesNotMatch is thrown if the two vectors don't have the same size.
         */

        static float value(const DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b);
        static double value(const DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b);

        /// \}
    };

    /**
     * \brief DotProduct of two vectors.
     *
     * DotProduct is the class template for the operation
     * \f[
     *     \texttt{DotProduct}(x, y): \quad r \leftarrow x \cdot y,
     * \f]
     * which yields the dot or inner product of the given vectors x and y.
     *
     * \ingroup grplaoperations
     * \ingroup grplavectoroperations
     */
    template <> struct DotProduct<tags::GPU::CUDA>
    {
        /**
         * \{
         *
         * Returns the dot product of two given vectors.
         *
         * \param x One of the vectors of which the dot product shall be computed.
         * \param y idem
         *
         * \retval r Will return an instance of the used data type containing the scalar product.
         *
         * \exception VectorSizeDoesNotMatch is thrown if the two vectors don't have the same size.
         */

        static float value(const DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b);

        static double value(const DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b);

        /// \}
    };

    template <> struct DotProduct<tags::GPU::MultiCore::CUDA>
    {
        /**
         * \{
         *
         * Returns the dot product of two given vectors.
         *
         * \param x One of the vectors of which the dot product shall be computed.
         * \param y idem
         *
         * \retval r Will return an instance of the used data type containing the scalar product.
         *
         * \exception VectorSizeDoesNotMatch is thrown if the two vectors don't have the same size.
         */

        static float value(const DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b);

        static double value(const DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b);

        /// \}
    };

    /**
     * \brief DotProduct of two vectors.
     *
     * DotProduct is the class template for the operation
     * \f[
     *     \texttt{DotProduct}(x, y): \quad r \leftarrow x \cdot y,
     * \f]
     * which yields the dot or inner product of the given vectors x and y.
     *
     * \ingroup grplaoperations
     * \ingroup grplavectoroperations
     */
    template <> struct DotProduct<tags::CPU::SSE>
    {
        /**
         * \{
         *
         * Returns the dot product of two given vectors.
         *
         * \param x One of the vectors of which the dot product shall be computed.
         * \param y idem
         *
         * \retval r Will return an instance of the used data type containing the scalar product.
         *
         * \exception VectorSizeDoesNotMatch is thrown if the two vectors don't have the same size.
         */

        static float value(const DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b);

        static double value(const DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b);

        template <typename DT_>
        static inline DT_ value(const DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y)
        {
            return MPIOps<tags::CPU::SSE>::dot_product(x, y);
        }

        /// \}
    };

    template <> struct DotProduct<tags::OpenCL::CPU>
    {
        template <typename DT_>
        static DT_ value(const DenseVectorContinuousBase<DT_> & a, const DenseVectorContinuousBase<DT_> & b);
    };

    template <> struct DotProduct<tags::OpenCL::GPU>
    {
        template <typename DT_>
        static DT_ value(const DenseVectorContinuousBase<DT_> & a, const DenseVectorContinuousBase<DT_> & b);
    };

    namespace mc
    {
        template <typename Tag_> struct DotProduct
        {
            template <typename DT1_, typename DT2_>
            static DT1_ value(const DenseVectorContinuousBase<DT1_> & x, const DenseVectorContinuousBase<DT2_> & y)
            {
                CONTEXT("When calculating DVCB, DBCB dot product using backend : " + Tag_::name);

                if (x.size() != y.size())
                    throw VectorSizeDoesNotMatch(y.size(), x.size());

                unsigned long min_part_size(Configuration::instance()->get_value("mc::dot_product(DVCB,DVCB)::min_part_size", 128));
                unsigned long max_count(Configuration::instance()->get_value("mc::dot_product(DVCB,DVCB)::max_count",
                            mc::ThreadPool::instance()->num_threads()));

                DT1_ result(0);

                Operation<honei::DotProduct<typename Tag_::DelegateTo> >::op(result, x, y, min_part_size, max_count);

                return result;
            }
        };
    }

    /**
     * \brief DotProduct of two vectors.
     *
     * DotProduct is the class template for the operation
     * \f[
     *     \texttt{DotProduct}(x, y): \quad r \leftarrow x \cdot y,
     * \f]
     * which yields the dot or inner product of the given vectors x and y.
     *
     * \ingroup grplaoperations
     * \ingroup grplavectoroperations
     */

    template <> struct DotProduct <tags::CPU::MultiCore> :
        public mc::DotProduct<tags::CPU::MultiCore>
    {
    };

    template <> struct DotProduct<tags::CPU::MultiCore::SSE> :
        public mc::DotProduct<tags::CPU::MultiCore::SSE>
    {
    };
}

#endif
