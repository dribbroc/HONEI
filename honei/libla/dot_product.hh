/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_DOT_PRODUCT_HH
#define LIBLA_GUARD_DOT_PRODUCT_HH 1

#include <honei/libla/vector.hh>
#include <honei/libla/dense_vector.hh>
#include <honei/libla/dot_product-mc.hh>
#include <honei/libla/sparse_vector.hh>
#include <honei/libla/vector_error.hh>
#include <honei/util/tags.hh>

#include <honei/libla/dense_vector_range.hh>
#include <honei/util/lock.hh>
#include <honei/util/pool_task.hh>
#include <honei/util/thread_pool.hh>
#include <honei/util/wrapper.hh>
#include <honei/util/benchmark_info.hh>
#include <tr1/functional>


namespace honei
{
    template <typename Tag_> struct MCDotProduct;

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

            for (typename Vector<DT1_>::ConstElementIterator l(x.begin_elements()),
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

            for (typename Vector<DT1_>::ConstElementIterator l(x.begin_non_zero_elements()),
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

            for (typename Vector<DT2_>::ConstElementIterator l(y.begin_non_zero_elements()),
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

            typename Vector<DT2_>::ConstElementIterator r(y.begin_non_zero_elements());
            for (typename Vector<DT1_>::ConstElementIterator l(x.begin_non_zero_elements()),
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
            typename Vector<DT2_>::ConstElementIterator r(b.begin_non_zero_elements());
            for (typename Vector<DT1_>::ConstElementIterator l(a.begin_non_zero_elements()), l_end(a.end_non_zero_elements()) ; l != l_end ; )
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
    template <> struct DotProduct <tags::CPU::MultiCore> : MCDotProduct <tags::CPU::MultiCore> {};
    template <> struct DotProduct <tags::CPU::MultiCore::SSE> : MCDotProduct <tags::CPU::MultiCore::SSE> {};
}

#endif
