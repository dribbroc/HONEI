/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_DOT_PRODUCT_HH
#define LIBLA_GUARD_DOT_PRODUCT_HH 1

#include <libla/vector.hh>
#include <libla/dense_vector.hh>
#include <libla/sparse_vector.hh>
#include <libla/vector_error.hh>
#include <libutil/tags.hh>

#include <libla/dense_vector_range.hh>
#include <libutil/lock.hh>
#include <libutil/pool_task.hh>
#include <libutil/thread_pool.hh>
#include <libutil/wrapper.hh>
#include <tr1/functional>


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

        #ifdef BENCHM 
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
        #endif
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

        static float value(const DenseVector<float> & a, const DenseVector<float> & b);

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

    template <typename Tag_> struct MCDotProduct
    { 
        template <typename DT1_, typename DT2_>
        static DT1_ value(SparseVector<DT1_> & a, const DenseVectorRange<DT2_> & b, unsigned long offset)
        {
            DT1_ result(0);
            typename Vector<DT1_>::ConstElementIterator r(a.begin_non_zero_elements());
            r += offset;
            offset = r.index();
            unsigned long limit = r.index() + b.size();
            while (r.index() < limit)
            {
                result += *r * b[r.index()-offset];
                ++r;
            }
            return result;
        }

        template <typename DT1_, typename DT2_>
        static DT1_ value(DenseVector<DT1_> & a, const DenseVector<DT2_> & b)
        {
            CONTEXT("When calculating DenseVector-DenseVector dot product (MultiCore):");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());
            DT1_ result(0);
            unsigned long parts(8);
            unsigned long div = a.size() / parts;
            if (div == 0)
            {
                result = DotProduct<typename Tag_::DelegateTo>::value(a,b);
            }
            else
            {
                Mutex mutex;
                unsigned long modulo = a.size() % parts;
                ThreadPool * p(ThreadPool::get_instance());
                PoolTask * pt[parts];
                for (int i(0); i < modulo; ++i)
                {
                    DenseVectorRange<DT1_> range_1(a, div+1, i*(div+1));
                    DenseVectorRange<DT2_> range_2(b, div+1, i*(div+1));
                    ResultTwoArgWrapper<DotProduct<typename Tag_::DelegateTo>, DT1_, const DenseVectorRange<DT1_>, const DenseVectorRange<DT2_> > mywrapper(result, range_1, range_2);
                    pt[i] = p->dispatch(std::tr1::bind(mywrapper, &mutex));
                }
                for (unsigned long i(modulo); i < parts; ++i)
                {
                    DenseVectorRange<DT1_> range_1(a, div, modulo+(i*div));
                    DenseVectorRange<DT2_> range_2(b, div, modulo+(i*div));
                    ResultTwoArgWrapper<DotProduct<typename Tag_::DelegateTo>, DT1_, const DenseVectorRange<DT1_>, const DenseVectorRange<DT2_> > mywrapper(result, range_1, range_2);
                    pt[i] = p->dispatch(std::tr1::bind(mywrapper, &mutex));
                }
                for (unsigned long i(0); i < parts; ++i)
                {
                    pt[i]->wait_on();
                }
            }
            return result;
        }

        template <typename DT1_, typename DT2_>
        static DT1_ value(SparseVector<DT1_> & a, const DenseVector<DT2_> & b)
        {
            CONTEXT("When calculating SparseVector-DenseVector dot product (MultiCore):");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());
            DT1_ result(0);
            unsigned long parts(8);
            unsigned long modulo = a.used_elements() % parts;
            unsigned long div = a.used_elements() / parts;
            if (div == 0) 
            {
                result = DotProduct<typename Tag_::DelegateTo>::value(a, b);
            }
            else
            {
                Mutex mutex;
                ThreadPool * p(ThreadPool::get_instance());
                PoolTask * pt[parts];
                typename Vector<DT2_>::ConstElementIterator r(a.begin_non_zero_elements());
                unsigned long offset;
                for (int i(0); i < modulo; ++i) 
                {
                    offset = r.index();
                    r += div;
                    DenseVectorRange<DT2_> range(b, r.index()-offset+1, offset);
                    ResultThreeArgWrapper<MCDotProduct<Tag_>, DT1_, SparseVector<DT1_>, const DenseVectorRange<DT2_>, const unsigned long > mywrapper(result, a, range, (i*(div+1)));
                    pt[i] = p->dispatch(std::tr1::bind(mywrapper, &mutex));
                    ++r;
                }
                for (unsigned long i(modulo); i < parts; ++i)
                {
                    offset = r.index();
                    r+= div-1;
                    DenseVectorRange<DT2_> range(b, r.index()-offset+1, offset);
                    ResultThreeArgWrapper<MCDotProduct<Tag_>, DT1_, SparseVector<DT1_>, const DenseVectorRange<DT2_>, const unsigned long > mywrapper(result, a, range, modulo + (i*div));
                    pt[i] = p->dispatch(std::tr1::bind(mywrapper, &mutex));
                    ++r;
                }
                for (unsigned long i = 0; i < parts;  ++i)
                {
                    pt[i]->wait_on();
                }
            }
            return result;
        }
    };
    template <> struct DotProduct <tags::CPU::MultiCore> : MCDotProduct <tags::CPU::MultiCore> {};
    template <> struct DotProduct <tags::CPU::MultiCore::SSE> : MCDotProduct <tags::CPU::MultiCore::SSE> {};
}

#endif
