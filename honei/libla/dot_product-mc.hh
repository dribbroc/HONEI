/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) Joachim Messer <joachim.messer@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_DOT_PRODUCT_MC_HH
#define LIBLA_GUARD_DOT_PRODUCT_MC_HH 1

#include <honei/libla/vector.hh>
#include <honei/libla/dense_vector.hh>
#include <honei/libla/dot_product.hh>
#include <honei/libla/sparse_vector.hh>
#include <honei/libla/vector_error.hh>
#include <honei/libutil/tags.hh>

#include <honei/libla/dense_vector_range.hh>
#include <honei/libutil/lock.hh>
#include <honei/libutil/pool_task.hh>
#include <honei/libutil/thread_pool.hh>
#include <honei/libutil/wrapper.hh>
#include <tr1/functional>


namespace honei
{
    template <typename Tag_> struct DotProduct;

    template <typename Tag_> struct MCDotProduct
    {
        template <typename DT1_, typename DT2_>
        static DT1_ value(const SparseVector<DT1_> & a, const DenseVectorRange<DT2_> & b, unsigned long offset)
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
        static DT1_ value(const DenseVectorContinuousBase<DT1_> & a, const DenseVectorContinuousBase<DT2_> & b)
        {
            CONTEXT("When calculating DenseVectorContinuousBase-DenseVectorContinuousBase dot product (MultiCore):");

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
                    DenseVectorRange<DT1_> range_1(a.range(div+1, i*(div+1)));
                    DenseVectorRange<DT2_> range_2(b.range(div+1, i*(div+1)));
                    ResultTwoArgWrapper<DotProduct<typename Tag_::DelegateTo>, DT1_, const DenseVectorRange<DT1_>, const DenseVectorRange<DT2_> > mywrapper(result, range_1, range_2);
                    pt[i] = p->dispatch(std::tr1::bind(mywrapper, &mutex));
                }
                for (unsigned long i(modulo); i < parts; ++i)
                {
                    DenseVectorRange<DT1_> range_1(a.range(div, modulo+(i*div)));
                    DenseVectorRange<DT2_> range_2(b.range(div, modulo+(i*div)));
                    ResultTwoArgWrapper<DotProduct<typename Tag_::DelegateTo>, DT1_, const DenseVectorRange<DT1_>, const DenseVectorRange<DT2_> > mywrapper(result, range_1, range_2);
                    pt[i] = p->dispatch(std::tr1::bind(mywrapper, &mutex));
                }
                for (unsigned long i(0); i < parts; ++i)
                {
                    pt[i]->wait_on();
                    delete pt[i];
                }
            }
            return result;
        }

        template <typename DT1_, typename DT2_>
        static DT1_ value(const SparseVector<DT1_> & a, const DenseVectorContinuousBase<DT2_> & b)
        {
            CONTEXT("When calculating SparseVector-DenseVectorContinuousBase dot product (MultiCore):");

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
                    DenseVectorRange<DT2_> range(b.range(r.index()-offset+1, offset));
                    ResultThreeArgWrapper<MCDotProduct<Tag_>, DT1_,const SparseVector<DT1_>, const DenseVectorRange<DT2_>, const unsigned long > mywrapper(result, a, range, (i*(div+1)));
                    pt[i] = p->dispatch(std::tr1::bind(mywrapper, &mutex));
                    ++r;
                }
                for (unsigned long i(modulo); i < parts; ++i)
                {
                    offset = r.index();
                    r+= div-1;
                    DenseVectorRange<DT2_> range(b.range(r.index()-offset+1, offset));
                    ResultThreeArgWrapper<MCDotProduct<Tag_>, DT1_,const SparseVector<DT1_>, const DenseVectorRange<DT2_>, const unsigned long > mywrapper(result, a, range, modulo + (i*div));
                    pt[i] = p->dispatch(std::tr1::bind(mywrapper, &mutex));
                    ++r;
                }
                for (unsigned long i = 0; i < parts;  ++i)
                {
                    pt[i]->wait_on();
                    delete pt[i];
                }
            }
            return result;
        }
    };
}

#endif
