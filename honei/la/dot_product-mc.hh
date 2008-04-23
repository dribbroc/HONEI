/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) Joachim Messer <joachim.messer@uni-dortmund.de>
 * Copyright (c) 2007, 2008 Volker Jung <volker.jung@honei.org>
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

#include <honei/la/dense_vector.hh>
#include <honei/la/dense_vector_range.hh>
#include <honei/la/dot_product.hh>
#include <honei/la/sparse_vector.hh>
#include <honei/la/vector.hh>
#include <honei/la/vector_error.hh>
#include <honei/util/lock.hh>
#include <honei/util/partitioner.hh>
#include <honei/util/pool_task.hh>
#include <honei/util/tags.hh>
#include <honei/util/thread_pool.hh>
#include <honei/util/wrapper.hh>

#include <list>
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

            unsigned long min_part_size(Configuration::instance()->get_value("mc::dotproduct[DVCB,DVCB]::min-part-size",
                                            Configuration::instance()->get_value("mc::default-min-size", 1024)));

            DT1_ result(0);
            if (a.size() < (min_part_size << 1))
            {
                result = DotProduct<typename Tag_::DelegateTo>::value(a,b);
            }
            else
            {
                Mutex mutex;
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;
                PartitionList partitions;

                Partitioner<typename tags::CPU::MultiCore>::Partitioner(Configuration::instance()->get_value("mc::dotproduct[DVCB,DVCB]::max-count",
                                2 * Configuration::instance()->get_value("mc::num-cores", 2)),
                            min_part_size, 16, a.size(), PartitionList::Filler(partitions));

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()) ; p != p_end ; ++p)
                {
                    DenseVectorRange<DT1_> range_1(a.range(p->size, p->start));
                    DenseVectorRange<DT2_> range_2(b.range(p->size, p->start));
                    ResultTwoArgWrapper<DotProduct<typename Tag_::DelegateTo>, DT1_, const DenseVectorRange<DT1_>, const DenseVectorRange<DT2_> > mywrapper(result, range_1, range_2);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(std::tr1::bind(mywrapper, &mutex)));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
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

            unsigned long min_part_size(Configuration::instance()->get_value("mc::dot-product[SV,DVCB]::min-part-size",
                        Configuration::instance()->get_value("mc::num-cores", 1024)));

            if (a.used_elements() < (min_part_size << 1))
            {
                result = DotProduct<typename Tag_::DelegateTo>::value(a, b);
            }
            else
            {
                Mutex mutex;
                std::list< PoolTask * > dispatched_tasks;

                PartitionList partitions;

                Partitioner<tags::CPU::MultiCore>::Partitioner(Configuration::instance()->get_value("mc::dot-product[SV,DVCB]::max-count",
                                2 * Configuration::instance()->get_value("mc::num-cores", 2)),
                        min_part_size, 16, a.used_elements(), PartitionList::Filler(partitions));

                typename Vector<DT2_>::ConstElementIterator r(a.begin_non_zero_elements());
                unsigned long offset(0);
                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()); p != p_end ; ++p)
                {
                    offset = r.index();
                    r += (p->size - 1);
                    DenseVectorRange<DT2_> range(b.range(r.index() - offset + 1, offset));
                    ResultThreeArgWrapper<MCDotProduct<Tag_>, DT1_,const SparseVector<DT1_>, const DenseVectorRange<DT2_>, const unsigned long >
                        mywrapper(result, a, range, p->start);
                    dispatched_tasks.push_back(ThreadPool::instance()->dispatch(std::tr1::bind(mywrapper, &mutex)));
                    ++r;
                }
                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    PoolTask * temp(dispatched_tasks.front());
                    dispatched_tasks.pop_front();
                    delete temp;
                }
            }
            return result;
        }

        template <typename DT1_, typename DT2_>
        static DT1_ value(const DenseVectorContinuousBase<DT1_> & a, const SparseVector<DT2_> & b)
        {
            return value(b, a);
        }
    };
}

#endif
