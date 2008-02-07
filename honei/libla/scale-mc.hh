/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Volker Jung <volker.jung@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_SCALE_MC_HH
#define LIBLA_GUARD_SCALE_MC_HH 1

#include <honei/libla/dense_matrix.hh>
#include <honei/libla/banded_matrix.hh>
#include <honei/libla/scale.hh>
#include <honei/libla/sparse_matrix.hh>
#include <honei/libutil/configuration.hh>
#include <honei/libutil/log.hh>
#include <honei/libutil/partitioner.hh>
#include <honei/libutil/stringify.hh>
#include <honei/libutil/thread_pool.hh>
#include <honei/libutil/wrapper.hh>
#include <honei/libutil/tags.hh>
#include <list>


namespace honei
{
    template <typename Tag_> struct Scale;

    template <typename Tag_> struct MCScale
    {
        template <typename DT_>
        static inline DT_ calculate_parts(const DT_ ref)
        {
            if (ref / 4096 < 8) return ref / 4096;
            return 8;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT2_> & value(DenseMatrix<DT2_> & x, const DT1_ a)
        {
            CONTEXT("When scaling DenseMatrix (MultiCore):");

            ThreadPool * tp(ThreadPool::get_instance());

            std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

            unsigned long min_part_size(Configuration::instance()->get_value("mc::scale[DM]::min-part-size", 1024));
            unsigned long overall_size(x.columns());


            if (overall_size >= (min_part_size << 1))
            {

                PartitionList partitions;

                Partitioner<tags::CPU::MultiCore>::Partitioner(
                    Configuration::instance()->get_value("mc::scale[DM]::max-count",
                        2 * Configuration::instance()->get_value("mc::num-cores", 2)),
                    min_part_size, 16, overall_size, PartitionList::Filler(partitions));

                for (unsigned long i(0); i < x.rows(); ++i)
                {
                    for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()) ;
                            p != p_end ; ++p)
                    {
                        DenseVectorRange<DT2_> range(x[i].range(p->size, p->start));
                        TwoArgWrapper< Scale<typename Tag_::DelegateTo>, DenseVectorRange<DT2_>, const DT1_ > wrapper(range, a);
                        std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(wrapper));
                        dispatched_tasks.push_back(ptr);
                    }
                }
            }
            else
            {
                for (unsigned long i(0); i < x.rows(); ++i)
                {
                    TwoArgWrapper< Scale<typename Tag_::DelegateTo>, DenseVectorRange<DT2_>, const DT1_ > wrapper(x[i], a);
                    std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(wrapper));
                    dispatched_tasks.push_back(ptr);
                }
            }

            while(! dispatched_tasks.empty())
            {
                dispatched_tasks.front()->wait_on();
                dispatched_tasks.pop_front();
            }
            return x;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT2_> & value(SparseMatrix<DT2_> & x, const DT1_ a)
        {
            CONTEXT("When scaling SparseMatrix (MultiCore):");

            ThreadPool * tp(ThreadPool::get_instance());

            std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

            unsigned long min_part_size(Configuration::instance()->get_value("mc::scale[SM]::min-part-size",
                                            Configuration::instance()->get_value("mc::default-min-size", 1024)));

            for (unsigned long i(0); i < x.rows(); ++i)
            {
                if (x[i].used_elements() >= (min_part_size << 1))
                {
                    PartitionList partitions;

                    Partitioner<tags::CPU::MultiCore>(
                            Configuration::instance()->get_value("mc::scale[SM]::max-count",
                                2 * Configuration::instance()->get_value("mc::num-cores", 2)),
                            min_part_size, 16, x[i].used_elements(), PartitionList::Filler(partitions));

                    for ( PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()) ; p != p_end ; ++p)
                    {
                        typename Vector<DT2_>::ElementIterator start(x[i].non_zero_element_at(p->start)),
                                stop(x[i].non_zero_element_at(p->start + p->size));
                        start += p->start;
                        stop += p->start + p->size;
                        ThreeArgWrapper< MCScale<Tag_>, typename Vector<DT2_>::ElementIterator, typename Vector<DT2_>::ElementIterator, const DT1_ >
                            wrapper(start, stop, a);
                        std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(wrapper));
                        dispatched_tasks.push_back(ptr);
                    }
                }
                else
                {
                    typename Vector<DT2_>::ElementIterator start(x[i].begin_non_zero_elements()), stop(x[i].end_non_zero_elements());
                        ThreeArgWrapper< MCScale<Tag_>, typename Vector<DT2_>::ElementIterator, typename Vector<DT2_>::ElementIterator, const DT1_ >
                            wrapper(start, stop, a);
                    std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(wrapper));
                    dispatched_tasks.push_back(ptr);
                }
            }

            while (! dispatched_tasks.empty())
            {
                dispatched_tasks.front()->wait_on();
                dispatched_tasks.pop_front();
            }

            return x;
        }

        template <typename DT1_, typename DT2_>
        static BandedMatrix<DT2_> & value(BandedMatrix<DT2_> & x, const DT1_ a)
        {
            CONTEXT("When scaling BandedMatrix (MultiCore):");

            ThreadPool * tp(ThreadPool::get_instance());

            std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

            unsigned long min_part_size(Configuration::instance()->get_value("mc::scale[BM]::min-part-size",
                                                Configuration::instance()->get_value("mc::default-min-size", 1024)));

            unsigned long max_count(Configuration::instance()->get_value("mc::scale[BM]::max-count",
                                        2 * Configuration::instance()->get_value("mc::num-cores", 2)));

            typename BandedMatrix<DT2_>::VectorIterator vi(x.begin_bands());

            // Calculating lower triangular matrix.
            for (typename BandedMatrix<DT2_>::VectorIterator vi_end(x.band_at(x.size() - 1)) ;
                    vi != vi_end ; ++vi)
            {
                unsigned long band_size(x.size() - vi.index() - 1);

                if (band_size >= (min_part_size << 1))
                {
                    PartitionList partitions;

                    Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, band_size, PartitionList::Filler(partitions));

                    for ( typename PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()); p != p_end ; ++p)
                    {
                        DenseVectorRange<DT2_> range((*vi).range(p->size, p->start + x.size() - band_size));
                        TwoArgWrapper< Scale<typename Tag_::DelegateTo>, DenseVectorRange<DT2_>, const DT1_ >
                            wrapper(range, a);
                        std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(wrapper));
                        dispatched_tasks.push_back(ptr);
                    }
                }
                else
                {
                    DenseVectorRange<DT2_> range((*vi).range(band_size, x.size() - band_size));
                    TwoArgWrapper< Scale<typename Tag_::DelegateTo>, DenseVectorRange<DT2_>, const DT1_ >  wrapper(range, a);
                    std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(wrapper));
                    dispatched_tasks.push_back(ptr);
                }
            }

            // Calculating diagonal band.
            {
                typename BandedMatrix<DT2_>::VectorIterator vi_diag(x.band_at(x.size() - 1));
                if (vi == vi_diag)
                {

                    if (x.size() >= (min_part_size << 1))
                    {
                        PartitionList partitions;

                        Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, x.size(), PartitionList::Filler(partitions));

                        for ( typename PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()) ; p != p_end; ++p)
                        {
                            DenseVectorRange<DT2_> range((*vi).range(p->size, p->start));
                            TwoArgWrapper< Scale<typename Tag_::DelegateTo>, DenseVectorRange<DT2_>, const DT1_ >
                                wrapper(range, a);
                            std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(wrapper));
                            dispatched_tasks.push_back(ptr);
                        }
                    }
                    else
                    {
                        TwoArgWrapper< Scale<typename Tag_::DelegateTo>, DenseVector<DT2_>, const DT1_ > wrapper((*vi), a);
                        std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(wrapper));
                        dispatched_tasks.push_back(ptr);
                    }
                }
            }

            ++vi;

            // Calculating upper traingular matrix.
            for (typename BandedMatrix<DT2_>::VectorIterator vi_end(x.end_bands()) ; vi != vi_end ; ++vi)
            {
                unsigned long band_size(2 * x.size() - vi.index() - 1);

                if (band_size >= (min_part_size << 1))
                {
                    PartitionList partitions;

                    Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, band_size, PartitionList::Filler(partitions));

                    for ( typename PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()) ; p != p_end ; ++p)
                    {
                        DenseVectorRange<DT2_> range((*vi).range(p->size, p->start));
                        TwoArgWrapper< Scale<typename Tag_::DelegateTo>, DenseVectorRange<DT2_>, const DT1_ >
                            wrapper(range, a);
                        std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(wrapper));
                        dispatched_tasks.push_back(ptr);
                    }
                }
                else
                {
                    DenseVectorRange<DT2_> range((*vi).range(band_size, 0));
                    TwoArgWrapper< Scale<typename Tag_::DelegateTo>, DenseVectorRange<DT2_>, const DT1_ >
                        wrapper(range, a);
                    std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(wrapper));
                    dispatched_tasks.push_back(ptr);
                }
            }
                // Wait until all jobs are done.
            while(! dispatched_tasks.empty())
            {
                dispatched_tasks.front()->wait_on();
                dispatched_tasks.pop_front();
            }

            return x;
        }

        template <typename DT1_, typename DT2_>
        static DenseVectorContinuousBase<DT2_> & value(DenseVectorContinuousBase<DT2_> & x, const DT1_ a)
        {
            CONTEXT("When scaling DenseVectorContinuousBase (MultiCore):");

            ThreadPool * tp(ThreadPool::get_instance());

            std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

            unsigned long min_part_size(Configuration::instance()->get_value("mc::scale[DVCB]::min-part-size",
                                            Configuration::instance()->get_value("mc::default-min-size", 1024)));

            unsigned long max_count(Configuration::instance()->get_value("mc::scale[DVCB]::max-count",
                                        Configuration::instance()->get_value("mc::num-cores", 2)));

            if (x.size() >= (min_part_size << 1))
            {
                PartitionList partitions;

                Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, x.size(), PartitionList::Filler(partitions));


                for ( typename PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()) ; p != p_end ; ++p)
                {
                    DenseVectorRange<DT2_> range(x.range(p->size, p->start));
                    TwoArgWrapper< Scale<typename Tag_::DelegateTo>, DenseVectorRange<DT2_>, const DT1_ >
                        wrapper(range, a);
                    std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(wrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while(! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            else
            {
                Scale<typename Tag_::DelegateTo>::value(x, a);
            }

            return x;
        }

        template <typename DT1_, typename DT2_>
        static DenseVectorSlice<DT2_> & value(DenseVectorSlice<DT2_> & x, const DT1_ a)
        {
            CONTEXT("When scaling DenseVectorSlice (MultiCore):");

            unsigned long min_part_size(Configuration::instance()->get_value("mc::scale[DVS]::min-part-size",
                                            Configuration::instance()->get_value("mc::default-min-size", 1024)));

            if (x.size() >= (min_part_size << 1))
            {
                ThreadPool * tp(ThreadPool::get_instance());

                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                PartitionList partitions;

                Partitioner<tags::CPU::MultiCore>(
                        Configuration::instance()->get_value("mc::scale[DVS]::max-count",
                            2 * Configuration::instance()->get_value("mc::num-cores", 2)),
                        min_part_size, 16, x.size(), PartitionList::Filler(partitions)
                        );

                for ( typename PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()); p != p_end ; ++p)
                {
                    typename Vector<DT2_>::ElementIterator start(x.element_at(p->start)), stop(x.element_at(p->size + p->start));
                    ThreeArgWrapper< MCScale<Tag_>, typename Vector<DT2_>::ElementIterator,
                        typename Vector<DT2_>::ElementIterator, const DT1_ > wrapper(start, stop, a);
                    std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(wrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            else
            {
                Scale<typename Tag_::DelegateTo>::value(x, a);
            }

            return x;
        }

        template <typename DT1_, typename DT2_>
        static SparseVector<DT2_> value(SparseVector<DT2_> & x, const DT1_ a)
        {
            CONTEXT("When scaling SparseVector (MultiCore):");

            unsigned long min_part_size(Configuration::instance()->get_value("mc::scale[SV]::min-part-size",
                                            Configuration::instance()->get_value("mc::default-min-size", 1024)));

            if (x.used_elements() >= (min_part_size << 1))
            {
                ThreadPool * tp(ThreadPool::get_instance());

                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                PartitionList partitions;

                Partitioner<tags::CPU::MultiCore>(
                    Configuration::instance()->get_value("mc::scale[SV]::max-count",
                        2 * Configuration::instance()->get_value("mc::num-cores", 2)),
                    min_part_size, 16, x.used_elements(), PartitionList::Filler(partitions));

                for ( typename PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()) ; p != p_end ; ++p)
                {
                    typename Vector<DT2_>::ElementIterator start(x.non_zero_element_at(p->start)),
                            stop(x.non_zero_element_at(p->start + p->size));
                    ThreeArgWrapper< MCScale<Tag_>, typename Vector<DT2_>::ElementIterator,
                        typename Vector<DT2_>::ElementIterator, const DT1_ > wrapper(start, stop, a);
                    std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(wrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            else
            {
                typename Vector<DT2_>::ElementIterator start(x.begin_non_zero_elements()), stop(x.end_non_zero_elements());
                MCScale<Tag_>::value(start, stop, a);
            }

            return x;
        }

        template <typename DT1_, typename IT1_, typename IT2_>
        static void value(IT1_ & x, const IT2_ & x_end, const DT1_ a)
        {
            CONTEXT("When calculating iterator-based scale:");

            for ( ; x < x_end; ++x)
            {
                *x *= a;
            }
        }
    };
}
#endif
