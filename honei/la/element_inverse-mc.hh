/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Mathias Kadolsky <mathkad@gmx.de>
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

#ifndef LIBLA_GUARD_ELEMENT_INVERSE_MC_HH
#define LIBLA_GUARD_ELEMENT_INVERSE_MC_HH 1

#include <honei/la/dense_vector_range.hh>
#include <honei/la/dense_vector_slice.hh>
#include <honei/la/banded_matrix.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/element_inverse.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/la/sparse_vector.hh>
#include <honei/la/vector_error.hh>
#include <honei/util/tags.hh>

#include <honei/la/dense_vector_range.hh>
#include <honei/util/pool_task.hh>
#include <honei/util/thread_pool.hh>
#include <honei/util/wrapper.hh>
#include <honei/util/configuration.hh>
#include <honei/util/partitioner.hh>

#include <algorithm>
#include <iostream>

#define PARTS 8

namespace honei
{
    template <typename Tag_> struct ElementInverse;

    /**
     * \brief Inversion of the elements of the given entity.
     *
     * ElementInverse is the template for the inversion of the elements
     * \f[
     *     \texttt{ElementInverse}(a): \quad a \leftarrow a[i]^{-1},
     * \f]
     *
     * of a given entity.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <typename Tag_> struct MCElementInverse
    {
        template <typename DT1_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & x)
        {
            CONTEXT("When calculating the inverse value of DenseMatrix elements (MultiCore):");

            unsigned long min_part_size(Configuration::instance()->get_value("mc::element_inverse[DM]::min-part-size", 64));
            unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
            unsigned long max_count(Configuration::instance()->get_value("mc::element_inverse[DM]::max-count", num_threads));
            unsigned long overall_size(x.rows());
            if ((overall_size < max_count) && ((x.columns() / min_part_size) > overall_size) && ((x.columns() / min_part_size) >= 2))
            {
                for (unsigned long i(0) ; i < x.rows() ; ++i)
                {
                    ElementInverse<Tag_>::value(x[i]);
                }
            }
            else if (overall_size < 2)
            {
                ElementInverse<typename Tag_::DelegateTo>::value(x);
            }
            else
            {
                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, 1, 1, overall_size, PartitionList::Filler(partitions));
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()) ; p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = p->size;
                    ThreeArgWrapper< ElementInverse<Tag_>, DenseMatrix<DT1_>, unsigned long, unsigned long> mywrapper(x, offset, part_size);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }
                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            return x;
        }

        template <typename DT_>
        static void value(DenseMatrix<DT_> & a, unsigned long start, unsigned long part_size)
        {
            for (unsigned long i(start) ; i < (start + part_size) ; ++i)
            {
                ElementInverse<typename Tag_::DelegateTo>::value(a[i]);
            }
        }

        template <typename DT1_>
        static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & x)
        {
            CONTEXT("When calculating the inverse value of SparseMatrix elements (MultiCore):");

            unsigned long min_part_size(Configuration::instance()->get_value("mc::element_inverse[SM]::min-part-size", 64));
            unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
            unsigned long max_count(Configuration::instance()->get_value("mc::element_inverse[SM]::max-count", num_threads));
            unsigned long overall_size(x.rows());
            if ((overall_size < max_count) && ((x.columns() / min_part_size) > overall_size) && ((x.columns() / min_part_size) >= 2))
            {
                for (unsigned long i(0) ; i < x.rows() ; ++i)
                {
                    ElementInverse<Tag_>::value(x[i]);
                }
            }
            else if (overall_size < 2)
            {
                ElementInverse<typename Tag_::DelegateTo>::value(x);
            }
            else
            {
                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, 1, 1, overall_size, PartitionList::Filler(partitions));
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()) ; p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = p->size;
                    ThreeArgWrapper< ElementInverse<Tag_>, SparseMatrix<DT1_>, unsigned long, unsigned long> mywrapper(x, offset, part_size);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }
                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            return x;
        }

        template <typename DT_>
        static void value(SparseMatrix<DT_> & a, unsigned long start, unsigned long part_size)
        {
            for (unsigned long i(start) ; i < (start + part_size) ; ++i)
            {
                ElementInverse<typename Tag_::DelegateTo>::value(a[i]);
            }
        }

        template <typename DT1_>
        static BandedMatrix<DT1_> & value(BandedMatrix<DT1_> & x)
        {
            CONTEXT("When calculating the inverse value of BandedMatrix elements (MultiCore):");

            ThreadPool * tp(ThreadPool::instance());

            std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

            typename BandedMatrix<DT1_>::VectorIterator vi(x.begin_bands());

            // Calculating lower triangular matrix.
            for (typename BandedMatrix<DT1_>::VectorIterator vi(x.begin_bands()), vi_end(x.band_at(x.size() - 1)) ;
                    vi != vi_end ; ++vi)
            {
                if (! vi.exists())
                    continue;

                unsigned long band_size(x.size() - vi.index() - 1);
                unsigned long chunk_size(band_size / PARTS);
                unsigned long rest(band_size % PARTS);
                unsigned long i(0);

                for ( ; i < rest; ++i)
                {
                    DenseVectorRange<DT1_> range((*vi).range(chunk_size + 1, i * (chunk_size + 1) + x.size() - band_size));
                    OneArgWrapper< ElementInverse<typename Tag_::DelegateTo>, DenseVectorRange<DT1_> >
                        wrapper(range);
                    std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(wrapper));
                    dispatched_tasks.push_back(ptr);
                }

                if (chunk_size != 0)
                {
                    for ( ; i < PARTS; ++i)
                    {
                        DenseVectorRange<DT1_> range((*vi).range(chunk_size, chunk_size + rest + x.size() - band_size));
                        OneArgWrapper< ElementInverse<typename Tag_::DelegateTo>, DenseVectorRange<DT1_> >
                            wrapper(range);
                        std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(wrapper));
                        dispatched_tasks.push_back(ptr);
                    }
                }
            }

            // Calculating diagonal band.
            if (vi.exists())
            {
                unsigned long chunk_size(x.size() / PARTS);
                unsigned long rest(x.size() % PARTS);
                unsigned long i(0);

                for ( ; i < rest; ++i)
                {
                    DenseVectorRange<DT1_> range((*vi).range(chunk_size + 1, i * (chunk_size + 1)));
                    OneArgWrapper< ElementInverse<typename Tag_::DelegateTo>, DenseVectorRange<DT1_> >
                        wrapper(range);
                    std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(wrapper));
                    dispatched_tasks.push_back(ptr);
                }

                if (chunk_size != 0)
                {
                    for ( ; i < PARTS; ++i)
                    {
                        DenseVectorRange<DT1_> range((*vi).range(chunk_size, i *  chunk_size + rest));
                        OneArgWrapper< ElementInverse<typename Tag_::DelegateTo>, DenseVectorRange<DT1_> >
                            wrapper(range);
                        std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(wrapper));
                        dispatched_tasks.push_back(ptr);
                    }
                }
            }

            ++vi;

            // Calculating upper traingular matrix.
            for (typename BandedMatrix<DT1_>::VectorIterator vi_end(x.end_bands()) ; vi != vi_end ; ++vi)
            {
                if (! vi.exists())
                    continue;

                unsigned long band_size(2 * x.size() - vi.index() - 1);
                unsigned long chunk_size(band_size / PARTS);
                unsigned long rest(band_size % PARTS);
                unsigned long i(0);

                for ( ; i < rest; ++i)
                {
                    DenseVectorRange<DT1_> range((*vi).range(chunk_size + 1, i * (chunk_size + 1)));
                    OneArgWrapper< ElementInverse<typename Tag_::DelegateTo>, DenseVectorRange<DT1_> >
                        wrapper(range);
                    std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(wrapper));
                    dispatched_tasks.push_back(ptr);
                }

                if (chunk_size != 0)
                {
                    for ( ; i < PARTS; ++i)
                    {
                        DenseVectorRange<DT1_> range((*vi).range(chunk_size, i *  chunk_size + rest));
                        OneArgWrapper< ElementInverse<typename Tag_::DelegateTo>, DenseVectorRange<DT1_> >
                            wrapper(range);
                        std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(wrapper));
                        dispatched_tasks.push_back(ptr);
                    }
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

        template <typename DT1_>
        static DenseVectorContinuousBase<DT1_> & value(DenseVectorContinuousBase<DT1_> & x)
        {
            CONTEXT("When calculating the inverse value of DenseVectorContinuousBase elements (MultiCore):");

            unsigned long min_part_size(Configuration::instance()->get_value("mc::element_inverse[DVCB]::min-part-size", 1024));
            unsigned long overall_size(x.size());

            if (overall_size < 2 * min_part_size)
            {
                ElementInverse<typename Tag_::DelegateTo>::value(x);
            }
            else
            {
                unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
                unsigned long max_count(Configuration::instance()->get_value("mc::element_inverse[DVCB]::max-count", num_threads));

                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, overall_size, PartitionList::Filler(partitions));
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()); p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = p->size;
                    DenseVectorRange<DT1_> range(x.range(part_size, offset));
                    OneArgWrapper< ElementInverse<typename Tag_::DelegateTo>, DenseVectorRange<DT1_> > mywrapper(range);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            return x;
        }


        template <typename DT1_>
        static DenseVectorSlice<DT1_> & value(DenseVectorSlice<DT1_> & x)
        {
            CONTEXT("When calculating the inverse value of DenseVectorSlice elements (MultiCore):");

            unsigned long min_part_size(Configuration::instance()->get_value("mc::element_inverse[DVS]::min-part-size", 64));
            unsigned long overall_size(x.used_elements());

            if (overall_size < 2 * min_part_size)
            {
                ElementInverse<typename Tag_::DelegateTo>::value(x);
            }
            else
            {
                unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
                unsigned long max_count(Configuration::instance()->get_value("mc::element_inverse[DVS]::max-count", num_threads));

                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, overall_size, PartitionList::Filler(partitions));
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;
                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()); p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = p->size;
                    typename Vector<DT1_>::ElementIterator start(x.begin_elements()), stop(x.begin_elements());
                    start += offset;
                    stop += (offset + part_size);
                    TwoArgWrapper< ElementInverse<Tag_>, typename Vector<DT1_>::ElementIterator, typename Vector<DT1_>::ElementIterator> mywrapper(start, stop);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            return x;
        }

        template <typename DT1_>
        static SparseVector<DT1_> value(SparseVector<DT1_> & x)
        {
            CONTEXT("When calculating the inverse value of SparseVector elements (MultiCore):");

            unsigned long min_part_size(Configuration::instance()->get_value("mc::element_inverse[SV]::min-part-size", 64));
            unsigned long overall_size(x.used_elements());

            if (overall_size < 2 * min_part_size)
            {
                ElementInverse<typename Tag_::DelegateTo>::value(x);
            }
            else
            {
                unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
                unsigned long max_count(Configuration::instance()->get_value("mc::element_inverse[SV]::max-count", num_threads));

                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, overall_size, PartitionList::Filler(partitions));
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;
                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()); p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = p->size;
                    typename Vector<DT1_>::ElementIterator start(x.begin_non_zero_elements()),
                             stop(x.begin_non_zero_elements());
                    start += offset;
                    stop += (offset + part_size);
                    TwoArgWrapper< ElementInverse<Tag_>, typename Vector<DT1_>::ElementIterator, typename Vector<DT1_>::ElementIterator> mywrapper(start, stop);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            return x;
        }

        template <typename IT1_, typename IT2_>
        static void value(IT1_ & x, const IT2_ & x_end)
        {
            CONTEXT("When calculating iterator-based element inverse:");

            typedef typeof(*x) DT1_;
            for ( ; x < x_end; ++x)
            {
                if (*x == DT1_(0))
                continue;
                *x = DT1_(1) / *x;
            }
        }
    };
}
#endif
