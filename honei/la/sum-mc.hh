/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 André Matuschek <andre@matuschek.org>
 * Copyright (c) 2007 Joachim Messer <joachim.messer@uni-dortmund.de>
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

#define PARTS 8

#ifndef LIBLA_GUARD_SUM_MC_HH
#define LIBLA_GUARD_SUM_MC_HH 1

#include <honei/la/banded_matrix.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/util/tags.hh>
#include <honei/util/thread_pool.hh>
#include <honei/util/wrapper.hh>

#include <honei/util/configuration.hh>
#include <honei/util/partitioner.hh>


namespace honei
{
    // Forward declaration.
    template <typename Tag_> struct Sum;

    template <typename Tag_> struct MCSum
    {
        template <typename DT1_, typename DT2_>
        static void value(DenseVectorRange<DT1_> & a, const SparseVector<DT2_> & b, unsigned long offset)
        {
            typename Vector<DT2_>::ConstElementIterator r(b.begin_non_zero_elements());
            r += offset;
            offset = r.index();
            unsigned long limit = r.index() + a.size();
            while (r.index() < limit)
            {
                a[r.index()-offset] += *r;
                ++r;
            }
        }

        template <typename DT1_, typename DT2_>
        static DenseVectorContinuousBase<DT1_> & value(DenseVectorContinuousBase<DT1_> & a, const DenseVectorContinuousBase<DT2_> & b)
        {
            CONTEXT("When adding DenseVectorContinuousBase to DenseVectorContinuousBase (MultiCore):");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());
            unsigned long min_part_size(Configuration::instance()->get_value("mc::sum[DVCB,DVCB]::min-part-size", 1024));
            unsigned long overall_size(a.size());

            if (overall_size < 2 * min_part_size)
            {
                Sum<typename Tag_::DelegateTo>::value(a,b);
            }
            else
            {
                unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
                unsigned long max_count(Configuration::instance()->get_value("mc::sum[DVCB,DVCB]::max-count", num_threads));

                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, overall_size, PartitionList::Filler(partitions));
                //ThreadPool * pool(ThreadPool::instance());
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()); p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = p->size;
                    DenseVectorRange<DT1_> range_1(a.range(part_size, offset));
                    DenseVectorRange<DT2_> range_2(b.range(part_size, offset));
                    TwoArgWrapper<Sum<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>, const DenseVectorRange<DT2_> > mywrapper(range_1, range_2);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            return a;
        }

        template <typename DT1_, typename DT2_>
        static DenseVectorContinuousBase<DT1_> & value(DenseVectorContinuousBase<DT1_> & a, const SparseVector<DT2_> & b)
        {
            CONTEXT("When adding DenseVectorContinuousBase to SparseVector (MultiCore):");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());
            unsigned long min_part_size(Configuration::instance()->get_value("mc::sum[DVCB,SV]::min-part-size", 1024));
            unsigned long overall_size(b.used_elements());
            if (overall_size < 2 * min_part_size)
            {
                Sum<typename Tag_::DelegateTo>::value(a, b);
            }
            else
            {
                unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
                unsigned long max_count(Configuration::instance()->get_value("mc::sum[DVCB,DVCB]::max-count", num_threads ));

                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, overall_size, PartitionList::Filler(partitions));
                //ThreadPool * pool(ThreadPool::instance());
                std::list<std::tr1::shared_ptr<PoolTask> > dispatched_tasks;
                typename Vector<DT2_>::ConstElementIterator r(b.begin_non_zero_elements());

                unsigned long offset, part_size;

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()) ; p != p_end ; ++p)
                {
                    part_size = p->size;
                    offset = r.index();
                    r += (p->size - 1);
                    DenseVectorRange<DT1_> range(a.range(r.index() - offset + 1, offset));
                    ThreeArgWrapper<Sum<Tag_>, DenseVectorRange<DT1_>, const SparseVector<DT2_>, const unsigned long > mywrapper(range, b, p->start);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                    ++r;
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            return a;
        }

        template <typename DT1_>
        static DenseVectorContinuousBase<DT1_> & value(DenseVectorContinuousBase<DT1_> & a, const DT1_ b)
        {
            CONTEXT("When adding Scalar to DenseVectorContinuousBase (MultiCore):");
            unsigned long min_part_size(Configuration::instance()->get_value("mc::sum[DVCB,SC]::min-part-size", 1024));
            unsigned long overall_size(a.size());

            if (overall_size < 2 * min_part_size)
            {
                Sum<typename Tag_::DelegateTo>::value(a,b);
            }
            else
            {
                unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
                unsigned long max_count(Configuration::instance()->get_value("mc::sum[DVCB,SC]::max-count", num_threads));

                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, overall_size, PartitionList::Filler(partitions));
                //ThreadPool * pool(ThreadPool::instance());
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()); p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = p->size;

                    DenseVectorRange<DT1_> range_1(a.range(part_size, offset));
                    TwoArgWrapper<Sum<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>, const DT1_ > mywrapper(range_1, b);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            return a;
        }

        template <typename DT1_, typename DT2_>
        static BandedMatrix<DT1_> & value(BandedMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
        {
            CONTEXT("When adding BandedMatrix to BandedMatrix (MultiCore):");

            if (a.size() != b.size())
            {
                throw MatrixSizeDoesNotMatch(b.size(), a.size());
            }

            //ThreadPool * p(ThreadPool::instance());
            PoolTask * pt[2*a.rows()-1];
            int taskcount(0);
            typename BandedMatrix<DT1_>::VectorIterator l(a.begin_bands()), l_end(a.end_bands());
            typename BandedMatrix<DT2_>::ConstVectorIterator r(b.begin_bands()), r_end(b.end_bands());
            for ( ; ((l != l_end) && (r != r_end)) ; ++l, ++r)
            {
                if (! r.exists())
                    continue;

                if (l.exists())
                {
                    TwoArgWrapper< Sum<typename Tag_::DelegateTo>, DenseVector<DT1_>, const DenseVector<DT2_> > mywrapper(*l, *r);
                    pt[taskcount] = ThreadPool::instance()->dispatch(mywrapper);
                    ++taskcount;
                }
                else
                    a.band(r.index()) = r->copy();
            }
            for (unsigned long j = 0; j < taskcount; ++j)
            {
                pt[j]->wait_on();
                delete pt[j];
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
        {
            CONTEXT("When adding DenseMatrix to DenseMatrix (MultiCore):");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }
            
            unsigned long min_part_size(Configuration::instance()->get_value("mc::sum[DM,DM]::min-part-size", 64));
            unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
            unsigned long max_count(Configuration::instance()->get_value("mc::sum[DM,DM]::max-count", num_threads));
            unsigned long overall_size(a.rows());
            if ((overall_size < max_count) && ((a.columns() / min_part_size) > overall_size) && ((a.columns() / min_part_size) >= 2))
            {
                for (int i(0) ; i < a.rows() ; ++i)
                {
                    Sum<Tag_>::value(a[i], b[i]);
                }
            }
            else if (overall_size < 2)
            {
                Sum<typename Tag_::DelegateTo>::value(a, b);
            }
            else
            {
                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, 1, 1, overall_size, PartitionList::Filler(partitions));
                //ThreadPool * pool(ThreadPool::instance());
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()) ; p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = p->size;

                    FourArgWrapper< Sum<Tag_>, DenseMatrix<DT1_>, const DenseMatrix<DT2_>, unsigned long, unsigned long> 
                        mywrapper(a, b, offset, part_size);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            return a;
        }
    
        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b, unsigned long offset, unsigned long part_size)
        {
            for (unsigned long i(offset) ; i < (offset + part_size) ; ++i)
            {
                Sum<typename Tag_::DelegateTo>::value(a[i], b[i]);
            }
        }

        template <typename DT1_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const DT1_ b)
        {
            CONTEXT("When adding Scalar to DenseMatrix (MultiCore):");

            unsigned long min_part_size(Configuration::instance()->get_value("mc::sum[DM,SC]::min-part-size", 64));
            unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
            unsigned long max_count(Configuration::instance()->get_value("mc::sum[DM,SC]::max-count", num_threads));
            unsigned long overall_size(a.rows());
            if ((overall_size < max_count) && ((a.columns() / min_part_size) > overall_size) && ((a.columns() / min_part_size) >= 2))
            {
                for (int i(0) ; i < a.rows() ; ++i)
                {
                    Sum<Tag_>::value(a[i], b);
                }
            }
            else if (overall_size < 2)
            {
                Sum<typename Tag_::DelegateTo>::value(a, b);
            }
            else
            {
                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, 1, 1, overall_size, PartitionList::Filler(partitions));
                //ThreadPool * pool(ThreadPool::instance());
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()) ; p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = p->size;

                    FourArgWrapper< Sum<Tag_>, DenseMatrix<DT1_>, const DT1_, unsigned long, unsigned long> 
                        mywrapper(a, b, offset, part_size);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            return a;
        }

        template <typename DT1_>
        static void value(DenseMatrix<DT1_> & a, const DT1_ b, unsigned long offset, unsigned long part_size)
        {
            for (unsigned long i(offset) ; i < (offset + part_size) ; ++i)
            {
                Sum<typename Tag_::DelegateTo>::value(a[i], b);
            }
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
        {
            CONTEXT("When adding DenseMatrix to SparseMatrix (MultiCore):");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            unsigned long min_part_size(Configuration::instance()->get_value("mc::sum[DM,SM]::min-part-size", 64));
            unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
            unsigned long max_count(Configuration::instance()->get_value("mc::sum[DM,SM]::max-count", num_threads));
            unsigned long overall_size(a.rows());
            if ((overall_size < 2) || (a.columns() < min_part_size))
            {
                Sum<typename Tag_::DelegateTo>::value(a,b);
            }
            else
            {
                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, 1, 1, overall_size, PartitionList::Filler(partitions));
                //ThreadPool * pool(ThreadPool::instance());
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()) ; p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = p->size;

                    FourArgWrapper< Sum<Tag_>, DenseMatrix<DT1_>, const SparseMatrix<DT2_>, unsigned long, unsigned long> 
                        mywrapper(a, b, offset, part_size);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            return a;
        }

        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b, unsigned long offset, unsigned long part_size)
        {
            for (unsigned long i(offset) ; i < (offset + part_size) ; ++i)
            {
                Sum<typename Tag_::DelegateTo>::value(a[i], b[i]);
            }
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
        {
            CONTEXT("When adding SparseMatrix to SparseMatrix (MultiCore):");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            unsigned long min_part_size(Configuration::instance()->get_value("mc::sum[SM,SM]::min-part-size", 64));
            unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
            unsigned long max_count(Configuration::instance()->get_value("mc::sum[SM,SM]::max-count", num_threads));
            unsigned long overall_size(a.rows());
            if ((overall_size < 2) || (a.columns() < min_part_size))
            {
                Sum<typename Tag_::DelegateTo>::value(a,b);
            }
            else
            {
                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, 1, 1, overall_size, PartitionList::Filler(partitions));
                //ThreadPool * pool(ThreadPool::instance());
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()) ; p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = p->size;

                    FourArgWrapper< Sum<Tag_>, SparseMatrix<DT1_>, const SparseMatrix<DT2_>, unsigned long, unsigned long> 
                        mywrapper(a, b, offset, part_size);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            return a;
        }

        template <typename DT1_, typename DT2_>
        static void value(SparseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b, unsigned long offset, unsigned long part_size)
        {
            for (unsigned long i(offset) ; i < (offset + part_size) ; ++i)
            {
                Sum<typename Tag_::DelegateTo>::value(a[i], b[i]);
            }
        }


        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
        {
            CONTEXT("When adding BandedMatrix to DenseMatrix (MutiCore):");

            if (a.columns() != a.rows())
            {
                throw MatrixIsNotSquare(a.rows(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            unsigned long min_part_size(Configuration::instance()->get_value("mc::sum[DM,BM]::min-part-size", 1024));
            unsigned long overall_size(b.size());

            if (overall_size < 2 * min_part_size)
            {
                Sum<typename Tag_::DelegateTo>::value(a, b);
            }

            else
            {
                unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
                unsigned long max_count(Configuration::instance()->get_value("mc::sum[DM,BM]::max-count", num_threads ));
                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, overall_size, PartitionList::Filler(partitions));
                //ThreadPool * pool(ThreadPool::instance());
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()) ; p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = offset + p->size - 1;

                    FourArgWrapper<Sum<Tag_>, DenseMatrix<DT2_>, const BandedMatrix<DT1_>, unsigned long, unsigned long>
                        mywrapper(a, b, offset, part_size);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(mywrapper));
                    dispatched_tasks.push_back(ptr);
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            return a;
        }

        template <typename DT1_, typename DT2_>
        static void value(DenseMatrix<DT1_> & a, const BandedMatrix<DT2_> & b, unsigned long start, unsigned long end)
        {
            CONTEXT("When partial adding BandedMatrix to DenseMatrix:");

            unsigned long size(b.size());
            for (typename BandedMatrix<DT2_>::ConstVectorIterator r(b.begin_non_zero_bands()), r_end(b.end_non_zero_bands()) ;
                    r != r_end ; ++r)
            {

                if (r.index() < b.size()-1)
                { // lower part.
                    unsigned long row_index(std::max(long(-(r.index() - size + 1)), long(0)));
                    unsigned long col_index(std::max(long(r.index() - size + 1), long(0)));
                    typename Vector<DT2_>::ConstElementIterator c(r->begin_elements()), c_end(r->end_elements());

                    c += ((size-1) - r.index());

                    for ( ; c != c_end ; ++c)
                    {
                        if (row_index > end)
                            break;

                        if (row_index >= start)
                        {
                            a[row_index][col_index] += *c;
                        }

                        ++row_index;
                        ++col_index;
                    }
                }
                else
                { // upper part.
                    unsigned long size(b.size());
                    unsigned long row_index(std::max(long(-(r.index() - size + 1)), long(0)));
                    unsigned long col_index(std::max(long(r.index() - size + 1), long(0)));
                    typename Vector<DT2_>::ConstElementIterator c(r->begin_elements()), c_end(r->end_elements());
                    c += start;
                    row_index += start;
                    col_index += start;

                    for ( ; c != c_end ; ++c)
                    {
                        if (row_index > end)
                            break;

                        if (row_index >= size)
                            break;

                        if (col_index >= size)
                            break;

                        if (row_index >= start)
                        {
                            a[row_index][col_index] += *c;
                        }

                        ++row_index;
                        ++col_index;
                    }
                }
            }
        }
    };
}
#endif
