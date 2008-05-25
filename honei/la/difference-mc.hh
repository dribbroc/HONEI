/* vim: set sw=4 sts=4 et nofoldenable : */

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

#ifndef LIBLA_GUARD_DIFFERENCE_MC_HH
#define LIBLA_GUARD_DIFFERENCE_MC_HH 1

#include <honei/la/banded_matrix.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/matrix_error.hh>
#include <honei/la/product.hh>
#include <honei/la/scale.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/vector.hh>
#include <honei/util/configuration.hh>
#include <honei/util/partitioner.hh>
#include <honei/util/tags.hh>

#include <list>
namespace honei
{
    template <typename Tag_> struct Difference;

    template <typename Tag_> struct MCDifference
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
                a[r.index()-offset] -= *r;
                ++r;
            }
        }

        template <typename DT1_, typename DT2_>
        static void value(const DenseVectorRange<DT1_> & a, SparseMatrix<DT2_> & b, unsigned long offset_y, unsigned long offset_x)
        {
            for(typename Vector<DT1_>::ConstElementIterator c(a.begin_elements()),
                        c_end(a.end_elements()) ; c != c_end ; ++c)
                {
                    b[offset_y][offset_x] += *c;
                    ++offset_y, ++offset_x;
                }
        }

        template <typename DT1_, typename DT2_>
        static DenseVectorContinuousBase<DT1_> & value(DenseVectorContinuousBase<DT1_> & a, const DenseVectorContinuousBase<DT2_> & b)
        {
            CONTEXT("When substracting DenseVectorContinuousBase from DenseVectorContinuousBase (MultiCore):");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());
            unsigned long min_part_size(Configuration::instance()->get_value("mc::difference[DVCB,DVCB]::min-part-size", 1024));
            unsigned long overall_size(a.size());

            if (overall_size < 2 * min_part_size)
            {
                Difference<typename Tag_::DelegateTo>::value(a,b);
            }
            else
            {
                unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
                unsigned long max_count(Configuration::instance()->get_value("mc::difference[DVCB,DVCB]::max-count", num_threads));

                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, overall_size, PartitionList::Filler(partitions));
                //ThreadPool * pool(ThreadPool::instance());
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()) ; p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = p->size;
                    DenseVectorRange<DT1_> range_1(a.range(part_size, offset));
                    DenseVectorRange<DT2_> range_2(b.range(part_size, offset));
                    TwoArgWrapper<Difference<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>, const DenseVectorRange<DT2_> > mywrapper(range_1, range_2);
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
            CONTEXT("When substracting DenseVectorContinuousBase from SparseVector (MultiCore):");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());
            unsigned long min_part_size(Configuration::instance()->get_value("mc::difference[DVCB,SV]::min-part-size", 1024));
            unsigned long overall_size(b.used_elements());
            if (overall_size < 2 * min_part_size)
            {
                Difference<typename Tag_::DelegateTo>::value(a, b);
            }
            else
            {
                unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
                unsigned long max_count(Configuration::instance()->get_value("mc::difference[DVCB,DVCB]::max-count", num_threads ));

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
                    ThreeArgWrapper<MCDifference<Tag_>, DenseVectorRange<DT1_>, const SparseVector<DT2_>, const unsigned long > mywrapper(range, b, p->start);
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

        template <typename DT1_, typename DT2_>
        static BandedMatrix<DT1_> & value(BandedMatrix<DT1_> & a, const BandedMatrix<DT2_> & b)
        {
            CONTEXT("When substracting BandedMatrix from BandedMatrix (MultiCore):");

            if (a.size() != b.size())
            {
                throw MatrixSizeDoesNotMatch(b.size(), a.size());
            }

            unsigned long min_part_size(Configuration::instance()->get_value("mc::difference[BM,BM]::min-part-size", 1024));
            if (a.size() < 2 * min_part_size)
            {
                Difference<typename Tag_::DelegateTo>::value(a, b);
            }
            else
            {
                //ThreadPool * p(ThreadPool::instance());
                PoolTask * pt[2*a.rows()-1];
                int taskcount(0);
                typename BandedMatrix<DT1_>::BandIterator l(a.begin_non_zero_bands()), l_end(a.end_non_zero_bands());
                typename BandedMatrix<DT2_>::ConstBandIterator r(b.begin_non_zero_bands()), r_end(b.end_non_zero_bands());
                unsigned long a_index, b_index;
                while((l != l_end) && (r != r_end))
                {
                    a_index = l.index();
                    b_index = r.index();
                    if (a_index < b_index)
                    {
                        ++l;
                    }

                    if (a_index == b_index)
                    {
                        DenseVector<DT1_> band(*l);
                        TwoArgWrapper< Difference<typename Tag_::DelegateTo>, DenseVector<DT1_>, const DenseVector<DT2_> > mywrapper(band, *r);
                        pt[taskcount] = ThreadPool::instance()->dispatch(mywrapper);
                        ++taskcount;
                        ++l;
                        ++r;
                    }
                    else
                    {
                        DenseVector<DT2_> band(r->copy());
                        Scale<typename Tag_::DelegateTo>::value(band, DT1_(-1));
                        a.band(r.index()-a.size()+1) = band;
                        ++r;
                    }
                }
                for (unsigned long j = 0; j < taskcount; ++j)
                {
                    pt[j]->wait_on();
                    delete pt[j];
                }
            }
            return a;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(const BandedMatrix<DT1_> & a, DenseMatrix<DT2_> & b)
        {
            CONTEXT("When subtracting DenseMatrix from BandedMatrix (MultiCore):");

            if (b.columns() != b.rows())
            {
                throw MatrixIsNotSquare(a.rows(), a.columns());
            }

            if (b.rows() != a.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            unsigned long min_part_size(Configuration::instance()->get_value("mc::difference[DM,BM]::min-part-size", 1024));
            unsigned long overall_size(a.size());

            if (overall_size < 2 * min_part_size)
            {
                Difference<typename Tag_::DelegateTo>::value(a, b);
            }

            else
            {
                Scale<Tag_>::value(b, -1);

                unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
                unsigned long max_count(Configuration::instance()->get_value("mc::difference[DM,BM]::max-count", num_threads ));
                PartitionList partitions;
                Partitioner<tags::CPU::MultiCore>(max_count, min_part_size, 16, overall_size, PartitionList::Filler(partitions));
                //ThreadPool * pool(ThreadPool::instance());
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

                unsigned long offset, part_size;

                for (PartitionList::ConstIterator p(partitions.begin()), p_end(partitions.end()) ; p != p_end ; ++p)
                {
                    offset = p->start;
                    part_size = offset + p->size - 1;

                    FourArgWrapper<MCDifference<typename Tag_::DelegateTo>, const BandedMatrix<DT1_>,
                        DenseMatrix<DT2_>, unsigned long, unsigned long>
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

            return b;
        }

        template <typename DT1_, typename DT2_>
        static void value(const BandedMatrix<DT1_> & a, DenseMatrix<DT2_> & b, unsigned long start, unsigned long end)
        {
            CONTEXT("When partial substracting BandedMatrix to DenseMatrix:");

            unsigned long size(a.size());
            for (typename BandedMatrix<DT1_>::ConstBandIterator r(a.begin_non_zero_bands()), r_end(a.end_non_zero_bands()) ;
                    r != r_end ; ++r)
            {

                if (r.index() < size - 1)
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
                            b[row_index][col_index] += *c;
                        }

                        ++row_index;
                        ++col_index;
                    }
                }
                else
                { // upper part.
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
                            b[row_index][col_index] += *c;
                        }

                        ++row_index;
                        ++col_index;
                    }
                }
            }
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const DenseMatrix<DT2_> & b)
        {
            CONTEXT("When substracting DenseMatrix from DenseMatrix (MultiCore):");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            unsigned long min_part_size(Configuration::instance()->get_value("mc::difference[DM,DM]::min-part-size", 64));
            unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
            unsigned long max_count(Configuration::instance()->get_value("mc::difference[DM,DM]::max-count", num_threads));
            unsigned long overall_size(a.rows());
            if ((overall_size < max_count) && ((a.columns() / min_part_size) > overall_size) && ((a.columns() / min_part_size) >= 2))
            {
                for (int i(0) ; i < a.rows() ; ++i)
                {
                    Difference<Tag_>::value(a[i], b[i]);
                }
            }
            else if (overall_size < 2)
            {
                Difference<typename Tag_::DelegateTo>::value(a, b);
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

                    FourArgWrapper< Difference<Tag_>, DenseMatrix<DT1_>, const DenseMatrix<DT2_>, unsigned long, unsigned long> 
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
                Difference<typename Tag_::DelegateTo>::value(a[i], b[i]);
            }
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(DenseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
        {
            CONTEXT("When substracting DenseMatrix from SparseMatrix (MultiCore):");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            unsigned long min_part_size(Configuration::instance()->get_value("mc::difference[DM,SM]::min-part-size", 64));
            unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
            unsigned long max_count(Configuration::instance()->get_value("mc::difference[DM,SM]::max-count", num_threads));
            unsigned long overall_size(a.rows());
            if ((overall_size < 2) || (a.columns() < min_part_size))
            {
                Difference<typename Tag_::DelegateTo>::value(a,b);
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

                    FourArgWrapper< Difference<Tag_>, DenseMatrix<DT1_>, const SparseMatrix<DT2_>, unsigned long, unsigned long> 
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
                Difference<typename Tag_::DelegateTo>::value(a[i], b[i]);
            }
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & a, const SparseMatrix<DT2_> & b)
        {
            CONTEXT("When substracting SparseMatrix from SparseMatrix (MultiCore):");

            if (a.columns() != b.columns())
            {
                throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            unsigned long min_part_size(Configuration::instance()->get_value("mc::difference[SM,SM]::min-part-size", 64));
            unsigned long num_threads(2 * Configuration::instance()->get_value("mc::num-cores", 2));
            unsigned long max_count(Configuration::instance()->get_value("mc::difference[SM,SM]::max-count", num_threads));
            unsigned long overall_size(a.rows());
            if ((overall_size < 2) || (a.columns() < min_part_size))
            {
                Difference<typename Tag_::DelegateTo>::value(a,b);
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

                    FourArgWrapper< Difference<Tag_>, SparseMatrix<DT1_>, const SparseMatrix<DT2_>, unsigned long, unsigned long> 
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
                Difference<typename Tag_::DelegateTo>::value(a[i], b[i]);
            }
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT1_> & value(const BandedMatrix<DT2_> & a, SparseMatrix<DT1_> & b)
        {
            CONTEXT("When subtracting SparseMatrix from BandedMatrix (MultiCore):");

            if (b.columns() != b.rows())
            {
                throw MatrixIsNotSquare(b.rows(), b.columns());
            }

            if (a.rows() != b.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            unsigned long min_part_size(Configuration::instance()->get_value("mc::difference[BM,SM]::min-part-size", 1024));
            unsigned long parts(8);
            unsigned long div = a.size() / parts;
            if (a.size() < min_part_size)
            {
                return Difference<typename Tag_::DelegateTo>::value(a,b);
            }
            else
            {
                Scale<Tag_>::value(b, -1);

                unsigned long modulo = a.size() % parts;
                //ThreadPool * tp(ThreadPool::instance());
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;
                Mutex mutex[parts];
                int middle_index(a.rows() -1);
                //if we are below the diagonal band
                for (typename BandedMatrix<DT2_>::ConstBandIterator vi(a.begin_non_zero_bands()),
                     vi_end(a.band_at(middle_index)) ; vi != vi_end ; ++vi)
                {
                    unsigned long i(parts), offset(a.size());
                    unsigned long start(middle_index - vi.index());
                    while(i > modulo && offset-div > start)
                    {
                        --i;
                        offset-=div;
                        DenseVectorRange<DT2_> range_1(*vi, div, offset);
                        FourArgWrapper<MCDifference<Tag_>, const DenseVectorRange<DT2_>, SparseMatrix<DT1_>, const unsigned long, const unsigned long > mywrapper(range_1, b, offset, (vi.index()+offset)-a.size()+1);
                        std::tr1::function<void ()> func = std::tr1::bind(mywrapper, &mutex[i]);
                        std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(func));
                        dispatched_tasks.push_back(ptr);
                    }
                    if (i == modulo)
                    {
                        while(i > 0 && offset-div-1 > start)
                        {
                            --i;
                            offset-=(div+1);
                            DenseVectorRange<DT2_> range_1(*vi, div+1, offset);
                            FourArgWrapper<MCDifference<Tag_>, DenseVectorRange<DT2_>, SparseMatrix<DT1_>, const unsigned long, const unsigned long > mywrapper(range_1, b, offset, (vi.index()+offset)-a.size()+1);
                            std::tr1::function<void ()> func = std::tr1::bind(mywrapper, &mutex[i]);
                            std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(func));
                            dispatched_tasks.push_back(ptr);
                        }
                    }
                    if (offset > start)
                    {
                        DenseVectorRange<DT2_> range_1(*vi, offset-start, start);
                        FourArgWrapper<MCDifference<Tag_>, DenseVectorRange<DT2_>, SparseMatrix<DT1_>, const unsigned long, const unsigned long > mywrapper(range_1, b, start, 0);
                        std::tr1::function<void ()> func = std::tr1::bind(mywrapper, &mutex[i-1]);
                        std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(func));
                        dispatched_tasks.push_back(ptr);
                    }
                }
                // If we are above or on the diagonal band
                for (typename BandedMatrix<DT2_>::ConstBandIterator vi(a.band_at(middle_index)),
                     vi_end(a.end_non_zero_bands()); vi != vi_end ; ++vi)
                {
                    unsigned long i(0), offset(0);
                    unsigned long index(vi.index() - middle_index);
                    unsigned long end(vi->size() - index);
                    while ((i < modulo) && (offset+div+1 < end))
                    {
                        DenseVectorRange<DT2_> range_1(*vi, div+1, offset);
                        FourArgWrapper<MCDifference<Tag_>, DenseVectorRange<DT2_>, SparseMatrix<DT1_>, const unsigned long, const unsigned long > mywrapper(range_1, b, offset, index + offset);
                        std::tr1::function<void ()> func = std::tr1::bind(mywrapper, &mutex[i]);
                        std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(func));
                        dispatched_tasks.push_back(ptr);
                        ++i;
                        offset+=div+1;
                    }
                    if (i == modulo)
                    {
                        while ((i < parts) && (offset+div  < end))
                        {
                            DenseVectorRange<DT2_> range_1(*vi, div, offset);
                            FourArgWrapper<MCDifference<Tag_>, DenseVectorRange<DT2_>, SparseMatrix<DT1_>, const unsigned long, const unsigned long > mywrapper(range_1, b, offset, index+offset);
                            std::tr1::function<void ()> func = std::tr1::bind(mywrapper, &mutex[i]);
                            std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(func));
                            dispatched_tasks.push_back(ptr);
                            ++i;
                            offset+=div;
                        }
                    }
                    if (offset < end)
                    {
                        DenseVectorRange<DT2_> range_1(*vi, end - offset, offset);
                        FourArgWrapper<MCDifference<Tag_>, DenseVectorRange<DT2_>, SparseMatrix<DT1_>, const unsigned long, const unsigned long > mywrapper(range_1, b, offset, index+offset);
                        std::tr1::function<void ()> func = std::tr1::bind(mywrapper, &mutex[i]);
                        std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(func));
                        dispatched_tasks.push_back(ptr);
                    }
                }

                while(! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
                return b;
            }
        }
    };
}
#endif
