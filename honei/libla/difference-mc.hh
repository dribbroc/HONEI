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

#include <honei/libla/banded_matrix.hh>
#include <honei/libla/dense_matrix.hh>
#include <honei/libla/matrix_error.hh>
#include <honei/libla/product.hh>
#include <honei/libla/scale.hh>
#include <honei/libla/sparse_matrix.hh>
#include <honei/libla/dense_vector.hh>
#include <honei/libla/vector.hh>
#include <honei/libutil/tags.hh>


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
        static DenseVector<DT1_> & value(DenseVector<DT1_> & a, const DenseVector<DT2_> & b)
        {
            CONTEXT("When substracting DenseVector from DenseVector (MultiCore):");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());
            unsigned long parts(8);
            unsigned long div = a.size() / parts;
            if (div == 0)
            {
                Difference<typename Tag_::DelegateTo>::value(a,b);
            }
            else
            {
                unsigned long modulo = a.size() % parts;
                ThreadPool * p(ThreadPool::get_instance());
                PoolTask * pt[parts];
                for (int i(0); i < modulo; ++i)
                {
                    DenseVectorRange<DT1_> range_1(a, div+1, i*(div+1));
                    DenseVectorRange<DT2_> range_2(b, div+1, i*(div+1));
                    TwoArgWrapper<Difference<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>, const DenseVectorRange<DT2_> > mywrapper(range_1, range_2);
                    pt[i] = p->dispatch(mywrapper);
                }
                for (unsigned long i(modulo); i < parts; ++i)
                {
                    DenseVectorRange<DT1_> range_1(a, div, modulo+(i*div));
                    DenseVectorRange<DT2_> range_2(b, div, modulo+(i*div));
                    TwoArgWrapper<Difference<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>, const DenseVectorRange<DT2_> > mywrapper(range_1, range_2);
                    pt[i] = p->dispatch(mywrapper);
                }
                for (unsigned long i(0); i < parts; ++i)
                {
                    pt[i]->wait_on();
                    delete pt[i];
                }
            }
            return a;
        }

        template <typename DT1_, typename DT2_>
        static DenseVector<DT1_> & value(DenseVector<DT1_> & a, const SparseVector<DT2_> & b)
        {
            CONTEXT("When substracting DenseVector from SparseVector (MultiCore):");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());
            unsigned long parts(8);
            unsigned long modulo = b.used_elements() % parts;
            unsigned long div = b.used_elements() / parts;
            if (div == 0)
            {
                Difference<typename Tag_::DelegateTo>::value(a, b);
            }
            else
            {
                ThreadPool * p(ThreadPool::get_instance());
                PoolTask * pt[parts];
                typename Vector<DT2_>::ConstElementIterator r(b.begin_non_zero_elements());
                unsigned long offset;
                for (int i(0); i < modulo; ++i)
                {
                    offset = r.index();
                    r += div;
                    DenseVectorRange<DT1_> range(a, r.index()-offset+1, offset);
                    ThreeArgWrapper<MCDifference<Tag_>, DenseVectorRange<DT1_>, const SparseVector<DT2_>, const unsigned long > mywrapper(range, b, (i*(div+1)));
                    pt[i] = p->dispatch(mywrapper);
                    ++r;
                }
                for (unsigned long i(modulo); i < parts; ++i)
                {
                    offset = r.index();
                    r+= div-1;
                    DenseVectorRange<DT1_> range(a, r.index()-offset+1, offset);
                    ThreeArgWrapper<MCDifference<Tag_>, DenseVectorRange<DT1_>, const SparseVector<DT2_>, const unsigned long > mywrapper(range, b, modulo + (i*div));
                    pt[i] = p->dispatch(mywrapper);
                    ++r;
                }
                for (unsigned long i = 0; i < parts;  ++i)
                {
                    pt[i]->wait_on();
                    delete pt[i];
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

            ThreadPool * p(ThreadPool::get_instance());
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
                    TwoArgWrapper< Difference<typename Tag_::DelegateTo>, DenseVector<DT1_>, const DenseVector<DT2_> > mywrapper(*l, *r);
                    pt[taskcount] = p->dispatch(mywrapper);
                    ++taskcount;
                }
                else
                {
                   DenseVector<DT2_> band(r->copy());
                   Scale<typename Tag_::DelegateTo>::value(band, DT1_(-1));
                   a.band(r.index()-a.size()+1) = band;
                }
            }
            for (unsigned long j = 0; j < taskcount; ++j)
            {
                pt[j]->wait_on();
                delete pt[j];
            }

            return a;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT1_> & value(const BandedMatrix<DT1_> & a, DenseMatrix<DT2_> & b)
        {
            CONTEXT("When subtracting DenseMatrix from BandedMatrix (MultiCore):");

            unsigned long parts(PARTS);
            unsigned long modulo(0);
            unsigned long div(1);
            unsigned long start(0);
            unsigned long end(0);

            if (b.columns() != b.rows())
            {
                throw MatrixIsNotSquare(a.rows(), a.columns());
            }

            if (b.rows() != a.rows())
            {
                throw MatrixRowsDoNotMatch(b.rows(), a.rows());
            }

            Scale<Tag_>::value(b, -1);

            if (a.size() < parts)
            {
                parts = a.size();
            }
            else
            {
                div = a.size() / parts;
                modulo = a.size() % parts;
            }

            ThreadPool * p(ThreadPool::get_instance());
            PoolTask * pt[parts];

            for (int i(0) ; i < modulo ; ++i)
            { // (div+1) elements, from (i * (div + 1)) to ((i+1) * (div + 1) - 1).
                start = (i * (div + 1));
                end   = ((i + 1) * (div + 1) - 1);
                FourArgWrapper<MCDifference<typename Tag_::DelegateTo>, const BandedMatrix<DT1_>,
                    DenseMatrix<DT2_>, unsigned long, unsigned long>
                    mywrapper(a, b, start, end);
                pt[i] = p->dispatch(mywrapper);
           }

            for (int i(modulo) ; i < parts ; ++i)
            { // div elements, from (modulo + div * i) to (modulo + div * (i+1) - 1).
                start = (modulo + (div * i));
                end   = (modulo + (div * (i +1) - 1));
                FourArgWrapper< MCDifference<typename Tag_::DelegateTo>, const BandedMatrix<DT1_>,
                    DenseMatrix<DT2_>, unsigned long, unsigned long>
                    mywrapper(a, b, start, end);
                pt[i] = p->dispatch(mywrapper);
            }

            for (int i(0) ; i < parts ; ++i)
            {
                pt[i]->wait_on();
                delete pt[i];
            }

            return b;
        }

        template <typename DT1_, typename DT2_>
        static void value(const BandedMatrix<DT1_> & a, DenseMatrix<DT2_> & b, unsigned long start, unsigned long end)
        {
            CONTEXT("When partial substracting BandedMatrix to DenseMatrix:");

            unsigned long size(a.size());
            for (typename BandedMatrix<DT1_>::ConstVectorIterator r(a.begin_non_zero_bands()), r_end(a.end_non_zero_bands()) ;
                    r != r_end ; ++r)
            {

                if (r.index() < a.size()-1)
                { // lower part.
                    unsigned long row_index(std::max(long(-(r.index() - size + 1)), long(0)));
                    unsigned long col_index(std::max(long(r.index() - size + 1), long(0)));
                    typename Vector<DT2_>::ConstElementIterator c(r->begin_elements()), c_end(r->end_elements());

                    c += ((size-1) - r.index());

                    for ( ; c != c_end ; ++c)
                    {
                        if (row_index > end)
                            break;

                        if ((row_index >= start) && (row_index <= end))
                        {
                            b[row_index][col_index] += *c;
                        }

                        ++row_index;
                        ++col_index;
                    }
                }
                else
                { // upper part.
                    unsigned long size(a.size());
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

                        if ((row_index >= start) && (row_index <= end))
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

            ThreadPool * p(ThreadPool::get_instance());
            PoolTask * pt[a.rows()];
            for (unsigned long i = 0 ; i < a.rows() ; ++i)
            {
                TwoArgWrapper< Difference<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>, const DenseVectorRange<DT2_> > mywrapper(a[i], b[i]);
                pt[i] = p->dispatch(mywrapper);
            }
            for (unsigned long i = 0; i < a.rows(); ++i)
            {
                pt[i]->wait_on();
                delete pt[i];
            }
            return a;
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

            ThreadPool * p(ThreadPool::get_instance());
            PoolTask * pt[a.rows()];
            for (unsigned long i = 0 ; i < a.rows() ; ++i)
            {
                TwoArgWrapper< Difference<typename Tag_::DelegateTo>, DenseVectorRange<DT1_>, const SparseVector<DT2_> > mywrapper(a[i], b[i]);
                pt[i] = p->dispatch(mywrapper);
            }
            for (unsigned long i = 0; i < a.rows(); ++i)
            {
                pt[i]->wait_on();
                delete pt[i];
            }
            return a;
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

            ThreadPool * p(ThreadPool::get_instance());
            PoolTask * pt[a.rows()];
            for (unsigned long i = 0 ; i < a.rows() ; ++i)
            {
                TwoArgWrapper< Difference<typename Tag_::DelegateTo>, SparseVector<DT1_>, const SparseVector<DT2_> > mywrapper(a[i], b[i]);
                pt[i] = p->dispatch(mywrapper);
            }
            for (unsigned long i = 0; i < a.rows(); ++i)
            {
                pt[i]->wait_on();
                delete pt[i];
            }
            return a;
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

            Scale<Tag_>::value(b, -1);

            unsigned long parts(8);
            unsigned long div = a.size() / parts;
            if (div == 0)
            {
                return Difference<typename Tag_::DelegateTo>::value(a,b);
            }
            else
            {
                unsigned long modulo = a.size() % parts;
                ThreadPool * tp(ThreadPool::get_instance());
                std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;
                Mutex mutex[parts];
                int middle_index(a.rows() -1);
                //if we are below the diagonal band
                for (typename BandedMatrix<DT2_>::ConstVectorIterator vi(a.begin_bands()),
                     vi_end(a.band_at(middle_index)) ; vi != vi_end ; ++vi)
                {
                    if (!vi.exists())
                        continue;
                    unsigned long i(parts), offset(a.size());
                    unsigned long start(middle_index - vi.index());
                    while(i > modulo && offset-div > start)
                    {
                        --i;
                        offset-=div;
                        DenseVectorRange<DT2_> range_1(*vi, div, offset);
                        FourArgWrapper<MCDifference<Tag_>, const DenseVectorRange<DT2_>, SparseMatrix<DT1_>, const unsigned long, const unsigned long > mywrapper(range_1, b, offset, (vi.index()+offset)-a.size()+1);
                        std::tr1::function<void ()> func = std::tr1::bind(mywrapper, &mutex[i]);
                        std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(func));
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
                            std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(func));
                            dispatched_tasks.push_back(ptr);
                        }
                    }
                    if (offset > start)
                    {
                        DenseVectorRange<DT2_> range_1(*vi, offset-start, start);
                        FourArgWrapper<MCDifference<Tag_>, DenseVectorRange<DT2_>, SparseMatrix<DT1_>, const unsigned long, const unsigned long > mywrapper(range_1, b, start, 0);
                        std::tr1::function<void ()> func = std::tr1::bind(mywrapper, &mutex[i-1]);
                        std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(func));
                        dispatched_tasks.push_back(ptr);
                    }
                }
                // If we are above or on the diagonal band
                for (typename BandedMatrix<DT2_>::ConstVectorIterator vi(a.band_at(middle_index)),
                     vi_end(a.end_bands()); vi != vi_end ; ++vi)
                {
                    if (!vi.exists())
                        continue;
                    unsigned long i(0), offset(0);
                    unsigned long index(vi.index() - middle_index);
                    unsigned long end(vi->size() - index);
                    while ((i < modulo) && (offset+div+1 < end))
                    {
                        DenseVectorRange<DT2_> range_1(*vi, div+1, offset);
                        FourArgWrapper<MCDifference<Tag_>, DenseVectorRange<DT2_>, SparseMatrix<DT1_>, const unsigned long, const unsigned long > mywrapper(range_1, b, offset, index + offset);
                        std::tr1::function<void ()> func = std::tr1::bind(mywrapper, &mutex[i]);
                        std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(func));
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
                            std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(func));
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
                        std::tr1::shared_ptr<PoolTask> ptr(tp->dispatch(func));
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
