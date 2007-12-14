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

#include <libla/dense_vector_range.hh>
#include <libla/dense_vector_slice.hh>
#include <libla/banded_matrix.hh>
#include <libla/dense_matrix.hh>
#include <libla/dense_vector.hh>
#include <libla/element_inverse.hh>
#include <libla/sparse_matrix.hh>
#include <libla/sparse_vector.hh>
#include <libla/vector_error.hh>
#include <libutil/tags.hh>

#include <libla/dense_vector_range.hh>
#include <libutil/pool_task.hh>
#include <libutil/thread_pool.hh>
#include <libutil/wrapper.hh>

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

            ThreadPool * tp(ThreadPool::get_instance());

            unsigned long rest(x.columns() % PARTS);
            unsigned long chunk_size(x.columns() / PARTS);

            std::list<PoolTask*> dispatched_tasks;

            for (unsigned long i(0); i < x.rows(); ++i)
            {
                unsigned long j(0);
                for ( ; j < rest; ++j)
                {
                    DenseVectorRange<DT1_> range(x[i].range(chunk_size + 1, j * (chunk_size + 1)));
                    OneArgWrapper< ElementInverse<typename Tag_::DelegateTo>, DenseVectorRange<DT1_> > wrapper(range);
                    dispatched_tasks.push_back(tp->dispatch(wrapper));

                }
                if (chunk_size != 0)
                {
                    for ( ; j < PARTS; ++j)
                    {
                        DenseVectorRange<DT1_> range(x[i].range(chunk_size, j * chunk_size + rest));
                        OneArgWrapper< ElementInverse<typename Tag_::DelegateTo>, DenseVectorRange<DT1_> > wrapper(range);
                        dispatched_tasks.push_back(tp->dispatch(wrapper));
                    }
                }
            }

            while(! dispatched_tasks.empty())
            {
                dispatched_tasks.front()->wait_on();
                dispatched_tasks.pop_front();
            }

            return x;
        }

        template <typename DT1_>
        static SparseMatrix<DT1_> & value(SparseMatrix<DT1_> & x)
        {
            CONTEXT("When calculating the inverse value of SparseMatrix elements (MultiCore):");

            ThreadPool * tp(ThreadPool::get_instance());

            std::list<PoolTask*> dispatched_tasks;

            for (unsigned long i(0); i < x.rows(); ++i)
            {
                unsigned long chunk_size(x[i].used_elements() / PARTS);
                unsigned long rest(x[i].used_elements() % PARTS);
                unsigned long j(0);
                for ( ; j < rest; ++j)
                {
                    typename Vector<DT1_>::ElementIterator start(x[i].begin_non_zero_elements()), stop(x[i].begin_non_zero_elements());
                    start += j * (chunk_size + 1);
                    stop += (j + 1) * (chunk_size + 1);
                    TwoArgWrapper< MCElementInverse<Tag_>, typename Vector<DT1_>::ElementIterator, typename Vector<DT1_>::ElementIterator >
                        wrapper(start, stop);
                    dispatched_tasks.push_back(tp->dispatch(wrapper));
                }

                if (chunk_size != 0)
                {
                    for ( ; j < PARTS; ++j)
                    {
                        typename Vector<DT1_>::ElementIterator start(x[i].begin_non_zero_elements()), stop(x[i].begin_non_zero_elements());
                        start += j * chunk_size + rest;
                        stop += (j + 1) * chunk_size + rest;
                        TwoArgWrapper< MCElementInverse<Tag_>, typename Vector<DT1_>::ElementIterator, typename Vector<DT1_>::ElementIterator >
                            wrapper(start, stop);
                        dispatched_tasks.push_back(tp->dispatch(wrapper));
                    }
                }
            }

            while (! dispatched_tasks.empty())
            {
                dispatched_tasks.front()->wait_on();
                dispatched_tasks.pop_front();
            }

            return x;
        }

        template <typename DT1_>
        static BandedMatrix<DT1_> & value(BandedMatrix<DT1_> & x)
        {
            CONTEXT("When calculating the inverse value of BandedMatrix elements (MultiCore):");

            ThreadPool * tp(ThreadPool::get_instance());

            std::list< PoolTask* > dispatched_tasks;

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
                    dispatched_tasks.push_back(tp->dispatch(wrapper));
                }

                if (chunk_size != 0)
                {
                    for ( ; i < PARTS; ++i)
                    {
                        DenseVectorRange<DT1_> range((*vi).range(chunk_size, chunk_size + rest + x.size() - band_size));
                        OneArgWrapper< ElementInverse<typename Tag_::DelegateTo>, DenseVectorRange<DT1_> >
                            wrapper(range);
                        dispatched_tasks.push_back(tp->dispatch(wrapper));
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
                    dispatched_tasks.push_back(tp->dispatch(wrapper));
                }

                if (chunk_size != 0)
                {
                    for ( ; i < PARTS; ++i)
                    {
                        DenseVectorRange<DT1_> range((*vi).range(chunk_size, i *  chunk_size + rest));
                        OneArgWrapper< ElementInverse<typename Tag_::DelegateTo>, DenseVectorRange<DT1_> >
                            wrapper(range);
                        dispatched_tasks.push_back(tp->dispatch(wrapper));
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
                    dispatched_tasks.push_back(tp->dispatch(wrapper));
                }

                if (chunk_size != 0)
                {
                    for ( ; i < PARTS; ++i)
                    {
                        DenseVectorRange<DT1_> range((*vi).range(chunk_size, i *  chunk_size + rest));
                        OneArgWrapper< ElementInverse<typename Tag_::DelegateTo>, DenseVectorRange<DT1_> >
                            wrapper(range);
                        dispatched_tasks.push_back(tp->dispatch(wrapper));
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

            unsigned long rest(x.size() % PARTS);

            unsigned long chunk_size(x.size() / PARTS);

            ThreadPool * tp(ThreadPool::get_instance());

            std::list<PoolTask *> dispatched_tasks;

            unsigned long i(0);
            for ( ; i < rest; ++i)
            {
                DenseVectorRange<DT1_> range(x.range(chunk_size + 1, i * (chunk_size + 1)));
                OneArgWrapper< ElementInverse<typename Tag_::DelegateTo>, DenseVectorRange<DT1_> > wrapper(range);
                dispatched_tasks.push_back(tp->dispatch(wrapper));
            }

            if (chunk_size != 0)
            {
                for ( ; i < PARTS; ++i)
                {
                    DenseVectorRange<DT1_> range(x.range(chunk_size, i * chunk_size + rest));
                    OneArgWrapper< ElementInverse<typename Tag_::DelegateTo>, DenseVectorRange<DT1_> > wrapper(range);
                    dispatched_tasks.push_back(tp->dispatch(wrapper));
                }
            }

            while(! dispatched_tasks.empty())
            {
                dispatched_tasks.front()->wait_on();
                dispatched_tasks.pop_front();
            }

            return x;
        }

        template <typename DT1_>
        static DenseVectorSlice<DT1_> & value(DenseVectorSlice<DT1_> & x)
        {
            CONTEXT("When calculating the inverse value of DenseVectorSlice elements (MultiCore):");

            unsigned long rest(x.size() % PARTS);

            unsigned long chunk_size(x.size() / PARTS);

            ThreadPool * tp(ThreadPool::get_instance());

            std::list<PoolTask *> dispatched_tasks;

            unsigned long i(0);
            for ( ; i < rest; ++i)
            {
                typename Vector<DT1_>::ElementIterator start(x.begin_elements()), stop(x.begin_elements());
                start += i * (chunk_size + 1);
                stop += (i + 1) * (chunk_size + 1);
                TwoArgWrapper< MCElementInverse<Tag_>, typename Vector<DT1_>::ElementIterator, typename Vector<DT1_>::ElementIterator > wrapper(start, stop);
                dispatched_tasks.push_back(tp->dispatch(wrapper));
            }

            if (chunk_size != 0)
            {
                for ( ; i < PARTS; ++i)
                {
                    typename Vector<DT1_>::ElementIterator start(x.begin_elements()), stop(x.begin_elements());
                    start += i * chunk_size + rest;
                    stop += (i + 1) * chunk_size + rest;
                    TwoArgWrapper< MCElementInverse<Tag_>, typename Vector<DT1_>::ElementIterator, typename Vector<DT1_>::ElementIterator > wrapper(start, stop);
                    dispatched_tasks.push_back(tp->dispatch(wrapper));
                }
            }

            while (! dispatched_tasks.empty())
            {
                dispatched_tasks.front()->wait_on();
                dispatched_tasks.pop_front();
            }

            return x;
        }

        template <typename DT1_>
        static SparseVector<DT1_> value(SparseVector<DT1_> & x)
        {
            CONTEXT("When calculating the inverse value of SparseVector elements (MultiCore):");

            unsigned long rest(x.used_elements() % PARTS);
            unsigned long chunk_size(x.used_elements() / PARTS);
            ThreadPool * tp(ThreadPool::get_instance());

            std::list<PoolTask *> dispatched_tasks;

            unsigned long i(0);
            for ( ; i < rest; ++i)
            {
                typename Vector<DT1_>::ElementIterator start(x.begin_non_zero_elements()), stop(x.begin_non_zero_elements());
                start += i * (chunk_size + 1);
                stop += (i + 1) * (chunk_size + 1);
                TwoArgWrapper< MCElementInverse<Tag_>, typename Vector<DT1_>::ElementIterator, typename Vector<DT1_>::ElementIterator > wrapper(start, stop);
                dispatched_tasks.push_back(tp->dispatch(wrapper));
            }

            if (chunk_size != 0)
            {
                for ( ; i < PARTS; ++i)
                {
                    typename Vector<DT1_>::ElementIterator start(x.begin_non_zero_elements()), stop(x.begin_non_zero_elements());
                    start += i * chunk_size + rest;
                    stop += (i + 1) * chunk_size + rest;
                    TwoArgWrapper< MCElementInverse<Tag_>, typename Vector<DT1_>::ElementIterator, typename Vector<DT1_>::ElementIterator > wrapper(start, stop);
                    dispatched_tasks.push_back(tp->dispatch(wrapper));
                }
            }

            while (! dispatched_tasks.empty())
            {
                dispatched_tasks.front()->wait_on();
                dispatched_tasks.pop_front();
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
