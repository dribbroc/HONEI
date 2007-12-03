/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Sven Mallach <sven.mallach@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_SCALE_HH
#define LIBLA_GUARD_SCALE_HH 1

#include <libla/dense_matrix.hh>
#include <libla/banded_matrix.hh>
#include <libla/sparse_matrix.hh>
#include <libutil/thread_pool.hh>
#include <libutil/wrapper.hh>
#include <libutil/tags.hh>

#include <iostream>
#include <list>


///\todo: Do not use define for setting size of multicore-partitions.
// For optimization purposes
#define MIN_CHUNK_SIZE 256
#define PARTS 8

namespace honei
{
    template <typename Tag_ = tags::CPU> struct Scale;

    /**
     * \brief Result of scaling an entity by a scalar factor.
     *
     * Scale is the class template for the operation
     * \f[
     *     \texttt{Scale}(a, x): \quad x[i] \leftarrow a \cdot x[i],
     * \f]
     * which yields the former entity scaled by a scalar factor.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <> struct Scale<tags::CPU>
    {
        /**
         * \name Scales
         * \{
         *
         * \brief Returns the product of a scalar factor and a given entity.
         *
         * \param a The scalar factor.
         * \param x The entity that shall be scaled.
         *
         * \retval x Will modify the entity x and return it.
         */

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT2_> & value(const DT1_ a, DenseMatrix<DT2_> & x)
        {
            CONTEXT("When scaling DenseMatrix");
            for (typename MutableMatrix<DT2_>::ElementIterator l(x.begin_elements()),
                    l_end(x.end_elements()) ; l != l_end ; ++l)
            {
                *l *= a;
            }

            return x;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT2_> & value(const DT1_ a, SparseMatrix<DT2_> & x)
        {
            CONTEXT("When scaling SparseMatrix");
            for (typename MutableMatrix<DT2_>::ElementIterator l(x.begin_non_zero_elements()),
                    l_end(x.end_non_zero_elements()) ; l != l_end ; ++l)
            {
                *l *= a;
            }

            return x;
        }

        template <typename DT1_, typename DT2_>
        static BandedMatrix<DT2_> & value(const DT1_ a, BandedMatrix<DT2_> & x)
        {
            CONTEXT("When scaling BandedMatrix");
            for (typename BandedMatrix<DT2_>::VectorIterator l(x.begin_bands()),
                    l_end(x.end_bands()) ; l != l_end ; ++l)
            {
                Scale<>::value(a, *l);
            }

            return x;
        }

        template <typename DT1_, typename DT2_>
        static DenseVectorBase<DT2_> & value(const DT1_ a, DenseVectorBase<DT2_> & x)
        {
            CONTEXT("When scaling DenseVectorBase");
            for (typename Vector<DT2_>::ElementIterator l(x.begin_elements()),
                    l_end(x.end_elements()) ; l != l_end ; ++l)
            {
                *l *= a;
            }

            return x;
        }

        template <typename DT1_, typename DT2_>
        static SparseVector<DT2_> & value(const DT1_ a, SparseVector<DT2_> & x)
        {
            CONTEXT("When scaling SparseVector");
            for (typename Vector<DT2_>::ElementIterator l(x.begin_non_zero_elements()),
                    l_end(x.end_non_zero_elements()) ; l != l_end ; ++l)
            {
                *l *= a;
            }

            return x;
        }

        /// \}

        #ifdef BENCHM
        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(DT1_ a, DenseMatrix<DT2_> & b)
        {
            BenchmarkInfo result;
            result.flops = b.rows() * b.columns();
            result.load = b.rows() * b.columns() * sizeof(DT2_) + sizeof(DT1_);
            result.store = b.rows() * b.columns() * sizeof(DT2_);
            result.size.push_back(b.rows() * b.columns());
            return result; 
        }
        #endif
    };

    /**
     * \brief Result of scaling an entity by a scalar factor.
     *
     * Scale is the class template for the operation
     * \f[
     *     \texttt{Scale}(a, x): \quad x[i] \leftarrow a \cdot x[i],
     * \f]
     * which yields the former entity scaled by a scalar factor.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <>
    struct Scale<tags::CPU::SSE>
    {
        /**
         * \name Scales
         * \{
         *
         * \brief Returns the product of a scalar factor and a given entity.
         *
         * \param a The scalar factor.
         * \param x The entity that shall be scaled.
         *
         * \retval x Will modify the entity x and return it.
         */

        static DenseVector<float> & value(const float a, DenseVector<float> & x);

        static DenseVector<double> & value(const double a, DenseVector<double> & x);

        static DenseMatrix<float> & value(const float a, DenseMatrix<float> & x);

        static DenseMatrix<double> & value(const double a, DenseMatrix<double> & x);

        static DenseVectorContinuousBase<float> & value(const float a, DenseVectorContinuousBase<float> & x);

        static DenseVectorContinuousBase<double> & value(const double a, DenseVectorContinuousBase<double> & x);

        /// \}
    };

    /**
     * \brief Result of scaling an entity by a scalar factor.
     *
     * Scale is the class template for the operation
     * \f[
     *     \texttt{Scale}(a, x): \quad x[i] \leftarrow a \cdot x[i],
     * \f]
     * which yields the former entity scaled by a scalar factor.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <>
    struct Scale<tags::Cell>
    {
        /**
         * \name Scales
         * \{
         *
         * \brief Returns the product of a scalar factor and a given entity.
         *
         * \param a The scalar factor.
         * \param x The entity that shall be scaled.
         *
         * \retval x Will modify the entity x and return it.
         */

        static DenseMatrix<float> & value(const float a, DenseMatrix<float> & x);

        static DenseVector<float> & value(const float a, DenseVector<float> & x);

        static SparseVector<float> & value(const float a, SparseVector<float> & x);

        /// \}
    };


    template <typename Tag_> struct MCScale
    {
        template <typename DT_>
        static inline DT_ calculate_parts(const DT_ ref)
        {
            if (ref / MIN_CHUNK_SIZE < PARTS) return ref / MIN_CHUNK_SIZE;
            return PARTS;
        }

        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT2_> & value(const DT1_ a, DenseMatrix<DT2_> & x)
        {
            CONTEXT("When scaling DenseMatrix (MultiCore):");

            ThreadPool * tp(ThreadPool::get_instance());

            std::list<PoolTask*> dispatched_tasks;

            unsigned long num_parts(MCScale<Tag_>::calculate_parts(x.columns()));

            if (num_parts)
            {

                unsigned long chunk_size(x.columns() / num_parts);
                unsigned long rest(x.columns() % chunk_size);

                for (unsigned long i(0); i < x.rows(); ++i)
                {
                    unsigned long j(0);
                    for ( ; j < rest; ++j)
                    {
                        DenseVectorRange<DT2_> range(x[i].range(chunk_size + 1, j * (chunk_size + 1)));
                        TwoArgWrapper< Scale<typename Tag_::DelegateTo>, const DT1_, DenseVectorRange<DT2_> > wrapper(a, range);
                        dispatched_tasks.push_back(tp->dispatch(wrapper));
                    }
                    for ( ; j < num_parts; ++j)
                    {
                        DenseVectorRange<DT2_> range(x[i].range(chunk_size, j * chunk_size + rest));
                        TwoArgWrapper< Scale<typename Tag_::DelegateTo>, const DT1_, DenseVectorRange<DT2_> > wrapper(a, range);
                        dispatched_tasks.push_back(tp->dispatch(wrapper));
                    }
                }
            }
            else
            {
                for (unsigned long i(0); i < x.rows(); ++i)
                {
                    TwoArgWrapper< Scale<typename Tag_::DelegateTo>, const DT1_, DenseVectorRange<DT2_> > wrapper(a, x[i]);
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

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT2_> & value(const DT1_ a, SparseMatrix<DT2_> & x)
        {
            CONTEXT("When scaling SparseMatrix (MultiCore):");

            ThreadPool * tp(ThreadPool::get_instance());

            std::list<PoolTask*> dispatched_tasks;

            for (unsigned long i(0); i < x.rows(); ++i)
            {
                unsigned long num_parts(MCScale<Tag_>::calculate_parts(x[i].used_elements()));
                if (num_parts)
                {
                    unsigned long chunk_size(x[i].used_elements() / num_parts);
                    unsigned long rest(x[i].used_elements() % chunk_size);
                    unsigned long j(0);
                    for ( ; j < rest; ++j)
                    {
                        typename Vector<DT2_>::ElementIterator start(x[i].begin_non_zero_elements()), stop(x[i].begin_non_zero_elements());
                        start += j * (chunk_size + 1);
                        stop += (j + 1) * (chunk_size + 1);
                        ThreeArgWrapper< MCScale<Tag_>, const DT1_, typename Vector<DT2_>::ElementIterator, typename Vector<DT2_>::ElementIterator >
                            wrapper(a, start, stop);
                        dispatched_tasks.push_back(tp->dispatch(wrapper));
                    }
                    for ( ; j < num_parts; ++j)
                    {
                        typename Vector<DT2_>::ElementIterator start(x[i].begin_non_zero_elements()), stop(x[i].begin_non_zero_elements());
                        start += j * chunk_size + rest;
                        stop += (j + 1) * chunk_size + rest;
                        ThreeArgWrapper< MCScale<Tag_>, const DT1_, typename Vector<DT2_>::ElementIterator, typename Vector<DT2_>::ElementIterator >
                            wrapper(a, start, stop);
                        dispatched_tasks.push_back(tp->dispatch(wrapper));
                    }
                }
                else
                {
                    typename Vector<DT2_>::ElementIterator start(x[i].begin_non_zero_elements()), stop(x[i].end_non_zero_elements());
                    ThreeArgWrapper< MCScale<Tag_>, const DT1_, typename Vector<DT2_>::ElementIterator, typename Vector<DT2_>::ElementIterator >
                        wrapper(a, start, stop);
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

        template <typename DT1_, typename DT2_>
        static BandedMatrix<DT2_> & value(const DT1_ a, BandedMatrix<DT2_> & x)
        {
            CONTEXT("When scaling BandedMatrix (MultiCore):");

            ThreadPool * tp(ThreadPool::get_instance());

            std::list< PoolTask* > dispatched_tasks;

            typename BandedMatrix<DT2_>::VectorIterator vi(x.begin_bands());

            // Calculating lower triangular matrix.
            for (typename BandedMatrix<DT2_>::VectorIterator vi(x.begin_bands()), vi_end(x.band_at(x.size() - 1)) ;
                    vi != vi_end ; ++vi)
            {
                if (! vi.exists())
                    continue;

                unsigned long band_size(x.size() - vi.index() - 1);
                unsigned long num_parts(MCScale<Tag_>::calculate_parts(band_size));

                if (num_parts)
                {
                    unsigned long chunk_size(band_size / num_parts);
                    unsigned long rest(band_size % chunk_size);
                    unsigned long i(0);

                    for ( ; i < rest; ++i)
                    {
                        DenseVectorRange<DT2_> range((*vi).range(chunk_size + 1, i * (chunk_size + 1) + x.size() - band_size));
                        TwoArgWrapper< Scale<typename Tag_::DelegateTo>, const DT1_, DenseVectorRange<DT2_> >
                            wrapper(a, range);
                        dispatched_tasks.push_back(tp->dispatch(wrapper));
                    }
                    for ( ; i < num_parts; ++i)
                    {
                        DenseVectorRange<DT2_> range((*vi).range(chunk_size, chunk_size + rest + x.size() - band_size));
                        TwoArgWrapper< Scale<typename Tag_::DelegateTo>, const DT1_, DenseVectorRange<DT2_> >
                            wrapper(a, range);
                        dispatched_tasks.push_back(tp->dispatch(wrapper));
                    }
                }
                else
                {
                    DenseVectorRange<DT2_> range((*vi).range(band_size, x.size() - band_size));
                    TwoArgWrapper< Scale<typename Tag_::DelegateTo>, const DT1_, DenseVectorRange<DT2_> > wrapper(a, range);
                    dispatched_tasks.push_back(tp->dispatch(wrapper));
                }
            }

            // Calculating diagonal band.
            if (vi.exists())
            {
                unsigned long num_parts(MCScale<Tag_>::calculate_parts(x.size()));

                if (num_parts)
                {
                    unsigned long chunk_size(x.size() / num_parts);
                    unsigned long rest(x.size() % chunk_size);
                    unsigned long i(0);

                    for ( ; i < rest; ++i)
                    {
                        DenseVectorRange<DT2_> range((*vi).range(chunk_size + 1, i * (chunk_size + 1)));
                        TwoArgWrapper< Scale<typename Tag_::DelegateTo>, const DT1_, DenseVectorRange<DT2_> >
                            wrapper(a, range);
                        dispatched_tasks.push_back(tp->dispatch(wrapper));
                    }
                    for ( ; i < num_parts; ++i)
                    {
                        DenseVectorRange<DT2_> range((*vi).range(chunk_size, i * chunk_size + rest));
                        TwoArgWrapper< Scale<typename Tag_::DelegateTo>, const DT1_, DenseVectorRange<DT2_> >
                            wrapper(a, range);
                        dispatched_tasks.push_back(tp->dispatch(wrapper));
                    }
                }
                else
                {
                    TwoArgWrapper< Scale<typename Tag_::DelegateTo>, const DT1_, DenseVector<DT2_> > wrapper(a, (*vi));
                    dispatched_tasks.push_back(tp->dispatch(wrapper));
                }
            }

            ++vi;

            // Calculating upper traingular matrix.
            for (typename BandedMatrix<DT2_>::VectorIterator vi_end(x.end_bands()) ; vi != vi_end ; ++vi)
            {
                if (! vi.exists())
                    continue;

                unsigned long band_size(2 * x.size() - vi.index() - 1);
                unsigned long num_parts(MCScale<Tag_>::calculate_parts(band_size));

                if (num_parts)
                {
                    unsigned long chunk_size(band_size / num_parts);
                    unsigned long rest(band_size % chunk_size);
                    unsigned long i(0);

                    for ( ; i < rest; ++i)
                    {
                        DenseVectorRange<DT2_> range((*vi).range(chunk_size + 1, i * (chunk_size + 1)));
                        TwoArgWrapper< Scale<typename Tag_::DelegateTo>, const DT1_, DenseVectorRange<DT2_> >
                            wrapper(a, range);
                        dispatched_tasks.push_back(tp->dispatch(wrapper));
                    }
                    for ( ; i < num_parts; ++i)
                    {
                        DenseVectorRange<DT2_> range((*vi).range(chunk_size, i *  chunk_size + rest));
                        TwoArgWrapper< Scale<typename Tag_::DelegateTo>, const DT1_, DenseVectorRange<DT2_> >
                            wrapper(a, range);
                        dispatched_tasks.push_back(tp->dispatch(wrapper));
                    }
                }
                else
                {
                    DenseVectorRange<DT2_> range((*vi).range(band_size, 0));
                    TwoArgWrapper< Scale<typename Tag_::DelegateTo>, const DT1_, DenseVectorRange<DT2_> >
                        wrapper(a, range);
                    dispatched_tasks.push_back(tp->dispatch(wrapper));
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
        static DenseVectorContinuousBase<DT2_> & value(const DT1_ a, DenseVectorContinuousBase<DT2_> & x)
        {
            CONTEXT("When scaling DenseVectorContinuousBase (MultiCore):");

            ThreadPool * tp(ThreadPool::get_instance());

            std::list<PoolTask *> dispatched_tasks;

            unsigned long num_parts(MCScale<Tag_>::calculate_parts(x.size()));

            if (num_parts)
            {
                unsigned long chunk_size(x.size() / num_parts);
                unsigned long rest(x.size() % chunk_size);

                unsigned long i(0);
                for ( ; i < rest; ++i)
                {
                    DenseVectorRange<DT2_> range(x.range(chunk_size + 1, i * (chunk_size + 1)));
                    TwoArgWrapper< Scale<typename Tag_::DelegateTo>, const DT1_, DenseVectorRange<DT2_> >
                        wrapper(a, range);
                    dispatched_tasks.push_back(tp->dispatch(wrapper));
                }
                for ( ; i < num_parts; ++i)
                {
                    DenseVectorRange<DT2_> range(x.range(chunk_size, i * chunk_size + rest));
                    TwoArgWrapper< Scale<typename Tag_::DelegateTo>, const DT1_,
                        DenseVectorRange<DT2_> > wrapper(a, range);
                    dispatched_tasks.push_back(tp->dispatch(wrapper));
                }

                while(! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            else
            {
                Scale<typename Tag_::DelegateTo>::value(a, x);
            }

            return x;
        }

        template <typename DT1_, typename DT2_>
        static DenseVectorSlice<DT2_> & value(const DT1_ a, DenseVectorSlice<DT2_> & x)
        {
            CONTEXT("When scaling DenseVectorSlice (MultiCore):");

            unsigned long num_parts(MCScale<Tag_>::calculate_parts(x.size()));
            if (num_parts)
            {
                unsigned long chunk_size(x.size() / num_parts);
                unsigned long rest(x.size() % chunk_size);

                ThreadPool * tp(ThreadPool::get_instance());

                std::list<PoolTask *> dispatched_tasks;

                unsigned long i(0);
                for ( ; i < rest; ++i)
                {
                    typename Vector<DT2_>::ElementIterator start(x.begin_elements()), stop(x.begin_elements());
                    start += i * (chunk_size + 1);
                    stop += (i + 1) * (chunk_size + 1);
                    ThreeArgWrapper< MCScale<Tag_>, const DT1_,
                        typename Vector<DT2_>::ElementIterator, typename Vector<DT2_>::ElementIterator > wrapper(a, start, stop);
                    dispatched_tasks.push_back(tp->dispatch(wrapper));
                }
                for ( ; i < num_parts; ++i)
                {
                    typename Vector<DT2_>::ElementIterator start(x.begin_elements()), stop(x.begin_elements());
                    start += i * chunk_size + rest;
                    stop += (i + 1) * chunk_size + rest;
                    ThreeArgWrapper< MCScale<Tag_>, const DT1_,
                        typename Vector<DT2_>::ElementIterator, typename Vector<DT2_>::ElementIterator > wrapper(a, start, stop);
                    dispatched_tasks.push_back(tp->dispatch(wrapper));
                }

                while (! dispatched_tasks.empty())
                {
                    dispatched_tasks.front()->wait_on();
                    dispatched_tasks.pop_front();
                }
            }
            else
            {
                Scale<>::value(a, x);
            }

            return x;
        }

        template <typename DT1_, typename DT2_>
        static SparseVector<DT2_> value(const DT1_ a, SparseVector<DT2_> & x)
        {
            CONTEXT("When scaling SparseVector (MultiCore):");

            unsigned long num_parts(MCScale<Tag_>::calculate_parts(x.used_elements()));

            if (num_parts)
            {
                unsigned long chunk_size(x.used_elements() / num_parts);
                unsigned long rest(x.used_elements() % chunk_size);
                ThreadPool * tp(ThreadPool::get_instance());

                std::list<PoolTask *> dispatched_tasks;

                unsigned long i(0);
                for ( ; i < rest; ++i)
                {
                    typename Vector<DT2_>::ElementIterator start(x.begin_non_zero_elements()), stop(x.begin_non_zero_elements());
                    start += i * (chunk_size + 1);
                    stop += (i + 1) * (chunk_size + 1);
                    ThreeArgWrapper< MCScale<Tag_>, const DT1_,
                        typename Vector<DT2_>::ElementIterator, typename Vector<DT2_>::ElementIterator > wrapper(a, start, stop);
                    dispatched_tasks.push_back(tp->dispatch(wrapper));
                }

                for ( ; i < num_parts; ++i)
                {
                    typename Vector<DT2_>::ElementIterator start(x.begin_non_zero_elements()), stop(x.begin_non_zero_elements());
                    start += i * chunk_size + rest;
                    stop += (i + 1) * chunk_size + rest;
                    ThreeArgWrapper< MCScale<Tag_>, const DT1_,
                        typename Vector<DT2_>::ElementIterator, typename Vector<DT2_>::ElementIterator > wrapper(a, start, stop);
                    dispatched_tasks.push_back(tp->dispatch(wrapper));
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
                MCScale<Tag_>::value(a, start, stop);
            }

            return x;
        }

        template <typename DT1_, typename IT1_, typename IT2_>
        static void value(const DT1_ a, IT1_ & x, const IT2_ & x_end)
        {
            CONTEXT("When calculating iterator-based scale:");

            for ( ; x < x_end; ++x)
            {
                *x *= a;
            }
        }
    };
    template <> struct Scale<tags::CPU::MultiCore> : public MCScale<tags::CPU::MultiCore> {};
    template <> struct Scale<tags::CPU::MultiCore::SSE> : public MCScale<tags::CPU::MultiCore::SSE> {};

}
#endif
