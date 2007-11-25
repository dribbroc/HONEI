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
#define PARTITION_SIZE 256

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

        static DenseVector<float> & value(const float a, DenseVector<float> & x);

        static DenseMatrix<float> & value(const float a, DenseMatrix<float> & x);

        /// \}
    };


    template <typename Tag_> struct MCScale
    {
        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT2_> & value(const DT1_ a, DenseMatrix<DT2_> & x)
        {
            CONTEXT("When scaling DenseMatrix (MultiCore):");

            ThreadPool * tp(ThreadPool::get_instance());

            PoolTask * pt[x.rows()];

            for (unsigned long i(0); i < x.rows(); ++i)
            {
                TwoArgWrapper< Scale<Tag_>, const DT1_, DenseVectorRange<DT2_> > wrapper(a, x[i]);
                pt[i] = tp->dispatch(wrapper);
            }

            for (unsigned long i(0); i < x.rows(); ++i)
            {
                pt[i]->wait_on();
            }

            return x;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT2_> & value(const DT1_ a, SparseMatrix<DT2_> & x)
        {
            CONTEXT("When scaling SparseMatrix (MultiCore):");

            ThreadPool * tp(ThreadPool::get_instance());

            PoolTask * pt[x.rows()];

            for (unsigned long i(0); i < x.rows(); ++i)
            {
                TwoArgWrapper< Scale<Tag_>, const DT1_, SparseVector<DT2_> > wrapper(a, x[i]);
                pt[i] = tp->dispatch(wrapper);
            }

            for (unsigned long i(0); i < x.rows(); ++i)
            {
                pt[i]->wait_on();
            }

            return x;
        }

        template <typename DT1_, typename DT2_>
        static BandedMatrix<DT2_> & value(const DT1_ a, BandedMatrix<DT2_> & x)
        {
            ///\todo: Avoid scaling padded zeros in bands.

            CONTEXT("When scaling BandedMatrix (MultiCore):");

            ThreadPool * tp(ThreadPool::get_instance());

            std::list< PoolTask* > dispatched_tasks;

            for (typename BandedMatrix<DT2_>::VectorIterator vi(x.begin_bands()), vi_end(x.end_bands());
                    vi != vi_end; ++vi)
            {
                if (! vi.exists())
                    continue;
                TwoArgWrapper< Scale<Tag_>, const DT1_, DenseVector<DT2_> > wrapper(a, *vi);
                dispatched_tasks.push_back(tp->dispatch(wrapper));
            }

            while(! dispatched_tasks.empty())
            {
                PoolTask * pt = dispatched_tasks.front();
                dispatched_tasks.pop_front();
                pt->wait_on();
            }

            return x;
        }

        template <typename DT1_, typename DT2_>
        static DenseVectorContinuousBase<DT2_> & value(const DT1_ a, DenseVectorContinuousBase<DT2_> & x)
        {
            CONTEXT("When scaling DenseVectorContinuousBase (MultiCore):");

            unsigned long partition_size(PARTITION_SIZE);

            unsigned long parts(x.size() / partition_size);

            ThreadPool * tp(ThreadPool::get_instance());

            PoolTask * pt[parts + (x.size() % partition_size ? 1 : 0)];

            for (unsigned long i(0); i < parts; ++i)
            {
                DenseVectorRange<DT2_> range(x.range(partition_size, i * partition_size));
                TwoArgWrapper< Scale<typename Tag_::DelegateTo>, const DT1_, DenseVectorRange<DT2_> >
                    wrapper(a, range);
                pt[i] = tp->dispatch(wrapper);
            }
            if (x.size() % partition_size)
            {
                DenseVectorRange<DT2_> range(x.range(x.size() % partition_size, parts * partition_size));
                TwoArgWrapper< Scale<typename Tag_::DelegateTo>, const DT1_,
                    DenseVectorRange<DT2_> > wrapper(a, range);
                pt[parts] = tp->dispatch(wrapper);
            }

            for (unsigned long i(0); i < parts; ++i)
            {
                pt[i]->wait_on();
            }

            if (x.size() % partition_size)
            {
                pt[parts]->wait_on();
            }

            return x;
        }

        template <typename DT1_, typename DT2_>
        static DenseVectorSlice<DT2_> & value(const DT1_ a, DenseVectorSlice<DT2_> & x)
        {
            CONTEXT("When scaling DenseVectorSlice (MultiCore):");

            unsigned long partition_size(PARTITION_SIZE);

            unsigned long parts(x.size() / partition_size);

            ThreadPool * tp(ThreadPool::get_instance());

            PoolTask * pt[parts + (x.size() % partition_size ? 1 : 0)];

            for (unsigned long i(0); i < parts; ++i)
            {
                typename Vector<DT2_>::ElementIterator start(x.begin_elements()), stop(x.begin_elements());
                start += i * partition_size;
                stop += (i + 1) * partition_size;
                ThreeArgWrapper< MCScale<Tag_>, const DT1_,
                    typename Vector<DT2_>::ElementIterator, typename Vector<DT2_>::ElementIterator > wrapper(a, start, stop);
                pt[i] = tp->dispatch(wrapper);
            }
            if (x.size() % partition_size)
            {
                typename Vector<DT2_>::ElementIterator start(x.begin_elements()), stop(x.end_elements());
                start += parts * partition_size;
                ThreeArgWrapper< MCScale<Tag_>, const DT1_,
                    typename Vector<DT2_>::ElementIterator, typename Vector<DT2_>::ElementIterator > wrapper(a, start, stop);
                pt[parts] = tp->dispatch(wrapper);
            }

            for (unsigned long i(0); i < parts; ++i)
            {
                pt[i]->wait_on();
            }

            if (x.size() % partition_size)
            {
                pt[parts]->wait_on();
            }

            return x;
        }

        template <typename DT1_, typename DT2_>
        static SparseVector<DT2_> value(const DT1_ a, SparseVector<DT2_> & x)
        {
            CONTEXT("When scaling SparseVector (MultiCore):");

            unsigned long partition_size(PARTITION_SIZE);
            unsigned long parts(x.used_elements() / partition_size);
            ThreadPool * tp(ThreadPool::get_instance());

            PoolTask * pt[(x.used_elements() % partition_size) ? parts + 1 : parts];

            for (unsigned long i(0); i < parts; ++i)
            {
                typename Vector<DT2_>::ElementIterator start(x.begin_non_zero_elements()), stop(x.begin_non_zero_elements());
                start += i * partition_size;
                stop += (i + 1) * partition_size;
                ThreeArgWrapper< MCScale<Tag_>, const DT1_,
                    typename Vector<DT2_>::ElementIterator, typename Vector<DT2_>::ElementIterator > wrapper(a, start, stop);
                pt[i] = tp->dispatch(wrapper);
            }

            if (x.used_elements() % partition_size)
            {
                typename Vector<DT2_>::ElementIterator start(x.begin_non_zero_elements()), stop(x.end_non_zero_elements());
                start += parts * partition_size;
                ThreeArgWrapper< MCScale<Tag_>, const DT1_,
                    typename Vector<DT2_>::ElementIterator, typename Vector<DT2_>::ElementIterator > wrapper(a, start, stop);
                pt[parts] = tp->dispatch(wrapper);
            }

            for (unsigned long i(0); i < parts; ++i)
            {
                pt[i]->wait_on();
            }

            if (x.used_elements() % partition_size) pt[parts]->wait_on();

            return x;
        }

        template <typename DT1_, typename IT1_, typename IT2_>
        static void value(const DT1_ a, IT1_ & x, const IT2_ & x_end)
        {
            CONTEXT("When calculating iterator-based scale (MultiCore):");

            for ( ; x < x_end; ++x)
            {
                *x *= a;
            }
        }

        #ifdef BENCHM
        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(DT1_ a, DenseMatrix<DT2_> & b)
        {
            BenchmarkInfo result;
            result.flops = b.rows() * b.columns();
            result.load = b.rows() * b.columns() * sizeof(DT2_) + sizeof(DT1_);
            result.store = b.rows() * b.columns() * sizeof(DT2_);
            return result; 
        }
        #endif
    };
    template <> struct Scale<tags::CPU::MultiCore> : public MCScale<tags::CPU::MultiCore> {};
    template <> struct Scale<tags::CPU::MultiCore::SSE> : public MCScale<tags::CPU::MultiCore::SSE> {};

}
#endif
