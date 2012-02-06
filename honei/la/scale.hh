/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007, 2008 Sven Mallach <sven.mallach@honei.org>
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

#pragma once
#ifndef LIBLA_GUARD_SCALE_HH
#define LIBLA_GUARD_SCALE_HH 1

#include <honei/la/banded_matrix.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/util/benchmark_info.hh>
#include <honei/util/tags.hh>
#include <honei/backends/multicore/operation.hh>
#include <honei/backends/multicore/thread_pool.hh>
#include <honei/util/configuration.hh>
#include <honei/util/operation_wrapper.hh>
#include <honei/util/partitioner.hh>
#include <honei/util/tags.hh>
#include <honei/mpi/operations.hh>
#include <honei/mpi/dense_vector_mpi-fwd.hh>

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
        static DenseMatrix<DT2_> & value(DenseMatrix<DT2_> & x, const DT1_ a)
        {
            CONTEXT("When scaling DenseMatrix:");
            for (typename DenseMatrix<DT2_>::ElementIterator l(x.begin_elements()),
                    l_end(x.end_elements()) ; l != l_end ; ++l)
            {
                *l *= a;
            }

            return x;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT2_> & value(SparseMatrix<DT2_> & x, const DT1_ a)
        {
            CONTEXT("When scaling SparseMatrix:");
            for (typename SparseMatrix<DT2_>::NonZeroElementIterator l(x.begin_non_zero_elements()),
                    l_end(x.end_non_zero_elements()) ; l != l_end ; ++l)
            {
                *l *= a;
            }

            return x;
        }

        template <typename DT1_, typename DT2_>
        static BandedMatrix<DT2_> & value(BandedMatrix<DT2_> & x, const DT1_ a)
        {
            CONTEXT("When scaling BandedMatrix:");
            for (typename BandedMatrix<DT2_>::BandIterator l(x.begin_non_zero_bands()),
                    l_end(x.end_non_zero_bands()) ; l != l_end ; ++l)
            {
                DenseVector<DT2_> band(*l);
                Scale<>::value(band, a);
            }

            return x;
        }

        template <typename DT1_, typename DT2_>
        static DenseVectorBase<DT2_> & value(DenseVectorBase<DT2_> & x, const DT1_ a)
        {
            CONTEXT("When scaling DenseVectorBase:");
            for (typename DenseVectorBase<DT2_>::ElementIterator l(x.begin_elements()),
                    l_end(x.end_elements()) ; l != l_end ; ++l)
            {
                *l *= a;
            }

            return x;
        }

        template <typename DT1_, typename DT2_>
        static inline DenseVectorContinuousBase<DT2_> & value(DenseVectorContinuousBase<DT2_> & x, const DT1_ a)
        {
            for (typename DenseVectorContinuousBase<DT2_>::ElementIterator l(x.begin_elements()),
                    l_end(x.end_elements()) ; l != l_end ; ++l)
            {
                *l *= a;
            }
            return x;
        }

        template <typename DT1_, typename DT2_>
        static inline DenseVectorRange<DT2_> & value(DenseVectorRange<DT2_> & x, const DT1_ a)
        {
            for (typename DenseVectorRange<DT2_>::ElementIterator l(x.begin_elements()),
                    l_end(x.end_elements()) ; l != l_end ; ++l)
            {
                *l *= a;
            }
            return x;
        }

        template <typename DT1_, typename DT2_>
        static inline DenseVectorSlice<DT2_> & value(DenseVectorSlice<DT2_> & x, const DT1_ a)
        {
            DenseVectorBase<DT2_> & temp = x;
            Scale<>::value(temp, a);
            return x;
        }

        template <typename DT1_, typename DT2_>
        static SparseVector<DT2_> & value(SparseVector<DT2_> & x, const DT1_ a)
        {
            CONTEXT("When scaling SparseVector:");
            for (typename SparseVector<DT2_>::NonZeroElementIterator l(x.begin_non_zero_elements()),
                    l_end(x.end_non_zero_elements()) ; l != l_end ; ++l)
            {
                *l *= a;
            }

            return x;
        }

        template <typename DT_>
        static inline DenseVectorMPI<DT_> & value(DenseVectorMPI<DT_> & x, DT_ a)
        {
            MPIOps<tags::CPU>::scale(x, a);
            return x;
        }

        /// \}

        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(DenseMatrix<DT2_> & b, HONEI_UNUSED DT1_ a)
        {
            BenchmarkInfo result;
            result.flops = b.rows() * b.columns();
            result.load = b.rows() * b.columns() * sizeof(DT2_) + sizeof(DT1_);
            result.store = b.rows() * b.columns() * sizeof(DT2_);
            result.size.push_back(b.rows() * b.columns());
            return result;
        }

        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(SparseMatrix<DT2_> & b, DT1_ a)
        {
            BenchmarkInfo result;
            for (typename SparseMatrix<DT2_>::NonZeroElementIterator l(b.begin_non_zero_elements()),
                    l_end(b.end_non_zero_elements()) ; l != l_end ; ++l)
            {
                result.flops += 1;
                result.load += sizeof(DT2_);
                result.store += sizeof(DT2_);
            }
            result.size.push_back(b.rows() * b.columns());
            result.load += sizeof(DT1_);
            result.scale = (double(b.rows() * b.columns()) / result.flops);
            return result;
        }

        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(DenseVectorBase<DT2_> & b, HONEI_UNUSED DT1_ a)
        {
            BenchmarkInfo result;
            result.flops = b.size();
            result.load = b.size() * sizeof(DT2_) + sizeof(DT1_);
            result.store = b.size() * sizeof(DT2_);
            result.size.push_back(b.size());
            return result;
        }

        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(DenseVectorContinuousBase<DT2_> & b, HONEI_UNUSED DT1_ a)
        {
            BenchmarkInfo result;
            result.flops = b.size();
            result.load = b.size() * sizeof(DT2_) + sizeof(DT1_);
            result.store = b.size() * sizeof(DT2_);
            result.size.push_back(b.size());
            return result;
        }
    };

    template <> struct Scale<tags::CPU::Generic>
    {
        template <typename DT_>
        static inline DenseVectorContinuousBase<DT_> & value(DenseVectorContinuousBase<DT_> & x, DT_ a)
        {
            BENCHADD(Scale<tags::CPU>::get_benchmark_info(x, a));

            DT_ * xe(x.elements());
            const unsigned long size(x.size());
            for (unsigned long i(0) ; i < size ; ++i)
            {
                xe[i] *= a;
            }
            return x;
        }
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
    template <> struct Scale<tags::GPU::CUDA>
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

        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & x, const float a);

        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & x, const double a);

        static DenseMatrix<float> & value(DenseMatrix<float> & x, const float a);

        /// \}
    };

    template <> struct Scale<tags::GPU::MultiCore::CUDA>
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

        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & x, const float a);

        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & x, const double a);

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
    template <> struct Scale<tags::CPU::SSE>
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

        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & x, const float a);

        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & x, const double a);

        static DenseMatrix<float> & value(DenseMatrix<float> & x, const float a);

        static DenseMatrix<double> & value(DenseMatrix<double> & x, const double a);

        static SparseVector<float> & value(SparseVector<float> & x, const float a);

        static SparseVector<double> & value(SparseVector<double> & x, const double a);

        /*static SparseMatrix<float> & value(SparseMatrix<float> & x, const float a);

        static SparseMatrix<double> & value(SparseMatrix<double> & x, const double a);*/

        template <typename DT_>
        static inline DenseVectorMPI<DT_> & value(DenseVectorMPI<DT_> & x, DT_ a)
        {
            MPIOps<tags::CPU::SSE>::scale(x, a);
            return x;
        }
        /// \}
    };

    template <> struct Scale<tags::OpenCL::CPU>
    {
        template <typename DT_>
        static DenseVectorContinuousBase<DT_> & value(DenseVectorContinuousBase<DT_> & x, const DT_ a);
    };

    template <> struct Scale<tags::OpenCL::GPU>
    {
        template <typename DT_>
        static DenseVectorContinuousBase<DT_> & value(DenseVectorContinuousBase<DT_> & x, const DT_ a);
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
    template <> struct Scale<tags::Cell>
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

        static DenseMatrix<float> & value(DenseMatrix<float> & x, const float a);

        static DenseMatrix<double> & value(DenseMatrix<double> & x, const double a);

        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & x, const float a);

        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & x, const double a);

        static SparseVector<float> & value(SparseVector<float> & x, const float a);

        static SparseMatrix<float> & value(SparseMatrix<float> & x, const float a);

        /// \}
    };

    namespace mc
    {
        template <typename Tag_> struct Scale
        {
            template <typename DT1_, typename DT2_>
            static DenseVectorContinuousBase<DT1_> & value(DenseVectorContinuousBase<DT1_> & x, const DT2_ a)
            {
                CONTEXT("When calculating Scale (DenseVectorContinuousBase) using backend : " + Tag_::name);

                unsigned long min_part_size(Configuration::instance()->get_value("mc::Scale(DVCB)::min_part_size", 128));
                unsigned long max_count(Configuration::instance()->get_value("mc::Scale(DVCB)::max_count",
                            mc::ThreadPool::instance()->num_threads()));

                Operation<honei::Scale<typename Tag_::DelegateTo> >::op(x, a, min_part_size, max_count);

                return x;
            }

            template <typename DT_>
            static inline DenseVectorMPI<DT_> & value(DenseVectorMPI<DT_> & x, const DT_ a)
            {
                MPIOps<Tag_>::scale(x, a);
                return x;
            }
        };
    }

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

    template <> struct Scale<tags::CPU::MultiCore> :
        public mc::Scale<tags::CPU::MultiCore>
    {
    };

    template <> struct Scale<tags::CPU::MultiCore::Generic> :
        public mc::Scale<tags::CPU::MultiCore::Generic>
    {
    };

    template <> struct Scale<tags::CPU::MultiCore::SSE> :
        public mc::Scale<tags::CPU::MultiCore::SSE>
    {
    };
}
#endif
