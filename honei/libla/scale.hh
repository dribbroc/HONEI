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

#ifndef LIBLA_GUARD_SCALE_HH
#define LIBLA_GUARD_SCALE_HH 1

#include <honei/libla/dense_matrix.hh>
#include <honei/libla/banded_matrix.hh>
#include <honei/libla/scale-mc.hh>
#include <honei/libla/sparse_matrix.hh>
#include <honei/libutil/thread_pool.hh>
#include <honei/libutil/wrapper.hh>
#include <honei/libutil/tags.hh>

#include <list>


///\todo: Do not use define for setting size of multicore-partitions.
// For optimization purposes
#define MIN_CHUNK_SIZE 256
#define PARTS 8

namespace honei
{
//    template <typename Tag_> struct MCScale;

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
            CONTEXT("When scaling DenseMatrix");
            for (typename MutableMatrix<DT2_>::ElementIterator l(x.begin_elements()),
                    l_end(x.end_elements()) ; l != l_end ; ++l)
            {
                *l *= a;
            }

            return x;
        }

        template <typename DT1_, typename DT2_>
        static SparseMatrix<DT2_> & value(SparseMatrix<DT2_> & x, const DT1_ a)
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
        static BandedMatrix<DT2_> & value(BandedMatrix<DT2_> & x, const DT1_ a)
        {
            CONTEXT("When scaling BandedMatrix");
            for (typename BandedMatrix<DT2_>::VectorIterator l(x.begin_bands()),
                    l_end(x.end_bands()) ; l != l_end ; ++l)
            {
                Scale<>::value(*l, a);
            }

            return x;
        }

        template <typename DT1_, typename DT2_>
        static DenseVectorBase<DT2_> & value(DenseVectorBase<DT2_> & x, const DT1_ a)
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
        static inline DenseVector<DT2_> & value(DenseVector<DT2_> & x, const DT1_ a)
        {
            DenseVectorBase<DT2_> & temp = x;
            Scale<>::value(temp, a);
            return x;
        }

        template <typename DT1_, typename DT2_>
        static inline DenseVectorRange<DT2_> & value(DenseVectorRange<DT2_> & x, const DT1_ a)
        {
            DenseVectorBase<DT2_> & temp = x;
            Scale<>::value(temp, a);
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

        template <typename DT1_, typename DT2_>
        static inline BenchmarkInfo get_benchmark_info(DT1_ a, SparseMatrix<DT2_> & b)
        {
            BenchmarkInfo result;
            for (typename MutableMatrix<DT2_>::ElementIterator l(b.begin_non_zero_elements()),
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
        static inline BenchmarkInfo get_benchmark_info(DT1_ a, DenseVector<DT2_> & b)
        {
            BenchmarkInfo result;
            result.flops = b.size();
            result.load = b.size() * sizeof(DT2_) + sizeof(DT1_);
            result.store = b.size() * sizeof(DT2_);
            result.size.push_back(b.size());
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

        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & x, const float a);

        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & x, const double a);

        static DenseMatrix<float> & value(DenseMatrix<float> & x, const float a);

        static DenseMatrix<double> & value(DenseMatrix<double> & x, const double a);

        static SparseVector<float> & value(SparseVector<float> & x, const float a);

        static SparseVector<double> & value(SparseVector<double> & x, const double a);

        static SparseMatrix<float> & value(SparseMatrix<float> & x, const float a);

        static SparseMatrix<double> & value(SparseMatrix<double> & x, const double a);
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

        static DenseMatrix<float> & value(DenseMatrix<float> & x, const float a);

        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & x, const float a);

        static SparseVector<float> & value(SparseVector<float> & x, const float a);

        static SparseMatrix<float> & value(SparseMatrix<float> & x, const float a);

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
    template <> struct Scale<tags::CPU::MultiCore> : public MCScale<tags::CPU::MultiCore> {};
    template <> struct Scale<tags::CPU::MultiCore::SSE> : public MCScale<tags::CPU::MultiCore::SSE> {};

}
#endif