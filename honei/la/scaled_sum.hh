/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Sven Mallach <sven.mallach@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_SCALED_SUM_HH
#define LIBLA_GUARD_SCALED_SUM_HH 1

#include <honei/la/dense_matrix.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/scaled_sum-mc.hh>
#include <honei/la/sparse_vector.hh>
#include <honei/la/vector_error.hh>
#include <honei/util/tags.hh>
#include <honei/util/benchmark_info.hh>

namespace honei
{
    template <typename Tag_ = tags::CPU> struct ScaledSum;

    /**
     * \brief Scaled sum of two given vectors and a given scalar.
     *
     * ScaledSum is the class template for the operation
     * \f[
     *     \texttt{ScaledSum}(x, a, y): \quad x \leftarrow a \cdot x + y,
     *     \texttt{ScaledSum}(x, y, b): \quad x \leftarrow x + b \cdot y,
     * \f]
     * which yields the scaled sum of x and y.
     *
     * \todo Implement variant using a.
     *
     * \ingroup grplaoperations
     * \ingroup grplavectoroperations
     */
    template <> struct ScaledSum<tags::CPU>
    {
        /**
         * \name Scaled sums
         * \{
         *
         * \brief Returns the vector x as the scaled sum of two given vectors.
         *
         * \param x The vector that shall not be scaled.
         * \param y The vector that shall be scaled.
         * \param b The scale factor.
         *
         * \retval x Will modify x and return it.
         *
         * \exception VectorSizeDoesNotMatch is thrown if the sizes of x and y do not match.
         */

        template <typename DT1_, typename DT2_>
        static DenseVectorBase<DT1_> & value(DenseVectorBase<DT1_> & x, const DenseVectorBase<DT2_> & y, DT2_ b)
        {
            CONTEXT("When calculating ScaledSum (DenseVectorBase, DenseVectorBase, scalar):");

            if (x.size() != y.size())
                throw VectorSizeDoesNotMatch(y.size(), x.size());

            typename Vector<DT2_>::ConstElementIterator r(y.begin_elements());
            for (typename Vector<DT1_>::ElementIterator l(x.begin_elements()),
                    l_end(x.end_elements()) ; l != l_end ; ++l, ++r)
            {
                *l += b * (*r);

            }

            return x;
        }

        template <typename DT1_, typename DT2_>
        static inline DenseVector<DT1_> & value(DenseVector<DT1_> & x, const DenseVectorBase<DT2_> & y, DT2_ b)
        {
            DenseVectorBase<DT1_> & temp = x;
            ScaledSum<>::value(temp, y, b);
            return x;
        }

        template <typename DT1_, typename DT2_>
        static inline DenseVectorRange<DT1_> & value(DenseVectorRange<DT1_> & x, const DenseVectorBase<DT2_> & y, DT2_ b)
        {
            DenseVectorBase<DT1_> & temp = x;
            ScaledSum<>::value(temp, y, b);
            return x;
        }

        template <typename DT1_, typename DT2_>
        static inline DenseVectorSlice<DT1_> & value(DenseVectorSlice<DT1_> & x, const DenseVectorBase<DT2_> & y, DT2_ b)
        {
            DenseVectorBase<DT1_> & temp = x;
            ScaledSum<>::value(temp, y, b);
            return x;
        }

        template <typename DT1_, typename DT2_, typename DT3_>
        static SparseVector<DT1_> & value(SparseVector<DT1_> & x, const SparseVector<DT2_> & y, DT3_ b)
        {
            CONTEXT("When calculating ScaledSum (SparseVector, SparseVector, scalar):");

            if (x.size() != y.size())
                throw VectorSizeDoesNotMatch(y.size(), x.size());

            typename Vector<DT2_>::ElementIterator l(x.begin_non_zero_elements());
            for (typename Vector<DT1_>::ConstElementIterator r(y.begin_non_zero_elements()),
                    r_end(y.end_non_zero_elements()) ; r != r_end ; )
            {
                if (r.index() < l.index())
                {
                    x[r.index()] = b * (*r);
                    ++r;
                }
                else if (l.index() < r.index())
                {
                    ++l;
                }
                else
                {
                    *l += b * (*r);
                    ++l; ++r;
                }
            }
            ///\todo: perhaps sparsify - written results may be zero.
            return x;
        }

        template <typename DT1_, typename DT2_, typename DT3_>
        static DenseVectorBase<DT1_> & value(DenseVectorBase<DT1_> & x, const SparseVector<DT2_> & y, DT3_ b)
        {
            CONTEXT("When calculating ScaledSum (DenseVectorBase, SparseVector, scalar):");

            if (x.size() != y.size())
                throw VectorSizeDoesNotMatch(y.size(), x.size());

            for (typename Vector<DT2_>::ConstElementIterator r(y.begin_non_zero_elements()),
                    r_end(y.end_non_zero_elements()) ; r != r_end ; ++r )
            {
                x[r.index()] += b * (*r);
            }
            ///\todo: perhaps sparsify - if senseless use with b == zero, zeros may be written.
            return x;
        }

        template <typename DT1_, typename DT2_, typename DT3_>
        static inline DenseVector<DT1_> & value(DenseVector<DT1_> & x, const SparseVector<DT2_> & y, DT3_ b)
        {
            DenseVectorBase<DT1_> & temp = x;
            ScaledSum<>::value(temp, y, b);
            return x;
        }

        template <typename DT1_, typename DT2_, typename DT3_>
        static inline DenseVectorRange<DT1_> & value(DenseVectorRange<DT1_> & x, const SparseVector<DT2_> & y, DT3_ b)
        {
            DenseVectorBase<DT1_> & temp = x;
            ScaledSum<>::value(temp, y, b);
            return x;
        }

        template <typename DT1_, typename DT2_, typename DT3_>
        static inline DenseVectorSlice<DT1_> & value(DenseVectorSlice<DT1_> & x, const SparseVector<DT2_> & y, DT3_ b)
        {
            DenseVectorBase<DT1_> & temp = x;
            ScaledSum<>::value(temp, y, b);
            return x;
        }

        template <typename DT1_, typename DT2_, typename DT3_>
        static DenseVectorBase<DT1_> & value(DenseVectorBase<DT1_> & a, const DenseVectorBase<DT2_> & b, const DenseVectorBase<DT3_> & c)
        {
            CONTEXT("When calculating ScaledSum (DenseVectorBase, DenseVectorBase, DenseVectorBase):");

            if (a.size() != b.size())
                throw VectorSizeDoesNotMatch(b.size(), a.size());
            if (a.size() != c.size())
                throw VectorSizeDoesNotMatch(c.size(), a.size());

            typename Vector<DT2_>::ConstElementIterator l(b.begin_elements());
            typename Vector<DT3_>::ConstElementIterator r(c.begin_elements());
            for (typename Vector<DT1_>::ElementIterator s(a.begin_elements()),
                    s_end(a.end_elements()) ; s != s_end ; ++l)
            {
                *s += *l * *r;
                ++r, ++s;
            }

            return a;
        }

        template <typename DT1_, typename DT2_, typename DT3_>
        static inline DenseVector<DT1_> & value(DenseVector<DT1_> & x, const DenseVectorBase<DT2_> & y, const DenseVectorBase<DT3_> & z)
        {
            DenseVectorBase<DT1_> & temp = x;
            ScaledSum<>::value(temp, y, z);
            return x;
        }

        template <typename DT1_, typename DT2_, typename DT3_>
        static inline DenseVectorRange<DT1_> & value(DenseVectorRange<DT1_> & x, const DenseVectorBase<DT2_> & y, const DenseVectorBase<DT3_> & z)
        {
            DenseVectorBase<DT1_> & temp = x;
            ScaledSum<>::value(temp, y, z);
            return x;
        }

        template <typename DT1_, typename DT2_, typename DT3_>
        static inline DenseVectorSlice<DT1_> & value(DenseVectorSlice<DT1_> & x, const DenseVectorBase<DT2_> & y, const DenseVectorBase<DT3_> & z)
        {
            DenseVectorBase<DT1_> & temp = x;
            ScaledSum<>::value(temp, y, z);
            return x;
        }
        /// \}

        template <typename DT1_, typename DT2_, typename DT3_>
        static inline BenchmarkInfo get_benchmark_info(DenseVectorBase<DT1_> & a, DenseVectorBase<DT2_> & b, DT3_ c)
        {
            BenchmarkInfo result;
            result.flops = a.size() * 2;
            result.load = a.size() * (sizeof(DT1_) + sizeof(DT2_)) + sizeof(DT3_);
            result.store = a.size() * sizeof(DT1_);
            result.size.push_back(a.size());
            result.size.push_back(b.size());
            return result; 
        }
    };

    /**
     * \brief Scaled sum of two given vectors and a given scalar.
     *
     * ScaledSum is the class template for the operation
     * \f[
     *     \texttt{ScaledSum}(x, a, y): \quad x \leftarrow a \cdot x + y,
     *     \texttt{ScaledSum}(x, y, b): \quad x \leftarrow x + b \cdot y,
     * \f]
     * which yields the scaled sum of x and y.
     *
     * \todo Implement variant using a.
     *
     * \ingroup grplaoperations
     * \ingroup grplavectoroperations
     */
    template <>
    struct ScaledSum<tags::CPU::SSE>
    {
        /**
         * \name Scaled sums
         * \{
         *
         * \brief Returns the vector x as the scaled sum of two given vectors.
         *
         * \param x The vector that shall not be scaled.
         * \param y The vector that shall be scaled.
         * \param b The scale factor.
         *
         * \retval x Will modify x and return it.
         *
         * \exception VectorSizeDoesNotMatch is thrown if the sizes of x and y do not match.
         */

        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & x, const DenseVectorContinuousBase<float> & y, float b);

        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & x, const DenseVectorContinuousBase<double> & y, double b);

        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b, const DenseVectorContinuousBase<float> & c );

        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b, const DenseVectorContinuousBase<double> & c);

        /// \}
    };

    /**
     * \brief Scaled sum of two given vectors and a given scalar.
     *
     * ScaledSum is the class template for the operation
     * \f[
     *     \texttt{ScaledSum}(x, a, y): \quad x \leftarrow a \cdot x + y,
     *     \texttt{ScaledSum}(x, y, b): \quad x \leftarrow x + b \cdot y,
     * \f]
     * which yields the scaled sum of x and y.
     *
     * \todo Implement variant using a.
     *
     * \ingroup grplaoperations
     * \ingroup grplavectoroperations
     */
    template <>
    struct ScaledSum<tags::GPU::CUDA>
    {
        /**
         * \name Scaled sums
         * \{
         *
         * \brief Returns the vector x as the scaled sum of two given vectors.
         *
         * \param x The vector that shall not be scaled.
         * \param y The vector that shall be scaled.
         * \param b The scale factor.
         *
         * \retval x Will modify x and return it.
         *
         * \exception VectorSizeDoesNotMatch is thrown if the sizes of x and y do not match.
         */

        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & x, const DenseVectorContinuousBase<float> & y, float b);

        /// \}
    };

    /**
     * \brief Scaled sum of two given vectors and a given scalar.
     *
     * ScaledSum is the class template for the operation
     * \f[
     *     \texttt{ScaledSum}(x, a, y): \quad x \leftarrow a \cdot x + y,
     *     \texttt{ScaledSum}(x, y, b): \quad x \leftarrow x + b \cdot y,
     * \f]
     * which yields the scaled sum of x and y.
     *
     * \todo Implement variant using a.
     *
     * \ingroup grplaoperations
     * \ingroup grplavectoroperations
     */
    template <>
    struct ScaledSum<tags::Cell>
    {
        /**
         * \name Scaled sums
         * \{
         *
         * \brief Returns the vector x as the scaled sum of two given vectors.
         *
         * \param x The vector that shall not be scaled.
         * \param y The vector that shall be scaled.
         * \param b The scale factor.
         *
         * \retval x Will modify x and return it.
         *
         * \exception VectorSizeDoesNotMatch is thrown if the sizes of x and y do not match.
         */

        static DenseVectorContinuousBase<float> & value(DenseVectorContinuousBase<float> & x, const DenseVectorContinuousBase<float> & y, const float & b);

        static DenseVectorContinuousBase<double> & value(DenseVectorContinuousBase<double> & x, const DenseVectorContinuousBase<double> & y, const double & b);

        static DenseMatrix<float> & value(DenseMatrix<float> & x, const DenseMatrix<float> & y, const float & b);

        /// \}
    };

    template <> struct ScaledSum<tags::CPU::MultiCore> : public MCScaledSum<tags::CPU::MultiCore> {};

    template <> struct ScaledSum<tags::CPU::MultiCore::SSE> : public MCScaledSum<tags::CPU::MultiCore::SSE> {};

}
#endif
