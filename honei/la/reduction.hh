/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Sven Mallach <sven.mallach@uni-dortmund.de>
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_REDUCTION_HH
#define LIBLA_GUARD_REDUCTION_HH 1

#include <honei/la/banded_matrix.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/reduction-fwd.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/la/sparse_vector.hh>
#include <honei/la/matrix_error.hh>
#include <honei/util/tags.hh>

#include <honei/util/configuration.hh>
#include <honei/util/partitioner.hh>
#include <honei/util/benchmark_info.hh>

namespace honei
{
    template <ReductionType type_, typename Tag_ = tags::CPU> struct Reduction;

    /**
     * \brief Reduction of an entity to the sum of its elements.
     *
     * Reduction is the class template for all types of reductions
     * to the sum of an entity's elements. Every Reduction means
     * the loss of one dimension to the entity, i.e. a Matrix is reduced
     * to a Vector, and a Vector is reduced to a scalar.
     *
     * \f[
     *     \texttt{Reduction}(a): \quad x \leftarrow \sum(a_0, \dots, a_{size - 1}),
     * \f]
     *
     * which yields the entity x, the reduction.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <> struct Reduction<rt_sum>
    {
        /**
         * \name Reductions
         * \{
         *
         * \brief Returns the scalar which is the sum of
         * \brief all elements of a given entity.
         *
         * \param a The entity to be reduced.
         *
         * \retval x If a is a matrix, a vector will be returned, else it will be a scalar.
         *
         */

        template <typename DT_>
        static DenseVector<DT_> value(const DenseMatrix<DT_> & a)
        {
            CONTEXT("When reducing DenseMatrix to DenseVector by sum:");
            DenseVector<DT_> result(a.rows());

            for (unsigned long i(0) ; i < a.rows() ; ++i) /// \todo BandIterator!
            {
                result[i] = Reduction<rt_sum>::value(a[i]);
            }

            return result;
        }

        template <typename DT_>
        static DenseVector<DT_> value(const SparseMatrix<DT_> & a)
        {
            CONTEXT("When reducing SparseMatrix to DenseVector by sum:");
            DenseVector<DT_> result(a.rows());

            for (unsigned long i(0) ; i < a.rows() ; ++i) /// \todo BandIterator!
            {
                /// \ todo use optimized Reduction(SparseVector) instead of Reduction(Vector)
                result[i] = Reduction<rt_sum>::value(a[i]);
            }

            return result;
        }

        template <typename DT_>
        static DenseVector<DT_> value(const BandedMatrix<DT_> & a)
        {
            CONTEXT("When reducing BandedMatrix to DenseVector by sum:");
            DenseVector<DT_> result(a.rows(), DT_(0));

            for (typename BandedMatrix<DT_>::ConstElementIterator i(a.begin_elements()),
                    i_end(a.end_elements()) ; i != i_end ; ++i)
            {
                result[i.row()] += *i;
            }

            return result;
        }

        template <typename DT_>
        static DT_ value(const DenseVectorBase<DT_> & vector)
        {
            CONTEXT("When reducing Vector to Scalar by sum:");

            DT_ result(0);

            for (typename DenseVectorBase<DT_>::ConstElementIterator i(vector.begin_elements()), i_end(vector.end_elements()) ;
                    i != i_end ; ++i)
            {
                result += *i;
            }

            return result;
        }

        template <typename DT_>
        static DT_ value(const SparseVector<DT_> & vector)
        {
            CONTEXT("When reducing SparseVector to Scalar by sum:");

            DT_ result(0);

            for (typename SparseVector<DT_>::NonZeroConstElementIterator i(vector.begin_non_zero_elements()), i_end(vector.end_non_zero_elements()) ;
                    i != i_end ; ++i)
            {
                result += *i;
            }

            return result;
        }

        /// \}

        template <typename DT1_>
        static inline BenchmarkInfo get_benchmark_info(DenseMatrix<DT1_> & a)
        {
            BenchmarkInfo result;
            DenseVector<DT1_> temp(a.columns());
            result = Reduction<rt_sum>::get_benchmark_info(temp) * a.rows();
            result.size.clear();
            result.size.push_back(a.rows() * a.columns());
            return result;
        }

        template <typename DT1_>
        static inline BenchmarkInfo get_benchmark_info(SparseMatrix<DT1_> & a)
        {
            BenchmarkInfo result;
            for (unsigned long i(0) ; i < a.rows() ; ++i)
            {
                result = result + Reduction<rt_sum>::get_benchmark_info(a[i]);
                result.scale += a[i].used_elements();
            }
            result.size.push_back(a.rows() * a.columns());
            result.scale = (double(a.rows() * a.columns()) / result.scale);
            return result;
        }

        template <typename DT1_>
        static inline BenchmarkInfo get_benchmark_info(DenseVectorBase<DT1_> & a)
        {
            BenchmarkInfo result;
            result.flops = a.size();
            result.load = sizeof(DT1_) * a.size();
            result.store = sizeof(DT1_);
            result.size.push_back(a.size());
            return result;
        }

        template <typename DT1_>
        static inline BenchmarkInfo get_benchmark_info(SparseVector<DT1_> & a)
        {
            BenchmarkInfo result;
            result.flops = a.used_elements();
            result.load = sizeof(DT1_) * a.used_elements();
            result.store = sizeof(DT1_);
            result.size.push_back(a.size());
            return result;
        }
    };

    /**
     * \brief Reduction of an entity to the maximum of its elements.
     *
     * Reduction is the class template for all types of reductions
     * to the maximum of an entity's elements. Every Reduction means
     * the loss of one dimension to the entity, i.e. a Matrix is reduced
     * to a Vector, and a Vector is reduced to a scalar.
     *
     * \f[
     *     \texttt{Reduction}(a): \quad x \leftarrow \max(a_0, \dots, a_{size - 1}),
     * \f]
     *
     * which yields the entity x, the reduction.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <> struct Reduction<rt_max>
    {
        /**
         * \brief Returns the scalar which is the maximum of
         * \brief all elements of a given entity.
         *
         * \param a The entity to be reduced.
         *
         * \retval x If a is a matrix, a vector will be returned, else it will be a scalar.
         */

        /// \{
        template <typename DT_>
        static DenseVector<DT_> value(const DenseMatrix<DT_> & matrix)
        {
            CONTEXT("When reducing DenseMatrix to Vector by max:");

            DenseVector<DT_> result(matrix.rows());

            for (unsigned int i(0); i < matrix.rows() ; ++i)
            {
                    result[i] = Reduction<rt_max>::value(matrix[i]);
            }

            return result;
        }

                template <typename DT_>
        static DenseVector<DT_> value(const SparseMatrix<DT_> & matrix)
        {
            CONTEXT("When reducing SparseMatrix to Vector by max:");

            DenseVector<DT_> result(matrix.rows());

            for (unsigned int i(0); i < matrix.rows() ; ++i)
            {
                    result[i] = Reduction<rt_max>::value(matrix[i]);
            }

            return result;
        }

                template <typename DT_>
        static DenseVector<DT_> value(const BandedMatrix<DT_> & matrix)
        {
            CONTEXT("When reducing BandedMatrix to Vector by max:");
            /// \todo Use band interator.
            DenseVector<DT_> result(matrix.rows());

            for (typename BandedMatrix<DT_>::ConstElementIterator i(matrix.begin_elements()),
                    i_end(matrix.end_elements()) ; i != i_end ; ++i)
            {
                if (i.column()==0)
                    result[i.row()] = *i;
                else if (*i > result[i.row()])
                    result[i.row()] = *i;
            }
            return result;
        }

                template <typename DT_>
        static DT_ value(const DenseVectorBase<DT_> & vector)
        {
            CONTEXT("When reducing DenseVectorBase to Scalar by max:");

            DT_ result(vector[0]);

            for (typename DenseVectorBase<DT_>::ConstElementIterator l(vector.begin_elements()),
                    l_end(vector.end_elements()) ; l != l_end ; ++l)
            {
                if (*l > result)
                {
                    result = *l;
                }
            }

            return result;
        }

        template <typename DT_>
        static DT_ value(const SparseVector<DT_> & vector)
        {
            CONTEXT("When reducing SparseVector to Scalar by max:");

            DT_ result(vector[0]);

            for (typename SparseVector<DT_>::ConstElementIterator l(vector.begin_elements()),
                    l_end(vector.end_elements()) ; l != l_end ; ++l)
            {
                if (*l > result)
                {
                    result = *l;
                }
            }

            return result;
        }

        /// \}
    };

    /**
     * \brief Reduction of an entity to the minimum of its elements.
     *
     * Reduction is the class template for all types of reductions
     * to the minimum of an entity's elements. Every Reduction means
     * the loss of one dimension to the entity, i.e. a Matrix is reduced
     * to a Vector, and a Vector is reduced to a scalar.
     *
     * \f[
     *     \texttt{Reduction}(a): \quad x \leftarrow \min(a_0, \dots, a_{size - 1}),
     * \f]
     *
     * which yields the entity x, the reduction.
     *
     * \ingroup grplaoperations
     * \ingroup grplamatrixoperations
     * \ingroup grplavectoroperations
     */
    template <> struct Reduction<rt_min>
    {
        /**
         * \name Reductions
         * \{
         *
         * \brief Returns the scalar which is the minimum of
         * \brief all elements of a given entity.
         *
         * \param a The entity to be reduced.
         *
         * \retval x Will return scalar x.
         */

        template <typename DT_>
        static DenseVector<DT_> value(const DenseMatrix<DT_> & matrix)
        {
            CONTEXT("When reducing DenseMatrix to Vector by min:");

            DenseVector<DT_> result(matrix.rows());

            for (unsigned int i(0); i < matrix.rows() ; ++i)
            {
                    result[i] = Reduction<rt_min>::value(matrix[i]);
            }

            return result;
        }

        template <typename DT_>
        static DenseVector<DT_> value(const SparseMatrix<DT_> & matrix)
        {
            CONTEXT("When reducing SparseMatrix to Vector by min:");

            DenseVector<DT_> result(matrix.rows());

            for (unsigned int i(0); i < matrix.rows() ; ++i)
            {
                    result[i] = Reduction<rt_min>::value(matrix[i]);
            }

            return result;
        }

        template <typename DT_>
        static DenseVector<DT_> value(const BandedMatrix<DT_> & matrix)
        {
            CONTEXT("When reducing BandedMatrix to Vector by max:");
            /// \todo Use band interator.

            DenseVector<DT_> result(matrix.rows());

            for (typename BandedMatrix<DT_>::ConstElementIterator i(matrix.begin_elements()),
                    i_end(matrix.end_elements()) ; i != i_end ; ++i)
            {
                if (i.column() == 0)
                {
                    result[i.row()] = *i;
                }
                else if (*i < result[i.row()])
                {
                    result[i.row()] = *i;
                }
            }
            return result;
        }

        template <typename DT_>
        static DT_ value(const DenseVectorBase<DT_> & vector)
        {
            CONTEXT("When reducing DenseVector to Scalar by min:");

            DT_ result(vector[0]);

            for (typename DenseVectorBase<DT_>::ConstElementIterator l(vector.begin_elements()),
                    l_end(vector.end_elements()) ; l != l_end ; ++l)
            {
                if (*l < result)
                {
                    result = *l;
                }
            }

            return result;
        }

        template <typename DT_>
        static DT_ value(const SparseVector<DT_> & vector)
        {
            CONTEXT("When reducing SparseVector to Scalar by min:");

            DT_ result(vector[0]);

            for (typename SparseVector<DT_>::ConstElementIterator l(vector.begin_elements()),
                    l_end(vector.end_elements()) ; l != l_end ; ++l)
            {
                if (*l < result)
                {
                    result = *l;
                }
            }

            return result;
        }

        /// \}
    };

    /**
     * \brief reduction of a vector to scalar (sum)
     *
     * Reduction is the class template for the operation
     * \f[
     *     \texttt{Reduction}(x): \quad r \leftarrow \sum x_i,
     * \f]
     * which yields the reduction of the given vector x .
     *
     * \ingroup grplaoperations
     * \ingroup grplavectoroperations
     */
    template <> struct Reduction<rt_sum,tags::CPU::SSE>
    {
        /**
         * \name Reduction
         * \{
         *
         * Returns the sum-reduction of a given vector.
         *
         * \param x One vectors of which reduction shall be computed.
         * 
         *
         * \retval r Will return an scalar instance of the used data type.
         *
         * 
         */

        static float value(const DenseVectorContinuousBase<float> & a);

        static double value(const DenseVectorContinuousBase<double> & a);

        static DenseVector<float> value(const DenseMatrix<float> & a);

        static DenseVector<double> value(const DenseMatrix<double> & a);

        static float value(const SparseVector<float> & a);

        static double value(const SparseVector<double> & a);

        static DenseVector<float> value(const SparseMatrix<float> & a);

        static DenseVector<double> value(const SparseMatrix<double> & a);

    };
    /**
     * \brief reduction of a vector to scalar (sum)
     *
     * Reduction is the class template for the operation
     * \f[
     *     \texttt{Reduction}(x): \quad r \leftarrow \sum x_i,
     * \f]
     * which yields the reduction of the given vector x .
     *
     * \ingroup grplaoperations
     * \ingroup grplavectoroperations
     */
    template <> struct Reduction<rt_sum,tags::Cell>
    {
        /**
         * \name Reduction
         * \{
         *
         * Returns the sum-reduction of a given vector.
         *
         * \param x One vectors of which reduction shall be computed.
         * 
         *
         * \retval r Will return an scalar instance of the used data type.
         *
         * 
         */

        static float value(const DenseVectorContinuousBase<float> & a);
        static DenseVector<float> value(const DenseMatrix<float> & a);
        static float value(const SparseVector<float> & a);
        static DenseVector<float> value(const SparseMatrix<float> & a);

        static double value(const DenseVectorContinuousBase<double> & a);
        static DenseVector<double> value(const DenseMatrix<double> & a);
        static double value(const SparseVector<double> & a);
        static DenseVector<double> value(const SparseMatrix<double> & a);


        /// \}
    };

        /**
        * \brief reduction of a vector to scalar (min)
        *
        * Reduction is the class template for the operation
        * \f[
        *     \texttt{Reduction}(x): \quad r \leftarrow \min x_i,
        * \f]
        * which yields the reduction of the given vector x .
        *
        * \ingroup grplaoperations
        * \ingroup grplavectoroperations
        */
    template <> struct Reduction<rt_min,tags::Cell>
    {
        /**
         * \name Reduction
         * \{
         *
         * Returns the min-reduction of a given vector.
         *
         * \param x One vectors of which reduction shall be computed.
         * 
         *
         * \retval r Will return an scalar instance of the used data type.
         *
         * 
         */

        static float value(const DenseVectorContinuousBase<float> & a);
        static DenseVector<float> value(const DenseMatrix<float> & a);
        static float value(const SparseVector<float> & a);
        static DenseVector<float> value(const SparseMatrix<float> & a);

        /// \}
    };

        /**
        * \brief reduction of a vector to scalar (max)
        *
        * Reduction is the class template for the operation
        * \f[
        *     \texttt{Reduction}(x): \quad r \leftarrow \max x_i,
        * \f]
        * which yields the reduction of the given vector x .
        *
        * \ingroup grplaoperations
        * \ingroup grplavectoroperations
        */
    template <> struct Reduction<rt_max,tags::Cell>
    {
        /**
         * \name Reduction
         * \{
         *
         * Returns the max-reduction of a given vector.
         *
         * \param x One vectors of which reduction shall be computed.
         * 
         *
         * \retval r Will return an scalar instance of the used data type.
         *
         * 
         */

        static float value(const DenseVectorContinuousBase<float> & a);
        static DenseVector<float> value(const DenseMatrix<float> & a);
        static float value(const SparseVector<float> & a);
        static DenseVector<float> value(const SparseMatrix<float> & a);

        /// \}
    };

    template <ReductionType type_> struct Reduction<type_, tags::CPU::MultiCore> : public Reduction<type_, tags::CPU>{};
    template <ReductionType type_> struct Reduction<type_, tags::CPU::MultiCore::SSE> : public Reduction<type_, tags::CPU::SSE>{};
}
#endif
