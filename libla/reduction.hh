/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Sven Mallach <sven.mallach@uni-dortmund.de>
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#include <libla/banded_matrix.hh>
#include <libla/dense_matrix.hh>
#include <libla/dense_vector.hh>
#include <libla/sparse_matrix.hh>
#include <libla/sparse_vector.hh>
#include <libla/matrix_error.hh>
#include <libutil/tags.hh>

namespace honei
{
    /**
     * ReductionType is a template tag parameter for the Reduction class
     * template. It governs the type of the (mathematical) reduction that shall
     * be computed.
     *
     * \ingroup grplaoperations
     */
    enum ReductionType
    {
        rt_sum = 0,
        rt_max,
        rt_min
    };

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

        template <typename DT1_>
        static inline BenchmarkInfo get_benchmark_info(unsigned long rows, unsigned long columns = 1, double nonzero = 1)
        {
            BenchmarkInfo result;
            result.flops = 0;
            result.load = 0;
            result.store = 0;
            cout << endl << "!! No detailed benchmark info available !!" << endl;
        }

        template <typename DT_>
        static DenseVector<DT_> value(const DenseMatrix<DT_> & a)
        {
            CONTEXT("When reducing DenseMatrix to DenseVector by sum:");
            DenseVector<DT_> result(a.rows());

            for (unsigned long i(0) ; i < a.rows() ; ++i) /// \todo VectorIterator!
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

            for (unsigned long i(0) ; i < a.rows() ; ++i) /// \todo VectorIterator!
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

            for (typename Matrix<DT_>::ConstElementIterator i(a.begin_elements()),
                    i_end(a.end_elements()) ; i != i_end ; ++i)
            {
                result[i.row()] += *i;
            }

            return result;
        }

        template <typename DT_>
        static DT_ value(const Vector<DT_> & vector)
        {
            CONTEXT("When reducing Vector to Scalar by sum:");

            DT_ result(0);

            for (typename Vector<DT_>::ConstElementIterator i(vector.begin_elements()), i_end(vector.end_elements()) ;
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

            for (typename Vector<DT_>::ConstElementIterator i(vector.begin_non_zero_elements()), i_end(vector.end_non_zero_elements()) ;
                    i != i_end ; ++i)
            {
                result += *i;
            }

            return result;
        }

        /// \}
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
            CONTEXT("When reducing DenseMatrix to Vector by max");

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
            CONTEXT("When reducing SparseMatrix to Vector by max");

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
            CONTEXT("When reducing BandedMatrix to Vector by max");
            /// \todo Use band interator.
            DenseVector<DT_> result(matrix.rows());

            for (typename Matrix<DT_>::ConstElementIterator i(matrix.begin_elements()),
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
        static DT_ value(const DenseVector<DT_> & vector)
        {
            CONTEXT("When reducing DenseVector to Scalar by max");

            DT_ result(vector[0]);

            for (typename Vector<DT_>::ConstElementIterator l(vector.begin_elements()),
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
            CONTEXT("When reducing SparseVector to Scalar by max");

            DT_ result(vector[0]);

            for (typename Vector<DT_>::ConstElementIterator l(vector.begin_elements()),
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
            CONTEXT("When reducing DenseMatrix to Vector by min");

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
            CONTEXT("When reducing SparseMatrix to Vector by min");

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
            CONTEXT("When reducing BandedMatrix to Vector by max");
            /// \todo Use band interator.

            DenseVector<DT_> result(matrix.rows());

            for (typename Matrix<DT_>::ConstElementIterator i(matrix.begin_elements()),
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
        static DT_ value(const DenseVector<DT_> & vector)
        {
            CONTEXT("When reducing DenseVector to Scalar by min");

            DT_ result(vector[0]);

            for (typename Vector<DT_>::ConstElementIterator l(vector.begin_elements()),
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
            CONTEXT("When reducing SparseVector to Scalar by min");

            DT_ result(vector[0]);

            for (typename Vector<DT_>::ConstElementIterator l(vector.begin_elements()),
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

        static float value(const DenseVector<float> & a);

        /// \}
    };
}
#endif
