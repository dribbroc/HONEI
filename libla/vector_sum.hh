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

#ifndef LIBLA_GUARD_VECTOR_SUM_HH
#define LIBLA_GUARD_VECTOR_SUM_HH 1

#include <libutil/tags.hh>
#include <libla/dense_vector.hh>
#include <libla/sparse_vector.hh>
#include <libla/vector_error.hh>

/**
 * \file
 *
 * Templatized definitions of operation VectorSum.
 *
 * \ingroup grpoperations
 */
namespace pg512 ///< \todo Namespace name?
{
    /**
     * \brief Sum of two vectors.
     *
     * VectorSum is the class template for the addition operation
     * \f[
     *     VectorSum(a, b): \quad r \leftarrow a + b,
     * \f]
     * which yields r, the sum of vectors a and b.
     *
     * The return value is the summand a after modification.
     *
     * \ingroup grpoperations
     */
    template <typename Tag_ = tags::CPU> struct VectorSum
    {
        /**
         * \brief Returns the the sum of two given vectors.
         *
         * \param a The vector that is the left-hand summand of the operation.
         * \param b The vector that is the right-hand summand of the operation.
         *
         * \retval Will modify the summand a and return it.
         *
         * \exception VectorSizeDoesNotMatch is thrown if the two vectors don't have the same size.
         */

        /// \{
        template <typename DataType1_, typename DataType2_> static DenseVector<DataType1_> & value(DenseVector<DataType1_> & left, const DenseVector<DataType2_> & right)
        {
            if (left.size() != right.size())
                throw VectorSizeDoesNotMatch(right.size(), left.size());

			typename Vector<DataType2_>::ConstElementIterator r(right.begin_elements());
            for (typename Vector<DataType1_>::ElementIterator l(left.begin_elements()),
                    l_end(left.end_elements()) ; l != l_end ; ++l)
            {
                *l += *r;
				++r;
            }

            return left;
        }

        /**
         * Returns the resulting vector of the sum of two given SparseVector instances.
         *
         * \param left Reference to sparse vector that will be also used as result vector.
         * \param right Reference to constant sparse vector to be added.
         **/
        template <typename DataType1_, typename DataType2_> static SparseVector<DataType1_> & value(SparseVector<DataType1_> & left, const SparseVector<DataType2_> & right)
        {
            if (left.size() != right.size())
                throw VectorSizeDoesNotMatch(right.size(), left.size());

            typename Vector<DataType1_>::ElementIterator l(left.begin_non_zero_elements());
            for (typename Vector<DataType2_>::ConstElementIterator r(right.begin_non_zero_elements()),
                    r_end(right.end_non_zero_elements()) ; r != r_end ; )
            {
                if (r.index() < l.index())
                {
                    left[r.index()] = *r;
                    ++r;
                }
                else if (l.index() < r.index())
                {
                    ++l;
                }
                else
                {
                    *l += *r;
                    ++l; ++r;
                }
            }
			///\todo: perhaps sparsify - i.e. addition of -7 and 7 possible
            return left;
        }

        /**
         * Returns the the resulting vector of the sum of a given dense and a given sparse vector.
         *
         * \param left Reference to dense vector that will be also used as result vector.
         * \param right Reference to constant sparse vector to be added.
         **/
        template <typename DataType1_, typename DataType2_> static DenseVector<DataType1_> & value(DenseVector<DataType1_> & left, const SparseVector<DataType2_> & right)
        {
            if (left.size() != right.size())
                throw VectorSizeDoesNotMatch(right.size(), left.size());

            for (typename Vector<DataType2_>::ConstElementIterator r(right.begin_non_zero_elements()),
                    r_end(right.end_non_zero_elements()) ; r != r_end ; ++r )
            {
                left[r.index()] += *r;
            }

            return left;
        }
        /// \}
    };
}
#endif
