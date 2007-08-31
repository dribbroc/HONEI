/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_DOT_PRODUCT_HH
#define LIBLA_GUARD_DOT_PRODUCT_HH 1

#include <libutil/tags.hh>
#include <libla/dense_vector.hh>
#include <libla/sparse_vector.hh>
#include <libla/vector_error.hh>

/**
 * \file
 *
 * Templatized definitions of operation DotProduct.
 *
 * \ingroup grpoperations
 */
namespace honei
{
    /**
     * DotProduct is the class template for the operation
     * \f[
     *     DotProduct(x, y): \quad r \leftarrow x \cdot y,
     * \f]
     * which yields the scalar or inner product of the given vectors x and y.
     *
     * \ingroup grpoperations
     * \ingroup grpvectoroperations
     * \ingroup grpreductions
     */
    template <typename Tag_ = tags::CPU> struct DotProduct
    {
        /**
         * Returns the dot- (or inner) product of two given vectors.
         *
         * \param x One of the vectors of which the scalar product shall be computed.
         * \param y idem
         *
         * \retval r Will return an instance of the used data type containing the scalar product.
         *
         * \exception VectorSizeDoesNotMatch is thrown if the two vectors don't have the same size.
         */

        /// \{
        template <typename DT1_, typename DT2_>
        static DT1_ value(const Vector<DT1_> & left, const Vector<DT2_> & right)
        {
            CONTEXT("When calculating Vector-Vector dot product:");

            if (left.size() != right.size())
                throw VectorSizeDoesNotMatch(right.size(), left.size());

            DT1_ result(0);

            for (typename Vector<DT1_>::ConstElementIterator l(left.begin_elements()),
                    l_end(left.end_elements()) ; l != l_end ; ++l )
            {
                result += (*l) * right[l.index()];
            }

            return result;
        }

        template <typename DT1_, typename DT2_>
        static DT1_ value(const SparseVector<DT1_> & left, const Vector<DT2_> & right)
        {
            CONTEXT("When calculating SparseVector-Vector dot product:");

            if (left.size() != right.size())
                throw VectorSizeDoesNotMatch(right.size(), left.size());

            DT1_ result(0);

            for (typename Vector<DT1_>::ConstElementIterator l(left.begin_non_zero_elements()),
                    l_end(left.end_non_zero_elements()) ; l != l_end ; ++l )
            {
                result += (*l) * right[l.index()];
            }

            return result;
        }

        template <typename DT1_, typename DT2_>
        static DT1_ value(const SparseVector<DT1_> & left, const SparseVector<DT2_> & right)
        {
            CONTEXT("When calculating SparseVector-SparseVector dot product:");

            if (left.size() != right.size())
                throw VectorSizeDoesNotMatch(right.size(), left.size());

            DT1_ result(0);

            typename Vector<DT2_>::ConstElementIterator r(right.begin_non_zero_elements());
            for (typename Vector<DT1_>::ConstElementIterator l(left.begin_non_zero_elements()),
                    l_end(left.end_non_zero_elements()) ; l != l_end ; )
            {
                if (l.index() == r.index())
                {
                    result += (*l) * (*r);
                    ++l; ++r;
                }
                else if (l.index() < r.index())
                {
                    ++l;
                }
                else
                {
                    ++r;
                }
            }

            return result;
        }

        template <typename IT1_, typename IT2_>
        static typename IT1_::value_type value(IT1_ & left, const IT1_ & left_end,
                IT2_ & right, const IT2_ & right_end)
        {
            CONTEXT("When calculating iterator-based dot product:");

            typename IT1_::value_type result(0);

            for ( ; (left != left_end) && (right != right_end) ; ++left, ++right)
            {
                result += (*left) * (*right);
            }

            return result;
        }
        /// \}
    };
}

#endif
