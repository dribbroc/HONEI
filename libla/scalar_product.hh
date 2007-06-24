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

#ifndef LIBLA_GUARD_SCALAR_PRODUCT_HH
#define LIBLA_GUARD_SCALAR_PRODUCT_HH 1

#include <libutil/tags.hh>
#include <libla/dense_vector.hh>
#include <libla/sparse_vector.hh>
#include <libla/vector_error.hh>

/**
 * \file
 *
 * Implementation of ScalarProduct.
 *
 * \ingroup grpvectoroperations
 */
namespace pg512 ///< \todo Namespace name?
{
    /**
     * \brief ScalarProduct yields the inner product of two descendants of type Vector.
     *
     * \ingroup grpvectoroperations
     **/
    template <typename Tag_ = tags::CPU> struct ScalarProduct
    {
        template <typename DataType1_, typename DataType2_> static DataType1_ value(const DenseVector<DataType1_> & left, const DenseVector<DataType2_> & right)
        {
            if (left.size() != right.size())
                throw VectorSizeDoesNotMatch(right.size(), left.size());

            DataType1_ result(0);

            for (typename Vector<DataType1_>::ConstElementIterator l(left.begin_elements()),
                    l_end(left.end_elements()) ; l != l_end ; ++l )
            {
                result += (*l) * right[l.index()];
            }

            return result;
        }

        template <typename DataType1_, typename DataType2_> static DataType1_ value(const SparseVector<DataType1_> & left, const Vector<DataType2_> & right)
        {
            if (left.size() != right.size())
                throw VectorSizeDoesNotMatch(right.size(), left.size());

            DataType1_ result(0);

            for (typename Vector<DataType1_>::ConstElementIterator l(left.begin_non_zero_elements()),
                    l_end(left.end_non_zero_elements()) ; l != l_end ; ++l )
            {
                result += (*l) * right[l.index()];
            }
        }
    };
/*
    /// Explicit instantiation for Cell-based float computation.
    template <>
    float ScalarProduct<tags::Cell>::value(const DenseVector<float> & left, const DenseVector<float> & right);
*/
}

#endif
