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

#ifndef LIBLA_GUARD_VECTOR_MASK_HH
#define LIBLA_GUARD_VECTOR_MASK_HH 1

#include <libutil/tags.hh>
#include <libla/dense_vector.hh>
#include <libla/sparse_vector.hh>
#include <libla/vector.hh>
#include <libla/vector_error.hh>

/**
 * \file
 *
 * Implementation of VectorMask.
 *
 * \ingroup grpvectoroperations
 **/
namespace pg512
{
    /**
     * \brief VectorMask masks a given vector by setting all elements to zero, for which the mask vector contains the value "false".
     * \brief Member functions use first parameter (reference) as return argument.
     *
     * \ingroup grpvectoroperations
     **/
    template <typename DataType_, typename Tag_ = tags::CPU> struct VectorMask
    {
         /**
         * Return a masked DenseVector.
         *
         * \param vector The DenseVector to be masked.
         * \param mask DenseVector<bool> to be used for masking.
         **/

        static DenseVector<DataType_> & value(DenseVector<DataType_> & vector, const Vector<bool> & mask)
        {
            if (vector.size() != mask.size())
                throw VectorSizeDoesNotMatch(mask.size(), vector.size());

            Vector<bool>::ConstElementIterator r(mask.begin_elements());
            for (typename Vector<DataType_>::ElementIterator l(vector.begin_elements()), l_end(vector.end_elements()) ; l != l_end ; ++l)
            {
                *l = (*r ? *l : static_cast<DataType_>(0));
                ++r;
            }

            return vector;
        }

        /**
         * Return a masked SparseVector. All elements that equal zero will
         * be invariant under this operation.
         *
         * \param vector The SparseVector to be masked.
         * \param mask SparseVector<bool> which is used as mask.
         **/
        static SparseVector<DataType_> & value(SparseVector<DataType_> & vector, const Vector<bool> & mask)
        {
            if (vector.size() != mask.size())
                throw VectorSizeDoesNotMatch(mask.size(), vector.size());

            Vector<bool>::ConstElementIterator r(mask.begin_elements());
            for (typename Vector<DataType_>::ElementIterator l(vector.begin_non_zero_elements()), l_end(vector.end_non_zero_elements()) ; l != l_end ; ++l)
            {
                *l = (*r ? *l : static_cast<DataType_>(0));
                ++r;
            }

            return vector;
        }

    };
}
#endif
