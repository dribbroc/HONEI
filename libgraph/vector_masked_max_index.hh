/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Sven Mallach <sven.mallach@uni-dortmund.de>
 *
 * This file is part of the Graph C++ library. LibGraph is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibGraph is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef LIBLA_GUARD_VECTOR_MASKED_MAX_INDEX_HH
#define LIBLA_GUARD_VECTOR_MASKED_MAX_INDEX_HH 1

#include <libutil/tags.hh>
#include <libla/dense_vector.hh>
#include <libla/dense_vector-impl.hh>
#include <libla/sparse_vector.hh>
#include <libla/vector.hh>
#include <libla/vector_error.hh>

/**
 * \file
 *
 * Implementation of VectorMaskedMaxIndex.
 *
 * \ingroup grpvectoroperations
 **/
namespace honei
{
    /**
     * \brief VectorMaskedMaxIndex retrieves the index of the maximum element of a vector.
     * \brief Every Element for which the mask vector has the value "false" will be ignored in search of maximum.
     * \brief All parameters will be invariant under this operation.
     * \ingroup grpvectoroperations
     **/
    template <typename DataType_, typename Tag_ = tags::CPU> struct VectorMaskedMaxIndex
    {
         /**
         * Return the index of the maximum element of a masked vector.
         *
         * \param vector The DenseVector to be used.
         * \param mask Vector<bool> to be used for masking.
         **/

        static unsigned long value(const DenseVector<DataType_> & vector, const Vector<bool> & mask)
        {
            CONTEXT("When calculating the maximum element of a masked DenseVector");
            if (vector.size() != mask.size())
                throw VectorSizeDoesNotMatch(mask.size(), vector.size());

            unsigned long result(-1);
            DataType_ temp(std::numeric_limits<DataType_>::min());
            Vector<bool>::ConstElementIterator r(mask.begin_elements());

            for (typename Vector<DataType_>::ConstElementIterator l(vector.begin_elements()), l_end(vector.end_elements()) ; l != l_end ; ++l)
            {
                if (*r && *l > temp)
                {
                    result = l.index();
                    temp = *l;
                }
                ++r;
            }

            return result;
        }

         /**
         * Return the index of the maximum element of a masked vector.
         *
         * \param vector The SparseVector to be used.
         * \param mask Vector<bool> to be used for masking.
         **/

        static unsigned long value(const SparseVector<DataType_> & vector, const Vector<bool> & mask)
        {
            CONTEXT("When calculating the maximum element of a masked SparseVector");
            if (vector.size() != mask.size())
                throw VectorSizeDoesNotMatch(mask.size(), vector.size());

            unsigned long result(-1);
            DataType_ temp(std::numeric_limits<DataType_>::min());
            Vector<bool>::ConstElementIterator r(mask.begin_elements());

            for (typename Vector<DataType_>::ConstElementIterator l(vector.begin_elements()), l_end(vector.end_elements()) ; l != l_end ; ++l)
            {
                if (*r && *l > temp)
                {
                    result = l.index();
                    temp = *l;
                }
                ++r;
            }

            return result;
        }

    };
}
#endif
