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
 * Templatized definitions of vector sums.<br/>
 *
 * \ingroup grpvectoroperations
 **/
namespace pg512
{

    /**
     * VectorSum is the class template for the sum of two vectors.
     * \brief The first referenced vector is changed under this operation.
     * \ingroup grpvectoroperations
     **/
    template <typename Tag_ = tags::CPU> struct VectorSum
    {
        /**
         * Returns the the resulting vector of the sum of two given DenseVector instances.
         *
         * \param left Reference to dense vector that will be also used as result vector.
         * \param right Reference to constant dense vector to be added.
         **/
        template <typename DataType1_, typename DataType2_> static DenseVector<DataType1_> & value(DenseVector<DataType1_> & left, const DenseVector<DataType2_> & right)
        {
            if (left.size() != right.size())
                throw VectorSizeDoesNotMatch(right.size(), left.size());

            for (typename Vector<DataType1_>::ElementIterator l(left.begin_elements()),
                    l_end(left.end_elements()) ; l != l_end ; ++l)
            {
                *l += right[l.index()];
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

            SparseVector<DataType1_> result(left.size(), right.used_elements() + left.used_elements());

            typename Vector<DataType2_>::ConstElementIterator r(right.begin_non_zero_elements()),
                    r_end(right.end_non_zero_elements());
            for (typename Vector<DataType1_>::ElementIterator l(left.begin_non_zero_elements()),
                    l_end(left.end_non_zero_elements()) ; l != l_end ;)
            {
                if (r.index() < l.index())
                {
                    result[r.index()] = *r;
                    ++r;
                }
                else if (l.index() < r.index())
                {
                    result[l.index()] = *l;
                    ++l;
                }
                else
                {
                    result[l.index()] = *l + *r;
                    ++l; ++r;
                }
            }

            return result;
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

            typename Vector<DataType1_>::ElementIterator l(left.begin_elements()),
                    l_end(left.end_elements());
            for (typename Vector<DataType2_>::ConstElementIterator r(right.begin_non_zero_elements()),
                    r_end(right.end_non_zero_elements()) ; r != r_end ; )
            {
                while (l.index() < r.index() && (l != l_end))
                {
                    ++l;
                }

                *l += *r;
                ++r;
            }
            return left;
        }

    };
}
#endif
