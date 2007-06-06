/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
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

#ifndef LIBLA_GUARD_VECTOR_DIFFERENCE_HH
#define LIBLA_GUARD_VECTOR_DIFFERENCE_HH 1

#include <libla/tags.hh>
#include <libla/dense_vector.hh>
#include <libla/sparse_vector.hh>
#include <libla/vector_error.hh>

/**
 * \file
 *
 * Templatized definitions of vector differences.<br/>
 *
 * \ingroup grpvectoroperations
 **/
namespace pg512
{

    /**
     * VectorSum is the class template for the difference of two vectors.
     *
     * \ingroup grpvectoroperations
     **/
    template <typename DataType_, typename Tag_ = tags::CPU> struct VectorDifference
    {
        /**
         * Returns the the resulting vector of the difference of two given DenseVector instances.
         **/
        static DenseVector<DataType_> value(const DenseVector<DataType_> & left, const DenseVector<DataType_> & right)
        {
            if (left.size() != right.size())
                throw VectorSizeDoesNotMatch(right.size(), left.size());

            DenseVector<DataType_> result(left.size(), 0, 0, 1);

            for (typename Vector<DataType_>::ConstElementIterator l(left.begin_elements()), l_end(left.end_elements()) ; l != l_end ; ++l)
            {
                result[l.index()] = *l - right[l.index()];
            }

            return result;
        }

        /**
         * Returns the resulting vector of the difference of two given SparseVector instances.
         **/
        static SparseVector<DataType_> value(const SparseVector<DataType_> & left, const SparseVector<DataType_> & right)
        {
            if (left.size() != right.size())
                throw VectorSizeDoesNotMatch(right.size(), left.size());

            SparseVector<DataType_> result(left.size(), right.used_elements() + left.used_elements());

            for (typename Vector<DataType_>::ConstElementIterator l(left.begin_non_zero_elements()),
                l_end(left.end_non_zero_elements()), r(right.begin_non_zero_elements()), r_end(right.end_non_zero_elements()) ;
                l != l_end ; ++l)
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
                    result[l.index()] = *l - *r;
                    ++l; ++r;
                }
            }
            return result;
        }

        /**
         * Returns the the resulting vector of the difference of a given dense and a given sparse vector.
         *
         * \param left The dense vector.
         * \param right The sparse vector.
         **/
         static DenseVector<DataType_> value(const DenseVector<DataType_> & left, const SparseVector<DataType_> & right)
        {
            if (left.size() != right.size())
                throw VectorSizeDoesNotMatch(right.size(), left.size());

            DenseVector<DataType_> result(left.size(),0, 0, 1);

            for (typename Vector<DataType_>::ConstElementIterator l(left.begin_elements()),
                l_end(left.end_elements()), r(right.begin_non_zero_elements()), r_end(right.end_non_zero_elements()) ;
                    r != r_end ; )
            {
                while (l.index() < r.index() && (l != l_end))
                {
                    result[l.index()] = *l;
                    ++l;
                }
                result[l.index()] = *l - *r;
                ++r;
            }
            return result;
        }

    };
}
#endif
