/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
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

#ifndef LIBLA_GUARD_VECTOR_ELEMENT_SUM_HH
#define LIBLA_GUARD_VECTOR_ELEMENT_SUM_HH 1

#include <libla/tags.hh>
#include <libla/dense_vector.hh>
#include <libla/sparse_vector.hh>

/**
 * \file
 *
 * Templatized definition of summing up all of a Vector's elements.<br/>
 *
 * \ingroup grpvectoroperations
 */
namespace pg512 ///< \todo Namespace name?
{
    /**
     * A VectorElementSum yields the sum of all elements of descendants of type
     * Vector.
     **/
    template <typename DataType_, typename Tag_ = tags::CPU> struct VectorElementSum
    {
        static DataType_ value(const DenseVector<DataType_> & vector)
        {
            DataType_ result(0);

            for (typename Vector<DataType_>::ConstElementIterator i(vector.begin_elements()), i_end(vector.end_elements()) ;
                    i != i_end ; ++i)
            {
                result += *i;
            }

            return result;
        }

        static DataType_ value(const SparseVector<DataType_> & vector)
        {
            DataType_ result(0);

            for (typename Vector<DataType_>::ConstElementIterator i(vector.begin_non_zero_elements()), i_end(vector.end_non_zero_elements()) ;
                    i != i_end ; ++i)
            {
                result += *i;
            }

            return result;
        }
    };
}

#endif
