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

#ifndef LIBLA_GUARD_VECTOR_ELEMENT_INVERSE_HH
#define LIBLA_GUARD_VECTOR_ELEMENT_INVERSE_HH 1

#include <libutil/tags.hh>
#include <libla/dense_vector.hh>
#include <libla/sparse_vector.hh>
#include <libla/vector_error.hh>

/**
 * \file
 *
 * Templatized definitions of vector inversion.
 *
 * \ingroup grpvectoroperations
 **/
namespace pg512
{

    /**
     * VectorElementInverse is the class template for the elementwise inversion of a vector.
     * \brief The referenced vector is changed under this operation.
     * \ingroup grpvectoroperations
     **/
    template <typename Tag_ = tags::CPU> struct VectorElementInverse
    {
        /**
         * Returns the inverse vector of a referenced DenseVector .
         *
         * \param vector Reference to a dense vector that will be inverted and used as return argument.
         **/
        template <typename DataType_> static DenseVector<DataType_> & value(DenseVector<DataType_> & vector)
        {

            for (typename Vector<DataType_>::ElementIterator l(vector.begin_elements()),
                    l_end(vector.end_elements()) ; l != l_end ; ++l)
            {
                 if (*l == static_cast<DataType_>(0))
                    continue;

                *l = DataType_(1) / *l;
            }
            return vector;
        }

        /**
         * Returns the inverse vector of a referenced SparseVector .
         *
         * \param vector Reference to a sparse vector that will be inverted and used as return argument.
         **/
        template <typename DataType_> static SparseVector<DataType_> & value(SparseVector<DataType_> & vector)
        {

            for (typename Vector<DataType_>::ElementIterator l(vector.begin_non_zero_elements()),
                    l_end(vector.end_non_zero_elements()) ; l != l_end ; ++l)
            {
                *l = DataType_(1) / *l;
            }
            return vector;
        }

    };
}
#endif
