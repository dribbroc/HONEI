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

#ifndef LIBLA_GUARD_SCALAR_VECTOR_SUM_HH
#define LIBLA_GUARD_SCALAR_VECTOR_SUM_HH 1

#include <libutil/tags.hh>
#include <libla/dense_vector.hh>
#include <libla/sparse_vector.hh>

/**
 * \file
 *
 * Templatized definitions of scalar-vector sums.
 *
 * \ingroup grpvectoroperations
 **/
namespace pg512
{

    /**
     * ScalarVectorSum is the class template for adding a scalar to a vector
     * \brief The referenced vector is changed by adding the given scalar to each of its elements.
     * \brief With a negative scalar, ScalarVectorSum can also be used to compute a ScalarVectorDifference.
     * \ingroup grpvectoroperations
     **/
    template <typename DataType_, typename Tag_ = tags::CPU> struct ScalarVectorSum
    {
        /**
         * Returns the resulting vector after multiplying a scalar to a given DenseVector instance.
         * \param vector DenseVector to be used.
         * \param scalar The scalar to be added.
         **/
        static DenseVector<DataType_> & value(const DataType_ scalar, DenseVector<DataType_> & vector)
        {

            for (typename Vector<DataType_>::ElementIterator l(vector.begin_elements()), l_end(vector.end_elements()) ; l != l_end ; ++l)
            {
                *l += scalar;
            }

            return vector;
        }

        /**
         * Returns the resulting vector after adding a scalar to a given SparseVector instance.
         * \param vector SparseVector to be used.
         * \param scalar The scalar to be added.
         **/
        static SparseVector<DataType_> & value(const DataType_ scalar, SparseVector<DataType_> & vector)
        {

            for (typename Vector<DataType_>::ElementIterator l(vector.begin_elements()), l_end(vector.end_elements()) ; l != l_end ; ++l)
            {
                *l += scalar;
            }

            return vector;
        }
    };
}
#endif
