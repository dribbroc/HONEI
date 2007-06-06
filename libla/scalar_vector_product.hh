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

#ifndef LIBLA_GUARD_SCALAR_VECTOR_PRODUCT_HH
#define LIBLA_GUARD_SCALAR_VECTOR_PRODUCT_HH 1

#include <libla/tags.hh>
#include <libla/dense_vector.hh>
#include <libla/sparse_vector.hh>

/**
 * \file
 *
 * Templatized definitions of scalar-vector products.<br/>
 *
 * \ingroup grpvectoroperations
 **/
namespace pg512
{

    /**
     * ScalarVectorProduct is the class template for multiplying a scalar to a vector
     *
     * \ingroup grpvectoroperations
     **/
    template <typename DataType_, typename Tag_ = tags::CPU> struct ScalarVectorProduct
    {
        /**
         * Returns the resulting vector after multiplying a scalar to a given DenseVector instance.
         **/
        static DenseVector<DataType_> value(DataType_ scalar, DenseVector<DataType_> & vector)
        {
            DenseVector<DataType_> result(vector.size(), 0, 0, 1);

            for (typename Vector<DataType_>::ConstElementIterator l(vector.begin_elements()), l_end(vector.end_elements()) ; l != l_end ; ++l)
            {
                result[l.index()] = scalar * *l;
            }

            return result;
        }

        /**
         * Returns the resulting vector after multiplying a scalar to a given SparseVector instance.
         **/
        static SparseVector<DataType_> value(DataType_ scalar, SparseVector<DataType_> & vector)
        {
            SparseVector<DataType_> result(vector.size(), 0, 0, 1);

            for (typename Vector<DataType_>::ConstElementIterator l(vector.begin_non_zero_elements()), l_end(vector.end_non_zero_elements()) ; l != l_end ; ++l)
            {
                result[l.index()] = scalar * *l;
            }

            return result;
        }
    };
}
#endif
