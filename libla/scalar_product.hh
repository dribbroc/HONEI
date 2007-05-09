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

#ifndef LIBLA_GUARD_SCALAR_PRODUCT_HH
#define LIBLA_GUARD_SCALAR_PRODUCT_HH 1

#include "tags.hh"
#include "vector.hh"

namespace pg512 ///< \todo Namespace name?
{
    /**
     * A ScalarProduct yields the inner product of two descendants of type Vector.
     **/
    template <typename DataType_, typename Tag_ = tags::CPU> struct ScalarProduct
    {
        static DataType_ value(const Vector<DataType_> & left, const Vector<DataType_> & right)
        {
            if (left.size() != right.size())
                throw std::string("Yikes. VectorSizesDoNotMatch and no exception classes yet.");

            DataType_ result(0);

            for (typename Vector<DataType_>::ElementIterator l(left.begin_elements()), l_end(left.end_elements()),
                    r(right.begin_elements()) ; l != l_end ; ++l )
            {
                result += (*l) * (*r);
                ++r;
            }

            return result;
        }

        /// \todo Inner product of SparseVector/Vector and
        /// SparseVector/SparseVector.
    };
}

#endif
