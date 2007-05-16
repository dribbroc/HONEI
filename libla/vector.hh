/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Michael Abshoff <michael.abshoff@fsmath.mathematik.uni-dortmund.de>
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

#ifndef LIBLA_GUARD_VECTOR_HH
#define LIBLA_GUARD_VECTOR_HH 1

#include <libutil/exception.hh>
#include <libutil/shared_array.hh>
#include <libla/element_iterator.hh>

#include <iterator>
#include <ostream>
#include <string.h>

namespace pg512 ///< \todo Namespace name?
{
    /**
     * A VectorError is thrown when Vector-private methods encounter an exception.
     **/
    class VectorError :
        public Exception
    {
        public:
            VectorError(const std::string & message) throw () :
                Exception(message)
            {
            }
    };

    /**
     * A Vector is the abstract baseclass for all vector-like types used.
     **/
    template <typename DataType_> class Vector
    {
        public:
            /// Type of the iterator over our elements.
            typedef ElementIteratorWrapper<Vector<DataType_>, DataType_> ElementIterator;

            /// Returns iterator pointing to the first element of the vector.
            virtual ElementIterator begin_elements() const = 0;

            /// Returns iterator pointing behind the last element of the vector.
            virtual ElementIterator end_elements() const = 0;

            /// Returns our size.
            virtual unsigned long size() const = 0;

            /// Retrieves element by index, zero-based, unassignable.
            virtual const DataType_ & operator[] (unsigned long index) const = 0;

            /// Retrieves element by index, zero-based, assignable
            virtual DataType_ & operator[] (unsigned long index) = 0;
    };

    /// Output our Vector to an ostream.
    template <typename DataType_> std::ostream & operator<< (std::ostream & lhs, const Vector<DataType_> & v)
    {
        lhs << "[ ";
        for (typename Vector<DataType_>::ElementIterator i(v.begin_elements()), i_end(v.end_elements()) ;
                i != i_end ; ++i)
        {
            lhs << *i << " ";
        }
        lhs << "]";

        return lhs;
    }
}

#endif
