/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * This file is part of the Utility C++ library. LibUtil is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibUtil is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#pragma once
#ifndef LIBUTIL_GUARD_SHARED_ARRAY_HH
#define LIBUTIL_GUARD_SHARED_ARRAY_HH 1

#include <honei/util/exception.hh>
#include <honei/util/private_implementation_pattern.hh>

#include <string>
#include <tr1/memory>

namespace honei
{
    /**
     * SharedArrayError is thrown by SharedArray and related classes.
     *
     * \ingroup grpexceptions
     * \ingroup grpsharedarray
     */
    class SharedArrayError :
        public Exception
    {
        public:
            SharedArrayError(const std::string & message) throw();
    };

    /**
     * SharedArray is a class template that allows sharing one array among several owners.
     */
    template <typename DataType_> class SharedArray :
        public PrivateImplementationPattern<SharedArray<DataType_>, Shared>
    {
        public:
            /// \name Constructors and destructor.
            /// \{

            /// (Explicit) constructor.
            explicit SharedArray(unsigned long size);

            /// Copy-constructor.
            SharedArray(const SharedArray<DataType_> & other);

            /// Destructor.
            ~SharedArray();

            /// \}

            /// \name Operators to mimic POA
            /// \{

            /// Subscript operator, return element at given index.
            DataType_ & operator[] (unsigned long index) const;

            /// Return whether we hold any data.
            bool operator! () const;

            /// \}

            /// Return our size.
            unsigned long size() const;

            /// Return our POA.
            DataType_ * get() const;

            /// Reset us with a new size and a new POA.
            void reset(unsigned long size, DataType_ * array) const;
    };

    extern template class SharedArray<float>;

    extern template class SharedArray<double>;

    extern template class SharedArray<int>;

    extern template class SharedArray<unsigned int>;

    extern template class SharedArray<long>;

    extern template class SharedArray<unsigned long>;
}

#endif
