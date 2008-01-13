/* vim: set sw=4 sts=4 et foldmethod=syntax : */

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

#ifndef LIBLA_GUARD_ELEMENT_ITERATOR_HH
#define LIBLA_GUARD_ELEMENT_ITERATOR_HH 1

#include <iterator>
#include <tr1/memory>

namespace honei
{
    template <typename DataType_> class Vector;

    template <typename DataType_> class Matrix;

    /**
     * IteratorTraits declares the interface for the iterators' parent classes.
     */
    template <typename Container_> class IteratorTraits;

    /// Specialisation of IteratorTraits for vector-like types.
    template <typename DataType_> class IteratorTraits<Vector<DataType_> >
    {
        public:
            /// Returns our index.
            virtual unsigned long index() const = 0;

            /// Returns a pointer to our parent container.
            virtual const Vector<DataType_> * parent() const = 0;

            /// Destructor.
            virtual ~IteratorTraits() {}
    };

    /// Specialisation of IteratorTraits for matrix-like types.
    template <typename DataType_> class IteratorTraits<Matrix<DataType_> >
    {
        public:
            /// Returns our index.
            virtual unsigned long index() const = 0;

            /// Returns our column index.
            virtual unsigned long column() const = 0;

            /// Returns our row index.
            virtual unsigned long row() const = 0;

            /// Returns a pointer to our parent container.
            virtual const Matrix<DataType_> * parent() const = 0;

            /// Destructor.
            virtual ~IteratorTraits() {}
    };

    /**
     * IteratorBase declares the minimal interface for the implementation of an ElementIterator.
     */
    template <typename DataType_, typename Container_> class IteratorBase :
        public std::iterator<std::forward_iterator_tag, DataType_>,
        public IteratorTraits<Container_>
    {
        public:
            /// Preincrement operator.
            virtual IteratorBase<DataType_, Container_> & operator++ () = 0;

            /// In-place-add operator.
            virtual IteratorBase<DataType_, Container_> & operator+= (const unsigned long step) = 0;

            /// Dereference operator that returns an assignable reference.
            virtual DataType_ & operator* () = 0;

            /// Dereference operator that returns an unassignable reference.
            virtual const DataType_ & operator* () const = 0;

            /// Comparison operator for less-than.
            virtual bool operator< (const IteratorBase<DataType_, Container_> & other) const = 0;

            /// Comparison operator for equality.
            virtual bool operator== (const IteratorBase<DataType_, Container_> & other) const = 0;

            /// Comparison operator for inequality.
            virtual bool operator!= (const IteratorBase<DataType_, Container_> & other) const = 0;

            /// Destructor.
            virtual ~IteratorBase() {}
    };
}

#endif
