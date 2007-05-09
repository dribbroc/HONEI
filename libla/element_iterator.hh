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

#include <tr1/memory>

namespace pg512 ///< \todo Namespace name?
{
    template <typename DataType_> class Vector;

    template <typename Tag_, typename DataType_> class ElementIteratorBase;

    template <typename DataType_> class ElementIteratorBase<Vector<DataType_>, DataType_> :
        public std::iterator<std::forward_iterator_tag, DataType_>
    {
        public:
            /// Preincrement operator.
            virtual ElementIteratorBase<Vector<DataType_>, DataType_> & operator++ () = 0;

            /// Returns our index.
            virtual const unsigned long index() const = 0;
    };

    template <typename Tag_, typename DataType_> class ElementIteratorImplBase :
        public ElementIteratorBase<Tag_, DataType_>
    {
        public:
            /// Equality operator.
            virtual bool operator== (const ElementIteratorImplBase<Tag_, DataType_> & other) const = 0;

            /// Inqquality operator.
            virtual bool operator!= (const ElementIteratorImplBase<Tag_, DataType_> & other) const = 0;

            /// Dereference operator.
            virtual DataType_ & operator* () const = 0;

            /// Returns pointer to our Vector.
            virtual const Vector<DataType_> * vector() const = 0;
    };

    template <typename Tag_, typename DataType_> class ElementIteratorWrapper;

    template <typename DataType_> class ElementIteratorWrapper<Vector<DataType_>, DataType_> :
        public ElementIteratorBase<Vector<DataType_>, DataType_>
    {
        private:
            /// Our wrapped iterator.
            ElementIteratorImplBase<Vector<DataType_>, DataType_> *_iterator;

        public:
            /// Constructor.
            ElementIteratorWrapper(ElementIteratorImplBase<Vector<DataType_>, DataType_> *iterator) :
                _iterator(iterator)
            {
            }

            /// Copy-constructor.
            ElementIteratorWrapper(const ElementIteratorWrapper<Vector<DataType_>, DataType_> & other) :
                _iterator(other._iterator)
            {
            }

            /// Destructor.
            ~ElementIteratorWrapper()
            {
                delete _iterator;
            }

            /// Preincrement operator.
            virtual ElementIteratorBase<Vector<DataType_>, DataType_> & operator++ ()
            {
                ++(*_iterator);
                return *this;
            }

            /// Postincrement operator.
            virtual ElementIteratorWrapper<Vector<DataType_>, DataType_> operator++ (int)
            {
                ElementIteratorWrapper<Vector<DataType_>, DataType_> result(*this);

                ++(*_iterator);

                return result;
            }

            /// Equality operator.
            virtual bool operator== (const ElementIteratorWrapper<Vector<DataType_>, DataType_> & other) const
            {
                return (*_iterator == *other._iterator);
            }

            /// Inequality operator.
            virtual bool operator!= (const ElementIteratorWrapper<Vector<DataType_>, DataType_> & other) const
            {
                return (*_iterator != *other._iterator);
            }

            /// Dereference operator 
            virtual DataType_ & operator* () const
            {
                return **_iterator;
            }

            /// Our index.
            virtual const unsigned long index() const
            {
                return _iterator->index();
            }
    };

}

#endif
