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

#ifndef LIBLA_GUARD_DENSE_VECTOR_HH
#define LIBLA_GUARD_DENSE_VECTOR_HH 1

#include <libutil/exception.hh>
#include <libutil/shared_array.hh>
#include <libla/element_iterator.hh>

#include <iterator>
#include <string.h>

namespace pg512 ///< \todo Namespace name?
{
    /**
     * A DenseVector is a vector with O(size) non-zero elements which keeps its data
     * sequential.
     **/
    template <typename DataType_> class DenseVector :
        public Vector<DataType_>
    {
        private:
            /// Pointer to our elements.
            SharedArray<DataType_> _elements;

            /// Our size.
            unsigned long _size;

            /// Our offset.
            unsigned long _offset;

            /// Our stepsize.
            unsigned long _stepsize;

        public:
            /// Our implementation of ElementIterator.
            template <typename ElementType_> class ElementIteratorImpl;
            friend class ElementIteratorImpl<DataType_>;

            /// Type of the iterator over our elements.
            typedef ElementIteratorWrapper<Vector<DataType_>, DataType_> ElementIterator;

            /// Constructor.
            DenseVector(const unsigned long size, unsigned long offset = 0, unsigned long stepsize = 1) :
                _elements(new DataType_[stepsize * size + offset]),
                _size(size),
                _offset(offset),
                _stepsize(stepsize)
            {
            }

            /// Constructor.
            DenseVector(const unsigned long size, DataType_ value, unsigned long offset = 0,
                    unsigned long stepsize = 1) :
                _elements(new DataType_[_stepsize * size + _offset]),
                _size(size),
                _offset(offset),
                _stepsize(stepsize)
            {
                for (unsigned long i(_offset) ; i < size ; i += _stepsize)
                    _elements[i] = value;
            }

            /// Constructor.
            DenseVector(const DenseVector<DataType_> & other) :
                _elements(other._elements),
                _size(other._size)
            {
            }

            /// Constructor.
            DenseVector(const unsigned long size, const SharedArray<DataType_> & elements, unsigned long offset = 0,
                    unsigned stepsize = 1) :
                _elements(elements),
                _size(size),
                _offset(offset),
                _stepsize(stepsize)
            {
            }

            /// Returns iterator pointing to the first element of the vector.
            virtual ElementIterator begin_elements() const
            {
                return ElementIterator(new ElementIteratorImpl<DataType_>(*this, 0));
            }

            /// Returns iterator pointing behind the last element of the vector.
            virtual ElementIterator end_elements() const
            {
                return ElementIterator(new ElementIteratorImpl<DataType_>(*this, this->size()));
            }

            /// Returns our size.
            virtual unsigned long size() const
            {
                return _size;
            }

            /// Retrieves element by index, zero-based, assignable.
            virtual const DataType_ & operator[] (unsigned long index) const
            {
                return _elements[_stepsize * index + _offset];
            }

            /// Retrieves element by index, zero-based, assignable.
            virtual DataType_ & operator[] (unsigned long index)
            {
                return _elements[_stepsize * index + _offset];
            }
    };

    /**
     * A DenseVector::ElementIteratorImpl is a simple iterator implementation for dense vectors.
     **/
    template <> template <typename DataType_> class DenseVector<DataType_>::ElementIteratorImpl<DataType_> :
        public ElementIteratorImplBase<Vector<DataType_>, DataType_>
    {
        private:
            const DenseVector<DataType_> & _vector;
            unsigned long _index;

        public:
            /// Constructor.
            ElementIteratorImpl(const DenseVector<DataType_> & vector, unsigned long index) :
                _vector(vector),
                _index(index)
            {
            }

            /// Constructor.
            ElementIteratorImpl(ElementIteratorImpl<DataType_> const & other) :
                _vector(other._vector),
                _index(other._index)
            {
            }

            /// Preincrement operator.
            virtual ElementIteratorImpl<DataType_> & operator++ ()
            {
                ++_index;
                return *this;
            }

            /// Postincrement operator.
            virtual ElementIteratorImpl<DataType_> operator++ (int)
            {
                ElementIteratorImpl<DataType_> result(*this);
                ++_index;
                return result;
            }

            /// Equality operator.
            virtual bool operator== (const ElementIteratorImplBase<Vector<DataType_>, DataType_> & other) const
            {
                return (&_vector == other.parent()) && (_index == other.index());
            }

            /// Inequality operator.
            virtual bool operator!= (const ElementIteratorImplBase<Vector<DataType_>, DataType_> & other) const
            {
                return ((&_vector != other.parent()) || (_index != other.index()));
            }

            /// Dereference operator 
            virtual DataType_ & operator* () const
            {
                return _vector._elements[_vector._stepsize * _index + _vector._offset];
            }

            /// Returns pointer to our vector.
            virtual const Vector<DataType_> * parent() const
            {
                return &_vector;
            }

            /// Returns our index.
            virtual const unsigned long index() const
            {
                return _index;
            }
    };
}

#endif
