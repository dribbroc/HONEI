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

#include "../libutil/shared_array.hh"
#include "element_iterator.hh"

#include <string.h>
#include <iterator>
#include <libwrapiter/libwrapiter_forward_iterator.hh>

namespace pg512 ///< \todo Namespace name?
{
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

        public:
            /// Our implementation of ElementIterator.
            template <typename ElementType_> class ElementIteratorImpl;
            friend class ElementIteratorImpl<DataType_>;

            /// Type of the iterator over our elements.
            typedef ElementIteratorWrapper<Vector<DataType_>, DataType_> ElementIterator;

            /// Constructor.
            DenseVector(unsigned long size) :
                _elements(new DataType_[size]),
                _size(size)
            {
            }

            /// Constructor.
            DenseVector(unsigned long size, DataType_ value) :
                _elements(new DataType_[size]),
                _size(size)
            {
                for (unsigned long i(0) ; i < size ; ++i)
                    _elements[i] = value;
            }

            /// Constructor.
            DenseVector(const DenseVector<DataType_> & other) :
                _elements(other._elements),
                _size(other._size)
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
                return _elements[index];
            }

            /// Retrieves element by index, zero-based, assignable.
            virtual DataType_ & operator[] (unsigned long index)
            {
                return _elements[index];
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
                return (&_vector == other.vector()) && (_index == other.index());
            }

            /// Inequality operator.
            virtual bool operator!= (const ElementIteratorImplBase<Vector<DataType_>, DataType_> & other) const
            {
                return ((&_vector != other.vector()) || (_index != other.index()));
            }

            /// Dereference operator 
            virtual DataType_ & operator* () const
            {
                return _vector._elements[_index];
            }

            /// Returns pointer to our vector.
            virtual const Vector<DataType_> * vector() const
            {
                return &_vector;
            }

            /// Returns our index.
            virtual const unsigned long index() const
            {
                return _index;
            }
    };

    /**
     * A SparseVector is a vector with O(1) non-zero elements which keeps its data
     * non-sequential.
     **/
    template <typename DataType_> class SparseVector :
        public Vector<DataType_>
    {
        private:
            /// Our non-zero elements.
            DataType_ *_elements;

            /// Out zero element.
            static DataType_ _zero_element;

            /// Indices of our non-zero elements.
            unsigned long *_indices;

            /// Out capacity of non-zero elements.
            unsigned long _capacity;

            /// Our number of current non-zero elements.
            unsigned long _used_elements;

            /// Our size, the maximal number of non-zero elements
            unsigned long _size;

            /// Resize our array of elements.
            void _resize(unsigned long new_capacity)
            {
                DataType_ *_new_elements(new DataType_[new_capacity]);

                memcpy(_new_elements, _elements, sizeof(DataType_) * _used_elements);
                delete[] _elements;
                _elements = _new_elements;
            }

        public:
            /// Iterator over our elements.
            typedef ElementIteratorWrapper<Vector<DataType_>, DataType_> ElementIterator;

            /// Our implementation of ElementIterator.
            template <typename ElementType_> class ElementIteratorImpl;
            friend class ElementIteratorImpl<DataType_>;

            /// Constructor.
            SparseVector(unsigned long size, unsigned long capacity) :
                _elements(new DataType_[capacity]),
                _indices(new unsigned long[capacity]),
                _used_elements(0),
                _capacity(capacity),
                _size(size)
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

            /// Returns out element capacity.
            virtual unsigned long capacity() const
            {
                return _capacity;
            }

            /// Returns our used element numer.
            virtual unsigned long used_elements() const
            {
                return _used_elements;
            }

            /// Returns our size.
            virtual unsigned long size() const
            {
                return _size;
            }

            /// Retrieves element by index, zero-based, unassignable.
            virtual const DataType_ & operator[] (unsigned long index) const
            {
                unsigned long i(0);

                for ( ; (i < _used_elements) && (_indices[i] < index) ; ++i)
                    ;

                if (_indices[i] == index)
                    return _elements[index];
                else
                    return _zero_element;
            }

            /// Retrieves (and inserts empty) element by index, zero-based, assignable.
            virtual DataType_ & operator[] (unsigned long index)
            {
                unsigned long i(0);

                for ( ; (i < _used_elements) && (_indices[i] < index) ; ++i)
                    ;

                if (_used_elements == _capacity)
                    _resize(_capacity + 10); /// \todo Implement intelligent resizing.

                if (_elements[i] != 0)
                {
                    memmove(_elements + 1 + i, _elements + i, sizeof(DataType_) * (_used_elements - i)); /// \todo Merge with resizing.
                    memmove(_indices + 1 + i, _indices + i, sizeof(unsigned long) * (_used_elements - i));
                }

                _elements[i] = DataType_(0);
                _indices[i] = index;
                ++_used_elements;

                return _elements[i];
            }
    };

    template <typename DataType_> DataType_ SparseVector<DataType_>::_zero_element = 0;

    /**
     * A SparseVector::ElementIteratorImpl is a smart iterator implementation for sparse vectors.
     **/
    template <> template <typename DataType_> class SparseVector<DataType_>::ElementIteratorImpl<DataType_> :
        public ElementIteratorImplBase<Vector<DataType_>, DataType_>
    {
        private:
            const SparseVector<DataType_> & _vector;
            unsigned long _pos;
            unsigned long _index;

        public:
            /// Constructor.
            ElementIteratorImpl(const SparseVector<DataType_> & vector, unsigned long index) :
                _vector(vector),
                _pos(0),
                _index(index)
            {
            }

            /// Constructor.
            ElementIteratorImpl(ElementIteratorImpl<DataType_> const & other) :
                _vector(other._vector),
                _pos(other._pos),
                _index(other._index)
            {
            }

            /// Preincrement operator.
            virtual ElementIteratorImpl<DataType_> & operator++ ()
            {
                ++_index;
                while (_vector._indices[_pos] < _index)
                    ++_pos;

                return *this;
            }

            /// Postincrement operator.
            virtual ElementIteratorImpl<DataType_> operator++ (int)
            {
                ElementIteratorImpl<DataType_> result(*this);

                ++_index;
                while (_vector._indices[_pos] < _index)
                    ++_pos;

                return result;
            }

            /// Equality operator.
            virtual bool operator== (const ElementIteratorImplBase<Vector<DataType_>, DataType_> & other) const
            {
                return (&_vector == other.vector()) && (_index == other.index());
            }

            /// Inequality operator.
            virtual bool operator!= (const ElementIteratorImplBase<Vector<DataType_>, DataType_> & other) const
            {
                return (&_vector != other.vector()) || (_index != other.index());
            }

            /// Dereference operator 
            virtual DataType_ & operator* () const
            {
                if (_vector._indices[_pos] > _index)
                    return _vector._zero_element;
                else if (_vector._indices[_pos] == _index)
                    return _vector._elements[_pos];
            }

            /// Returns pointer to our vector.
            virtual const Vector<DataType_> * vector() const
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
