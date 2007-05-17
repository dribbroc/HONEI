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

#ifndef LIBLA_GUARD_SPARSE_VECTOR_HH
#define LIBLA_GUARD_SPARSE_VECTOR_HH 1

#include <libutil/exception.hh>
#include <libutil/shared_array.hh>
#include <libla/element_iterator.hh>

#include <iterator>
#include <ostream>
#include <string.h>

/**
 * \file
 *
 * Implementation of SparseVector and related classes.
 *
 * \ingroup grpvector
 **/
namespace pg512 ///< \todo Namespace name?
{
    /**
     * A SparseVector is a vector with O(1) non-zero elements which keeps its data
     * non-sequential.
     *
     * \ingroup grpvector
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
                _capacity = new_capacity;
            }

        public:
            /// Iterator over our elements.
            typedef ElementIteratorWrapper<Vector<DataType_>, DataType_> ElementIterator;

            /// Our implementation of ElementIterator.
            template <typename ElementType_> class ElementIteratorImpl;
            friend class ElementIteratorImpl<DataType_>;

            /// Our implementation of NonZeroElementIterator.
            template <typename ElementType_> class NonZeroElementIteratorImpl;
            friend class NonZeroElementIteratorImpl<DataType_>;

            /**
             * Constructor.
             *
             * \param size Size of the new sparse vector.
             * \param capacity Number of elements for which the vector shall reserve space.
             **/
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
                return ElementIterator(new ElementIteratorImpl<DataType_>(*this, _size));
            }

            /// Returns iterator pointing to the first non-zero element of the vector.
            virtual ElementIterator begin_non_zero_elements() const
            {
                return ElementIterator(new NonZeroElementIteratorImpl<DataType_>(*this, 0));
            }

            /// Returns iterator pointing behind the last element of the vector.
            virtual ElementIterator end_non_zero_elements() const
            {
                return ElementIterator(new NonZeroElementIteratorImpl<DataType_>(*this, _used_elements));
            }

            /// Returns out element capacity.
            virtual unsigned long capacity() const
            {
                return _capacity;
            }

            /// Returns our used element number.
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
     *
     * \ingroup grpvector
     **/
    template <> template <typename DataType_> class SparseVector<DataType_>::ElementIteratorImpl<DataType_> :
        public ElementIteratorImplBase<Vector<DataType_>, DataType_>
    {
        private:
            const SparseVector<DataType_> & _vector;
            unsigned long _pos;
            unsigned long _index;

        public:
            /**
             * Constructor.
             *
             * \param vector The parent vector that is referenced by the iterator.
             * \param index The index into the vector.
             **/
            ElementIteratorImpl(const SparseVector<DataType_> & vector, unsigned long index) :
                _vector(vector),
                _pos(0),
                _index(index)
            {
            }

            /// Copy-constructor.
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
                return (&_vector == other.parent()) && (_index == other.index());
            }

            /// Inequality operator.
            virtual bool operator!= (const ElementIteratorImplBase<Vector<DataType_>, DataType_> & other) const
            {
                return (&_vector != other.parent()) || (_index != other.index());
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

    /**
     * A SparseVector::NonZeroElementIteratorImpl is a smart iterator implementation that iterates over non-zero
     * elements of a sparse vectors.
     *
     * \ingroup grpvector
     **/
    template <> template <typename DataType_> class SparseVector<DataType_>::NonZeroElementIteratorImpl<DataType_> :
        public ElementIteratorImplBase<Vector<DataType_>, DataType_>
    {
        private:
            const SparseVector<DataType_> & _vector;
            unsigned long _pos;
            unsigned long _index;

        public:
            /**
             * Constructor.
             *
             * \param vector The parent vector that is referenced by the iterator.
             * \param pos The index of a non-zero element into the vectors internal index table.
             **/
            NonZeroElementIteratorImpl(const SparseVector<DataType_> & vector, unsigned long pos) :
                _vector(vector),
                _pos(pos),
                _index(_vector._indices[pos])
            {
            }

            /// Copy-cnstructor.
            NonZeroElementIteratorImpl(NonZeroElementIteratorImpl<DataType_> const & other) :
                _vector(other._vector),
                _pos(other._pos),
                _index(other._index)
            {
            }

            /// Preincrement operator.
            virtual NonZeroElementIteratorImpl<DataType_> & operator++ ()
            {
                ++_pos;
                _index = _vector._indices[_pos];

                return *this;
            }

            /// Postincrement operator.
            virtual NonZeroElementIteratorImpl<DataType_> operator++ (int)
            {
                NonZeroElementIteratorImpl<DataType_> result(*this);

                ++_pos;
                _index = _vector._indices[_pos];

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
                return (&_vector != other.parent()) || (_index != other.index());
            }

            /// Dereference operator 
            virtual DataType_ & operator* () const
            {
                return _vector._elements[_pos];
            }

            /// Returns pointer to our parent.
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
