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

#include <libla/element_iterator.hh>
#include <libla/vector.hh>
#include <libutil/exception.hh>
#include <libutil/shared_array.hh>

#include <iterator>
#include <ostream>
#include <string>
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
     * SparseVector is a vector with O(1) non-zero elements which keeps its data
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
            static const DataType_ _zero_element;

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

            /// Our normal implementation of ElementIteratorBase.
            template <typename ElementType_> class SparseElementIterator;

            /// Our smart implementation of ElementIteratorBase.
            template <typename ElementType_> class NonZeroElementIterator;

            typedef typename Vector<DataType_>::VectorElementIterator VectorElementIterator;

        public:
            friend class SparseElementIterator<DataType_>;

            /// Type of the const iterator over our elements.
            typedef typename Vector<DataType_>::ConstElementIterator ConstElementIterator;

            /// Type of the iterator over our elements.
            typedef typename Vector<DataType_>::ElementIterator ElementIterator;

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

            /// Returns const iterator pointing to the first element of the vector.
            virtual ConstElementIterator begin_elements() const
            {
                return ConstElementIterator(new SparseElementIterator<DataType_>(*this, 0));
            }

            /// Returns iterator pointing behind the last element of the vector.
            virtual ConstElementIterator end_elements() const
            {
                return ConstElementIterator(new SparseElementIterator<DataType_>(*this, _size));
            }

            /// Returns iterator pointing to the first element of the vector.
            virtual ElementIterator begin_elements()
            {
                return ElementIterator(new SparseElementIterator<DataType_>(*this, 0));
            }

            /// Returns iterator pointing behind the last element of the vector.
            virtual ElementIterator end_elements()
            {
                return ElementIterator(new SparseElementIterator<DataType_>(*this, _size));
            }

            /// Returns const iterator pointing to the first non-zero element of the vector.
            virtual ConstElementIterator begin_non_zero_elements() const
            {
                return ConstElementIterator(new NonZeroElementIterator<DataType_>(*this, 0));
            }

            /// Returns const iterator pointing behind the last element of the vector.
            virtual ConstElementIterator end_non_zero_elements() const
            {
                return ConstElementIterator(new NonZeroElementIterator<DataType_>(*this, _used_elements));
            }

            /// Returns iterator pointing to the first non-zero element of the vector.
            virtual ElementIterator begin_non_zero_elements()
            {
                return ElementIterator(new NonZeroElementIterator<DataType_>(*this, 0));
            }

            /// Returns iterator pointing behind the last element of the vector.
            virtual ElementIterator end_non_zero_elements()
            {
                return ElementIterator(new NonZeroElementIterator<DataType_>(*this, _used_elements));
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

    template <typename DataType_> const DataType_ SparseVector<DataType_>::_zero_element = 0;

    /**
     * SparseVector::SparseElementIterator is a plain iterator implementation for sparse vectors.
     *
     * \ingroup grpvector
     **/
    template <> template <typename DataType_> class SparseVector<DataType_>::SparseElementIterator<DataType_> :
        public VectorElementIterator
    {
        private:
            /// Our parent vector.
            const SparseVector<DataType_> & _vector;

            /// Our position in the index table.
            unsigned long _pos;

            /// Our index.
            unsigned long _index;

        public:
            /// Constructors
            /// \{

            /**
             * Constructor.
             *
             * \param vector The parent vector that is referenced by the iterator.
             * \param index The index into the vector.
             **/
            SparseElementIterator(const SparseVector<DataType_> & vector, unsigned long index) :
                _vector(vector),
                _pos(0),
                _index(index)
            {
            }

            /// Copy-constructor.
            SparseElementIterator(SparseElementIterator<DataType_> const & other) :
                _vector(other._vector),
                _pos(other._pos),
                _index(other._index)
            {
            }

            /// \}

            /// Forward iterator interface
            /// \{

            /// Preincrement operator.
            virtual SparseElementIterator<DataType_> & operator++ ()
            {
                ++_index;
                while (_vector._indices[_pos] < _index)
                    ++_pos;

                return *this;
            }

            /// Dereference operator that returns an assignable reference.
            virtual DataType_ & operator* ()
            {
                if (_vector._indices[_pos] > _index)
                    throw std::string("NOT YET IMPLEMENTED!");
                else if (_vector._indices[_pos] == _index)
                    return _vector._elements[_pos];
            }

            /// Dereference operator that returns an unassignable reference.
            virtual const DataType_ & operator* () const
            {
                if (_vector._indices[_pos] > _index)
                    return _vector._zero_element;
                else if (_vector._indices[_pos] == _index)
                    return _vector._elements[_pos];
            }

            /// Comparison operator for equality.
            virtual bool operator== (const IteratorBase<DataType_, Vector<DataType_> > & other) const
            {
                return ((&_vector == other.parent()) && (_index == other.index()));
            }

            /// Comparison operator for inequality.
            virtual bool operator!= (const IteratorBase<DataType_, Vector<DataType_> > & other) const
            {
                return ((&_vector != other.parent()) || (_index != other.index()));
            }

            /// \}

            /// IteratorTraits interface
            /// \{

            /// Returns our index.
            virtual unsigned long index() const
            {
                return _index;
            }

            /// Returns a pointer to our parent container.
            virtual const Vector<DataType_> * parent() const
            {
                return &_vector;
            }

            /// \}
    };

    /**
     * SparseVector::NonZeroElementIterator is a smart iterator implementation that iterates over non-zero
     * elements of sparse vectors.
     *
     * \ingroup grpvector
     **/
    template <> template <typename DataType_> class SparseVector<DataType_>::NonZeroElementIterator<DataType_> :
        public VectorElementIterator
    {
        private:
            /// Our parent vector.
            const SparseVector<DataType_> & _vector;

            /// Our position in the index table.
            unsigned long _pos;

            /// Our index.
            unsigned long _index;

        public:
            /// Constructors
            /// \{

            /**
             * Constructor.
             *
             * \param vector The parent vector that is referenced by the iterator.
             * \param pos The index of a non-zero element into the vectors internal index table.
             **/
            NonZeroElementIterator(const SparseVector<DataType_> & vector, unsigned long pos) :
                _vector(vector),
                _pos(pos),
                _index(_vector._indices[pos])
            {
            }

            /// Copy-cnstructor.
            NonZeroElementIterator(NonZeroElementIterator<DataType_> const & other) :
                _vector(other._vector),
                _pos(other._pos),
                _index(other._index)
            {
            }

            /// \}

            /// Forward iterator interface
            /// \{

            /// Preincrement operator.
            virtual NonZeroElementIterator<DataType_> & operator++ ()
            {
                ++_pos;
                _index = _vector._indices[_pos];

                return *this;
            }

            /// Dereference operator that returns an assignable reference.
            virtual DataType_ & operator* ()
            {
                return _vector._elements[_pos];
            }

            /// Dereference operator that returns an unassignable reference.
            virtual const DataType_ & operator* () const
            {
                return _vector._elements[_pos];
            }

            /// Comparison operator for equality.
            virtual bool operator== (const IteratorBase<DataType_, Vector<DataType_> > & other) const
            {
                return ((&_vector == other.parent()) && (_index == other.index()));
            }

            /// Comparison operator for inequality.
            virtual bool operator!= (const IteratorBase<DataType_, Vector<DataType_> > & other) const
            {
                return ((&_vector != other.parent()) || (_index != other.index()));
            }

            /// \}

            /// IteratorTraits interface
            /// \{

            /// Returns our index.
            virtual unsigned long index() const
            {
                return _index;
            }

            /// Returns a pointer to our parent container.
            virtual const Vector<DataType_> * parent() const
            {
                return &_vector;
            }

            /// \}
    };
}

#endif
