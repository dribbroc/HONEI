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

#include <libla/element_iterator.hh>
#include <libla/vector.hh>
#include <libutil/shared_array.hh>

#include <iterator>
#include <string.h>

/**
 * \file
 *
 * Implementation of DenseVector and related classes.
 *
 * \ingroup grpvector
 **/
namespace pg512 ///< \todo Namespace name?
{
    /**
     * A DenseVector is a vector with O(size) non-zero elements which keeps its data
     * sequential.
     *
     * \ingroup grpvector
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
            friend class ElementIteratorImpl<const DataType_>;

            /// Type of the const iterator over our elements.
            typedef ElementIteratorWrapper<Vector<DataType_>, DataType_, const DataType_> ConstElementIterator;

            /// Type of the iterator over our elements.
            typedef ElementIteratorWrapper<Vector<DataType_>, DataType_, DataType_> ElementIterator;

            /**
             * Constructor.
             *
             * \param size Size of the new dense vector.
             * \param offset Offset of the vector's data inside the shared array.
             * \param stepsize Stepsize between two of the vector's elements inside the shared array.
             **/
            DenseVector(const unsigned long size, unsigned long offset = 0, unsigned long stepsize = 1) :
                _elements(new DataType_[stepsize * size + offset]),
                _size(size),
                _offset(offset),
                _stepsize(stepsize)
            {
            }

            /**
             * Constructor.
             *
             * \param size Size of the new dense vector.
             * \param value Default value for all of the vector's elements.
             * \param offset Offset of the vector's data inside the shared array.
             * \param stepsize Stepsize between two of the vector's elements inside the shared array.
             **/
            DenseVector(const unsigned long size, DataType_ value, unsigned long offset = 0,
                    unsigned long stepsize = 1) :
                _elements(new DataType_[stepsize * size + offset]),
                _size(size),
                _offset(offset),
                _stepsize(stepsize)
            {
                for (unsigned long i(_offset) ; i < (_stepsize * _size + _offset) ; i += _stepsize)
                    _elements[i] = value;
            }

            /// Copy-constructor.
            DenseVector(const DenseVector<DataType_> & other) :
                _elements(other._elements),
                _size(other._size),
                _offset(other._offset),
                _stepsize(other._stepsize)
            {
            }

            /**
             * Constructor.
             *
             * \param size Size of the new dense vector.
             * \param elements SharedArray of the vector's elements.
             * \param offset Offset of the vector's data inside the shared array.
             * \param stepsize Stepsize between two of the vector's elements inside the shared array.
             **/
            DenseVector(const unsigned long size, const SharedArray<DataType_> & elements, unsigned long offset = 0,
                    unsigned stepsize = 1) :
                _elements(elements),
                _size(size),
                _offset(offset),
                _stepsize(stepsize)
            {
            }

            /// Returns const iterator pointing to the first element of the vector.
            virtual ConstElementIterator begin_elements() const
            {
                return ConstElementIterator(new ElementIteratorImpl<const DataType_>(*this, 0));
            }

            /// Returns const iterator pointing behind the last element of the vector.
            virtual ConstElementIterator end_elements() const
            {
                return ConstElementIterator(new ElementIteratorImpl<const DataType_>(*this, this->size()));
            }

            /// Returns iterator pointing to the first element of the vector.
            virtual ElementIterator begin_elements()
            {
                return ElementIterator(new ElementIteratorImpl<DataType_>(*this, 0));
            }

            /// Returns iterator pointing behind the last element of the vector.
            virtual ElementIterator end_elements() 
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
     *
     * \ingroup grpvector
     **/
    template <> template <typename DataType_> class DenseVector<DataType_>::ElementIteratorImpl<DataType_> :
        public ElementIteratorImplBase<Vector<DataType_>, DataType_>
    {
        private:
            const DenseVector<DataType_> & _vector;
            unsigned long _index;

        public:
            /**
             * Constructor.
             *
             * \param vector The parent vector that is referenced by the iterator.
             * \param index The index into the vector.
             **/
            ElementIteratorImpl(const DenseVector<DataType_> & vector, unsigned long index) :
                _vector(vector),
                _index(index)
            {
            }

            /// Copy-constructor.
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

            /// Returns our parent vector.
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

    template <> template <typename DataType_> class DenseVector<DataType_>::ElementIteratorImpl<const DataType_> :
        public ElementIteratorImplBase<Vector<DataType_>, DataType_, const DataType_>
    {
        private:
            const DenseVector<DataType_> & _vector;
            unsigned long _index;

        public:
            /**
             * Constructor.
             *
             * \param vector The parent vector that is referenced by the iterator.
             * \param index The index into the vector.
             **/
            ElementIteratorImpl(const DenseVector<DataType_> & vector, unsigned long index) :
                _vector(vector),
                _index(index)
            {
            }

            /// Copy-constructor.
            ElementIteratorImpl(ElementIteratorImpl<const DataType_> const & other) :
                _vector(other._vector),
                _index(other._index)
            {
            }

            /// Preincrement operator.
            virtual ElementIteratorImpl<const DataType_> & operator++ ()
            {
                ++_index;
                return *this;
            }

            /// Postincrement operator.
            virtual ElementIteratorImpl<const DataType_> operator++ (int)
            {
                ElementIteratorImpl<const DataType_> result(*this);
                ++_index;
                return result;
            }

            /// Equality operator.
            virtual bool operator== (const ElementIteratorImplBase<Vector<DataType_>, DataType_, const DataType_> & other) const
            {
                return (&_vector == other.parent()) && (_index == other.index());
            }

            /// Inequality operator.
            virtual bool operator!= (const ElementIteratorImplBase<Vector<DataType_>, DataType_, const DataType_> & other) const
            {
                return ((&_vector != other.parent()) || (_index != other.index()));
            }

            /// Dereference operator 
            virtual const DataType_ & operator* () const
            {
                return _vector._elements[_vector._stepsize * _index + _vector._offset];
            }

            /// Returns our parent vector.
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
