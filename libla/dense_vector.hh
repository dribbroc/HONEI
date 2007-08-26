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
#include <libutil/assertion.hh>
#include <libutil/shared_array.hh>

#include <iterator>

/**
 * \file
 *
 * Implementation of DenseVector and related classes.
 *
 * \ingroup grpvector
 **/
namespace pg512 ///< \todo Namespace name?
{
    template <typename DataType_> class DenseMatrix;

    /**
     * \brief DenseVector is a vector with O(size) non-zero elements which keeps its data
     * \brief sequential.
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

            /// Our implementation of ElementIteratorBase.
            template <typename ElementType_> class DenseElementIterator;

            typedef typename Vector<DataType_>::VectorElementIterator VectorElementIterator;

            /**
             * Constructor.
             *
             * For use by DenseMatrix.
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
                CONTEXT("When creating DenseVector:");
                ASSERT(size > 0, "size is zero!");
            }

        public:
            friend class DenseElementIterator<DataType_>;
            friend class DenseMatrix<DataType_>;

            /// Type of the const iterator over our elements.
            typedef typename Vector<DataType_>::ConstElementIterator ConstElementIterator;

            /// Type of the iterator over our elements.
            typedef typename Vector<DataType_>::ElementIterator ElementIterator;

            /// Constructors
            /// \{

            /**
             * Constructor.
             *
             * \param size Size of the new dense vector.
             * \param offset Offset of the vector's data inside the shared array.
             * \param stepsize Stepsize between two of the vector's elements inside the shared array.
             **/
            DenseVector(const unsigned long size, const unsigned long offset = 0, const unsigned long stepsize = 1) :
                _elements(stepsize * size + offset),
                _size(size),
                _offset(offset),
                _stepsize(stepsize)
            {
                CONTEXT("When creating DenseVector:");
                ASSERT(size > 0, "size is zero!");
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
                _elements(stepsize * size + offset),
                _size(size),
                _offset(offset),
                _stepsize(stepsize)
            {
                CONTEXT("When creating DenseVector:");
                ASSERT(size > 0, "size is zero!");

                for (unsigned long i(_offset) ; i < (_stepsize * _size + _offset) ; i += _stepsize)
                    _elements[i] = value;
            }

            /**
             * Constructor.
             *
             * Create a sub vector from a given source vector.
             * \param source The source vector.
             * \param start The starting point of the new vector in the old vector.
             * \param size Size of the new vector.
             **/
            DenseVector(const DenseVector<DataType_> & source, unsigned long start, unsigned long size) :
                _elements(size),
                _size(size),
                _offset(0),
                _stepsize(1)
            {
                CONTEXT("When creating DenseVector:");
                ASSERT(size > 0, "size is zero!");

                if  (start + size > source.size())
                {
                    throw VectorSizeDoesNotMatch(start + size, source.size());
                }

                for (int i = 0 ; i < size ; ++i)
                {
                    _elements[i] = source._elements[i + start];
                }

            }

            /// Copy-constructor.
            DenseVector(const DenseVector<DataType_> & other) :
                _elements(other._elements),
                _size(other._size),
                _offset(other._offset),
                _stepsize(other._stepsize)
            {
            }

            /// \}

            /// Returns const iterator pointing to the first element of the vector.
            virtual ConstElementIterator begin_elements() const
            {
                return ConstElementIterator(new DenseElementIterator<DataType_>(*this, 0));
            }

            /// Returns const iterator pointing behind the last element of the vector.
            virtual ConstElementIterator end_elements() const
            {
                return ConstElementIterator(new DenseElementIterator<DataType_>(*this, this->size()));
            }

            /// Returns const iterator pointing to a given element of the vector.
            virtual ConstElementIterator element_at(unsigned long index) const
            {
                return ConstElementIterator(new DenseElementIterator<DataType_>(*this, index));
            }

            /// Returns iterator pointing to the first element of the vector.
            virtual ElementIterator begin_elements()
            {
                return ElementIterator(new DenseElementIterator<DataType_>(*this, 0));
            }

            /// Returns iterator pointing behind the last element of the vector.
            virtual ElementIterator end_elements()
            {
                return ElementIterator(new DenseElementIterator<DataType_>(*this, this->size()));
            }

            /// Returns iterator pointing to a given element of the vector.
            virtual ElementIterator element_at(unsigned long index)
            {
                return ElementIterator(new DenseElementIterator<DataType_>(*this, index));
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

            /// Return a copy to the Vector.
            virtual DenseVector * copy() const
            {
                DenseVector * result(new DenseVector(_size));

                std::copy(begin_elements(), end_elements(), result->_elements.get());

                return result;
            }
    };

    /**
     * \brief DenseVector::DenseElementIterator is a simple iterator implementation for dense vectors.
     *
     * \ingroup grpvector
     **/
    template <> template <typename DataType_> class DenseVector<DataType_>::DenseElementIterator<DataType_> :
        public VectorElementIterator
    {
        private:
            /// Our parent vector.
            const DenseVector<DataType_> & _vector;

            /// Our index.
            unsigned long _index;

        public:
            /// \name Constructors
            /// \{

            /**
             * Constructor.
             *
             * \param vector The parent vector that is referenced by the iterator.
             * \param index The index into the vector.
             **/
            DenseElementIterator(const DenseVector<DataType_> & vector, unsigned long index) :
                _vector(vector),
                _index(index)
            {
            }

            /// Copy-constructor.
            DenseElementIterator(DenseElementIterator<DataType_> const & other) :
                _vector(other._vector),
                _index(other._index)
            {
            }

            /// Destructor.
            virtual ~DenseElementIterator()
            {
            }

            /// \}

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            virtual DenseElementIterator<DataType_> & operator++ ()
            {
                CONTEXT("When incrementing iterator by one:");

                ++_index;

                return *this;
            }

            /// In-place-add operator.
            virtual DenseElementIterator<DataType_> & operator+= (const unsigned long step)
            {
                CONTEXT("When incrementing iterator by '" + stringify(step) + "':");

                _index += step;

                return *this;
            }

            /// Dereference operator that returns an assignable reference.
            virtual DataType_ & operator* ()
            {
                CONTEXT("When accessing assignable element at index '" + stringify(_index) + "':");

                return _vector._elements[_vector._stepsize * _index + _vector._offset];
            }

            /// Dereference operator that returns an unassignable reference.
            virtual const DataType_ & operator* () const
            {
                CONTEXT("When accessing unassignable element at index '" + stringify(_index) + "':");

                return _vector._elements[_vector._stepsize * _index + _vector._offset];
            }

            /// Less-than operator.
            virtual bool operator< (const IteratorBase<DataType_, Vector<DataType_> > & other) const
            {
                return _index < other.index();
            }

            /// Equality operator.
            virtual bool operator== (const IteratorBase<DataType_, Vector<DataType_> > & other) const
            {
                return ((&_vector == other.parent()) && (_index == other.index()));
            }

            /// Inequality operator.
            virtual bool operator!= (const IteratorBase<DataType_, Vector<DataType_> > & other) const
            {
                return ((&_vector != other.parent()) || (_index != other.index()));
            }

            /// \}

            /// \name IteratorTraits interface
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
