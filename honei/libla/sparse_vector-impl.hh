/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Michael Abshoff <michael.abshoff@fsmath.mathematik.uni-dortmund.de>
 * Copyright (c) 2007 Sven Mallach <sven.mallach@honei.org>
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_SPARSE_VECTOR_IMPL_HH
#define LIBLA_GUARD_SPARSE_VECTOR_IMPL_HH 1

#include <honei/libla/element_iterator.hh>
#include <honei/libla/sparse_vector.hh>
#include <honei/libutil/assertion.hh>
#include <honei/libutil/exception.hh>
#include <honei/libutil/log.hh>
#include <honei/libutil/private_implementation_pattern-impl.hh>
#include <honei/libutil/shared_array-impl.hh>
#include <honei/libutil/stringify.hh>
#include <honei/libutil/type_traits.hh>

#include <iterator>
#include <string>

namespace honei
{
    template <typename DataType_> const DataType_ SparseVector<DataType_>::_zero_element = 0;

    /**
     * \brief SparseVector::Implementation is the private implementation class for SparseVector.
     *
     * \ingroup grpvector
     */
    template <typename DataType_> struct Implementation<SparseVector<DataType_> >
    {
        /// Our non-zero elements.
        SharedArray<DataType_> _elements;

        /// Our indices of non-zero elements.
        SharedArray<unsigned long> _indices;

        /// Out capacity of non-zero elements.
        unsigned long _capacity;

        /// Our size, the maximal number of non-zero elements
        const unsigned long _size;

        /// Our number of current non-zero elements.
        unsigned long _used_elements;

        /// \name Constructors
        /// \{

        /**
         * Constructor.
         *
         * \param size Size of the new SparseVector.
         * \param capacity Capacity of elements that can be held without resizing.
         */
        Implementation(unsigned long size, unsigned long capacity) :
            _elements(capacity),
            _indices(capacity),
            _capacity(capacity),
            _size(size),
            _used_elements(0)
        {
            CONTEXT("When creating SparseVector::Implementation:");
            ASSERT(capacity > 0, "capacity is zero!");

            // Sneak in 'terminating elements', as index can never be size.
            _elements[0] = DataType_(0);
            TypeTraits<unsigned long>::fill(_indices.get(), capacity, size);
        }

        /**
         * Destructor.
         */
        ~Implementation()
        {
        }

        /// \}
    };

    template <typename DataType_>
    void SparseVector<DataType_>::_insert_element(unsigned long position, unsigned long index) const
    {
        CONTEXT("When inserting element at position '" + stringify(position) + "' with index '" +
                stringify(index) + "':");

        bool realloc(this->_imp->_capacity <= this->_imp->_used_elements + 1);
        unsigned long capacity(realloc ? std::min(this->_imp->_capacity + 10, this->_imp->_size + 1) : this->_imp->_capacity);
        DataType_ * elements(realloc ? new DataType_[capacity] : this->_imp->_elements.get());
        unsigned long * indices(realloc ? new unsigned long[capacity] : this->_imp->_indices.get());

        ASSERT(position < capacity, "position '" + stringify(position) + "' out of bounds!");
        ASSERT(index < this->_imp->_size, "index '" + stringify(index) + "' out of bounds!");

        if (realloc)
        {
            // Write out the terminating elements.
            TypeTraits<unsigned long>::fill(indices, capacity, this->_imp->_size);

            TypeTraits<DataType_>::copy(this->_imp->_elements.get(), elements, position + 1);
            TypeTraits<unsigned long>::copy(this->_imp->_indices.get(), indices, position + 1);
        }

        // Relies on capactiy >= used_elements + 1.
        std::copy_backward(this->_imp->_elements.get() + position, this->_imp->_elements.get() + this->_imp->_used_elements,
                elements + this->_imp->_used_elements + 1);
        std::copy_backward(this->_imp->_indices.get() + position, this->_imp->_indices.get() + this->_imp->_used_elements,
                indices + this->_imp->_used_elements + 1);

        ++this->_imp->_used_elements;

        if (realloc)
        {
            this->_imp->_elements.reset(capacity, elements);
            this->_imp->_indices.reset(capacity, indices);
            this->_imp->_capacity = capacity;
        }

        // Set new element's index and reset it to zero.
        this->_imp->_indices[position] = index;
        this->_imp->_elements[position] = DataType_(0);
    }

    template <typename DataType_>
    SparseVector<DataType_>::SparseVector(unsigned long size, unsigned long capacity) :
        PrivateImplementationPattern<SparseVector<DataType_>, Shared>(new Implementation<SparseVector<DataType_> >(size, capacity))
    {
        CONTEXT("When creating SparseVector:");
        ASSERT(size >= capacity, "capacity '" + stringify(capacity) + "' exceeds size '" +
                stringify(size) + "'!");
        ASSERT(size > 0, "size is zero!");
    }

    template <typename DataType_>
    SparseVector<DataType_>::~SparseVector()
    {
    }

    template <typename DataType_>
    typename Vector<DataType_>::ConstElementIterator SparseVector<DataType_>::begin_elements() const
    {
        return ConstElementIterator(new SparseElementIterator(*this, 0));
    }

    template <typename DataType_>
    typename Vector<DataType_>::ConstElementIterator SparseVector<DataType_>::end_elements() const
    {
        return ConstElementIterator(new SparseElementIterator(*this, this->_imp->_size));
    }

    template <typename DataType_>
    typename Vector<DataType_>::ElementIterator SparseVector<DataType_>::begin_elements()
    {
        return ElementIterator(new SparseElementIterator(*this, 0));
    }

    template <typename DataType_>
    typename Vector<DataType_>::ElementIterator SparseVector<DataType_>::end_elements()
    {
        return ElementIterator(new SparseElementIterator(*this, this->_imp->_size));
    }

    template <typename DataType_>
    typename Vector<DataType_>::ConstElementIterator SparseVector<DataType_>::begin_non_zero_elements() const
    {
        return ConstElementIterator(new NonZeroElementIterator(*this, 0));
    }

    template <typename DataType_>
    typename Vector<DataType_>::ConstElementIterator SparseVector<DataType_>::end_non_zero_elements() const
    {
        return ConstElementIterator(new NonZeroElementIterator(*this, this->_imp->_used_elements));
    }

    template <typename DataType_>
    typename Vector<DataType_>::ElementIterator SparseVector<DataType_>::begin_non_zero_elements()
    {
        return ElementIterator(new NonZeroElementIterator(*this, 0));
    }

    template <typename DataType_>
    typename Vector<DataType_>::ElementIterator SparseVector<DataType_>::end_non_zero_elements()
    {
        return ElementIterator(new NonZeroElementIterator(*this, this->_imp->_used_elements));
    }

    template <typename DataType_>
    typename Vector<DataType_>::ElementIterator SparseVector<DataType_>::non_zero_element_at(unsigned long pos)
    {
        return ElementIterator(new NonZeroElementIterator(*this, pos));
    }

    template <typename DataType_>
    typename Vector<DataType_>::ConstElementIterator SparseVector<DataType_>::non_zero_element_at(unsigned long pos) const
    {
        return ConstElementIterator(new NonZeroElementIterator(*this, pos));
    }

    template <typename DataType_>
    unsigned long SparseVector<DataType_>::capacity() const
    {
        return this->_imp->_capacity;
    }

    template <typename DataType_>
    unsigned long SparseVector<DataType_>::used_elements() const
    {
        return this->_imp->_used_elements;
    }

    template <typename DataType_>
    unsigned long SparseVector<DataType_>::size() const
    {
        return this->_imp->_size;
    }

    template <typename DataType_>
    inline DataType_ * SparseVector<DataType_>::elements() const
    {
        return this->_imp->_elements.get();
    }

    template <typename DataType_>
    inline unsigned long * SparseVector<DataType_>::indices() const
    {
        return this->_imp->_indices.get();
    }

    template <typename DataType_>
    const DataType_ & SparseVector<DataType_>::operator[] (unsigned long index) const
    {
        CONTEXT("When accessing unassignable element at index '" + stringify(index) + "':");
        ASSERT(index < this->_imp->_size, "index '" + stringify(index) + "' exceeds size '" +
                stringify(this->_imp->_size) + "'!");

        unsigned long i(0);

        for ( ; (i < this->_imp->_used_elements) && (this->_imp->_indices[i] < index) ; ++i)
            ;

        if (this->_imp->_indices[i] == index)
            return this->_imp->_elements[i];
        else
            return _zero_element;
    }

    template <typename DataType_>
    DataType_ & SparseVector<DataType_>::operator[] (unsigned long index)
    {
        CONTEXT("When accessing assignable element at index '" + stringify(index) + "':");
        ASSERT(index < this->_imp->_size, "index '" + stringify(index) + "' exceeds size '" +
                stringify(this->_imp->_size) + "'!");

        unsigned long i(0);

        for ( ; (i < this->_imp->_used_elements) && (this->_imp->_indices[i] < index) ; ++i)
            ;

        if (this->_imp->_indices[i] != index)
            _insert_element(i, index);

        return this->_imp->_elements[i];
    }

    template <typename DataType_>
    SparseVector<DataType_> SparseVector<DataType_>::copy() const
    {
        CONTEXT("When creating a copy:'");
        SparseVector result(this->_imp->_size, this->_imp->_capacity);

        result._imp->_used_elements = this->_imp->_used_elements;
        TypeTraits<DataType_>::copy(this->_imp->_elements.get(), result._imp->_elements.get(), this->_imp->_used_elements);
        TypeTraits<unsigned long>::copy(this->_imp->_indices.get(), result._imp->_indices.get(), this->_imp->_used_elements);

        return result;
    }

    /**
     * \brief SparseVector::SparseElementIterator is a plain iterator implementation for sparse vectors.
     *
     * \ingroup grpvector
     */
    template <> template <typename DataType_> class SparseVector<DataType_>::SparseElementIterator :
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
            /// \name Constructors
            /// \{

            /**
             * Constructor.
             *
             * \param vector The parent vector that is referenced by the iterator.
             * \param index The index into the vector.
             */
            SparseElementIterator(const SparseVector<DataType_> & vector, unsigned long index) :
                _vector(vector),
                _pos(0),
                _index(index)
            {
            }

            /// Copy-constructor.
            SparseElementIterator(const SparseElementIterator & other) :
                _vector(other._vector),
                _pos(other._pos),
                _index(other._index)
            {
            }

            /// Destructor
            virtual ~SparseElementIterator()
            {
            }

            /// \}

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            virtual SparseElementIterator & operator++ ()
            {
                CONTEXT("When incrementing iterator by one:");

                ++_index;
                while ((_pos < _vector._imp->_used_elements) && (_vector._imp->_indices[_pos] < _index))
                    ++_pos;

                return *this;
            }

            /// In-place-add operator.
            virtual SparseElementIterator & operator+= (const unsigned long step)
            {
                CONTEXT("When incrementing iterator by '" + stringify(step) + "':");

                _index += step;
                while ((_pos < _vector._imp->_used_elements) && (_vector._imp->_indices[_pos] < _index))
                    ++_pos;

                return *this;
            }

            /// Dereference operator that returns an assignable reference.
            virtual DataType_ & operator* ()
            {
                CONTEXT("When accessing assignable element at position '" + stringify(_pos) + "':");

                if (_vector._imp->_indices[_pos] != _index)
                    _vector._insert_element(_pos, _index);

                return _vector._imp->_elements[_pos];
            }

            /// Dereference operator that returns an unassignable reference.
            virtual const DataType_ & operator* () const
            {
                CONTEXT("When accessing unassignable element at position '" + stringify(_pos) + "':");

                if (_vector._imp->_indices[_pos] != _index)
                    return _vector._zero_element;
                else
                    return _vector._imp->_elements[_pos];
            }

            /// Less-than operator.
            virtual bool operator< (const IteratorBase<DataType_, Vector<DataType_> > & other) const
            {
                return _index < other.index();
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

    /**
     * \brief SparseVector::NonZeroElementIterator is a smart iterator implementation that iterates over non-zero
     * \brief elements of sparse vectors.
     *
     * \ingroup grpvector
     */
    template <> template <typename DataType_> class SparseVector<DataType_>::NonZeroElementIterator :
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
            /// \name Constructors
            /// \{

            /**
             * Constructor.
             *
             * \param vector The parent vector that is referenced by the iterator.
             * \param pos The index of a non-zero element into the vectors internal index table.
             */
            NonZeroElementIterator(const SparseVector<DataType_> & vector, unsigned long pos) :
                _vector(vector),
                _pos(pos),
                _index(_vector._imp->_indices[pos])
            {
            }

            /// Copy-cnstructor.
            NonZeroElementIterator(NonZeroElementIterator const & other) :
                _vector(other._vector),
                _pos(other._pos),
                _index(other._index)
            {
            }

            /// Destructor.
            virtual ~NonZeroElementIterator()
            {
            }

            /// \}

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            virtual NonZeroElementIterator & operator++ ()
            {
                CONTEXT("When incrementing iterator by one:");
                if (_vector._imp->_indices[_pos] == _index)
                {
                    ++_pos;
                }
                else if (_vector._imp->_indices[_pos] < _index)
                {
                    while (_vector._imp->_indices[_pos] <= _index)
                    {
                        ++_pos;
                    }
                }
                else if (_vector._imp->_indices[_pos] > _index)
                {
                    while (_vector._imp->_indices[_pos] > _index)
                    {
                        --_pos;
                    }
                    ++_pos;
                }

                _index = _vector._imp->_indices[_pos];
                return *this;
            }

            /// In-place-add operator.
            virtual NonZeroElementIterator & operator+= (const unsigned long step)
            {
                CONTEXT("When incrementing iterator by '" + stringify(step) + "':");

                //Restore position if the vector was modified

                if (_vector._imp->_indices[_pos] < _index)
                {
                    while (_vector._imp->_indices[_pos] < _index)
                    {
                        _pos++;
                    }
                }
                else if (_vector._imp->_indices[_pos] > _index)
                {
                    while (_vector._imp->_indices[_pos] > _index)
                    {
                        _pos--;
                    }
                }

                _pos += step;
                if (_pos < _vector.capacity())
                    _index = _vector._imp->_indices[_pos];
                else
                    _index = _vector.size();

                return *this;
            }

            /// Dereference operator that returns an assignable reference.
            virtual DataType_ & operator* ()
            {
                CONTEXT("When accessing assignable element at position '" + stringify(_pos) + "':");

                //Restore position if the vector was modified
                if (_vector._imp->_indices[_pos] < _index)
                {
                    while (_vector._imp->_indices[_pos] < _index)
                    {
                        ++_pos;
                    }
                }
                else if (_vector._imp->_indices[_pos] > _index)
                {
                    while (_vector._imp->_indices[_pos] > _index)
                    {
                        --_pos;
                    }
                }
                return _vector._imp->_elements[_pos];
            }

            /// Dereference operator that returns an unassignable reference.
            virtual const DataType_ & operator* () const
            {
                CONTEXT("When accessing unassignable element at position '" + stringify(_pos) + "':");

                //Restore position if the vector was modified
                unsigned long temp_pos(_pos);
                if (_vector._imp->_indices[temp_pos] < _index)
                {
                    while (_vector._imp->_indices[temp_pos] < _index)
                    {
                        ++temp_pos;
                    }
                }
                else if (_vector._imp->_indices[temp_pos] > _index)
                {
                    while (_vector._imp->_indices[temp_pos] > _index)
                    {
                        --temp_pos;
                    }
                }
                return _vector._imp->_elements[temp_pos];
            }

            /// Comparison operator for less-than.
            virtual bool operator< (const IteratorBase<DataType_, Vector<DataType_> > & other) const
            {
                return _index < other.index();
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
