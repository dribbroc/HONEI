/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Michael Abshoff <michael.abshoff@fsmath.mathematik.uni-dortmund.de>
 * Copyright (c) 2007 Sven Mallach <sven.mallach@honei.org>
 * Copyright (c) 2007, 2009 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/la/element_iterator.hh>
#include <honei/la/sparse_vector.hh>
#include <honei/util/assertion.hh>
#include <honei/util/exception.hh>
#include <honei/util/log.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <honei/util/shared_array-impl.hh>
#include <honei/util/stringify.hh>
#include <honei/util/type_traits.hh>
#include <honei/la/vector_error.hh>

#include <string>
#include <limits>
#include <cmath>

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
        ASSERT(size + 1 >= capacity, "capacity '" + stringify(capacity) + "' exceeds size '" +
                stringify(size) + "'!");
        ASSERT(size > 0, "size is zero!");
    }

    template <typename DataType_>
    SparseVector<DataType_>::SparseVector(const SparseVector<DataType_> & other) :
        PrivateImplementationPattern<SparseVector<DataType_>, Shared>(other._imp)
    {
    }

    template <typename DataType_>
    SparseVector<DataType_>::~SparseVector()
    {
    }

    template <typename DataType_>
    typename SparseVector<DataType_>::ConstElementIterator
    SparseVector<DataType_>::begin_elements() const
    {
        return ConstElementIterator(*this, 0);
    }

    template <typename DataType_>
    typename SparseVector<DataType_>::ConstElementIterator
    SparseVector<DataType_>::end_elements() const
    {
        return ConstElementIterator(*this, this->_imp->_size);
    }

    template <typename DataType_>
    typename SparseVector<DataType_>::ElementIterator
    SparseVector<DataType_>::begin_elements()
    {
        return ElementIterator(*this, 0);
    }

    template <typename DataType_>
    typename SparseVector<DataType_>::ElementIterator
    SparseVector<DataType_>::end_elements()
    {
        return ElementIterator(*this, this->_imp->_size);
    }

    template <typename DataType_>
    typename SparseVector<DataType_>::NonZeroConstElementIterator
    SparseVector<DataType_>::begin_non_zero_elements() const
    {
        return NonZeroConstElementIterator(*this, 0);
    }

    template <typename DataType_>
    typename SparseVector<DataType_>::NonZeroConstElementIterator
    SparseVector<DataType_>::end_non_zero_elements() const
    {
        return NonZeroConstElementIterator(*this, this->_imp->_used_elements);
    }

    template <typename DataType_>
    typename SparseVector<DataType_>::NonZeroElementIterator
    SparseVector<DataType_>::begin_non_zero_elements()
    {
        return NonZeroElementIterator(*this, 0);
    }

    template <typename DataType_>
    typename SparseVector<DataType_>::NonZeroElementIterator
    SparseVector<DataType_>::end_non_zero_elements()
    {
        return NonZeroElementIterator(*this, this->_imp->_used_elements);
    }

    template <typename DataType_>
    typename SparseVector<DataType_>::NonZeroElementIterator
    SparseVector<DataType_>::non_zero_element_at(unsigned long pos)
    {
        return NonZeroElementIterator(*this, pos);
    }

    template <typename DataType_>
    typename SparseVector<DataType_>::NonZeroConstElementIterator
    SparseVector<DataType_>::non_zero_element_at(unsigned long pos) const
    {
        return NonZeroConstElementIterator(*this, pos);
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
        CONTEXT("When creating a copy:");
        SparseVector result(this->_imp->_size, this->_imp->_capacity);

        result._imp->_used_elements = this->_imp->_used_elements;
        TypeTraits<DataType_>::copy(this->_imp->_elements.get(), result._imp->_elements.get(), this->_imp->_used_elements);
        TypeTraits<unsigned long>::copy(this->_imp->_indices.get(), result._imp->_indices.get(), this->_imp->_used_elements);

        return result;
    }

    template <typename DataType_> struct Implementation<ElementIterator<storage::Sparse, container::Vector, DataType_> >
    {
        /// Our parent vector.
        SparseVector<DataType_> _vector;

        /// Our index.
        unsigned long _index;

        /// Our position in the index table.
        unsigned long _pos;

        Implementation(SparseVector<DataType_> & vector, unsigned long index, unsigned long pos) :
            _vector(vector),
            _index(index),
            _pos(pos)
        {
        }
    };

    template <typename DataType_>
    ElementIterator<storage::Sparse, container::Vector, DataType_>::ElementIterator(SparseVector<DataType_> & vector,
            unsigned long index) :
        PrivateImplementationPattern<ElementIterator<storage::Sparse, container::Vector, DataType_>, Single>(
                new Implementation<ElementIterator<storage::Sparse, container::Vector, DataType_> >(vector, index, 0))
    {
    }

    template <typename DataType_>
    ElementIterator<storage::Sparse, container::Vector, DataType_>::ElementIterator(
            const ElementIterator<storage::Sparse, container::Vector, DataType_> & other) :
        PrivateImplementationPattern<ElementIterator<storage::Sparse, container::Vector, DataType_>, Single>(
                new Implementation<ElementIterator<storage::Sparse, container::Vector, DataType_> >(*other._imp))
    {
    }

    template <typename DataType_>
    ElementIterator<storage::Sparse, container::Vector, DataType_>::~ElementIterator()
    {
    }

    template <typename DataType_>
    ElementIterator<storage::Sparse, container::Vector, DataType_> &
    ElementIterator<storage::Sparse, container::Vector, DataType_>::operator= (
            const ElementIterator<storage::Sparse, container::Vector, DataType_> & other)
    {
        if (&other == this)
            return *this;

        this->_imp->_vector = other._imp->_vector;
        this->_imp->_index = other._imp->_index;
        this->_imp->_pos = other._imp->_pos;

        return *this;
    }

    template <typename DataType_>
    ElementIterator<storage::Sparse, container::Vector, DataType_> &
    ElementIterator<storage::Sparse, container::Vector, DataType_>::operator++ ()
    {
        CONTEXT("When incrementing ElementIterator<Sparse, Vector> by one:");

        ++this->_imp->_index;
        while ((this->_imp->_pos < this->_imp->_vector._imp->_used_elements) &&
                (this->_imp->_vector._imp->_indices[this->_imp->_pos] < this->_imp->_index))
            ++this->_imp->_pos;

        return *this;
    }

    template <typename DataType_>
    ElementIterator<storage::Sparse, container::Vector, DataType_> &
    ElementIterator<storage::Sparse, container::Vector, DataType_>::operator+= (const unsigned long step)
    {
        CONTEXT("When incrementing ElementIterator<Sparse, Vector> by '" + stringify(step) + "':");

        this->_imp->_index += step;
        while ((this->_imp->_pos < this->_imp->_vector._imp->_used_elements)
                && (this->_imp->_vector._imp->_indices[this->_imp->_pos] < this->_imp->_index))
            ++this->_imp->_pos;

        return *this;
    }

    template <typename DataType_>
    DataType_ &
    ElementIterator<storage::Sparse, container::Vector, DataType_>::operator* () const
    {
        CONTEXT("When accessing assignable element at index '" + stringify(this->_imp->_index) + "':");

        if (this->_imp->_vector._imp->_indices[this->_imp->_pos] != this->_imp->_index)
            this->_imp->_vector._insert_element(this->_imp->_pos, this->_imp->_index);

        return this->_imp->_vector._imp->_elements[this->_imp->_pos];
    }

    template <typename DataType_>
    bool
    ElementIterator<storage::Sparse, container::Vector, DataType_>::operator< (
            const ElementIterator<storage::Sparse, container::Vector, DataType_> & other) const
    {
        return this->_imp->_index < other._imp->_index;
    }

    template <typename DataType_>
    bool
    ElementIterator<storage::Sparse, container::Vector, DataType_>::operator== (
            const ElementIterator<storage::Sparse, container::Vector, DataType_> & other) const
    {
        return ((this->_imp->_vector._imp == other._imp->_vector._imp) && (this->_imp->_index == other._imp->_index));
    }

    template <typename DataType_>
    bool
    ElementIterator<storage::Sparse, container::Vector, DataType_>::operator!= (
            const ElementIterator<storage::Sparse, container::Vector, DataType_> & other) const
    {
        return ((this->_imp->_vector._imp != other._imp->_vector._imp) || (this->_imp->_index != other._imp->_index));
    }

    template <typename DataType_>
    unsigned long
    ElementIterator<storage::Sparse, container::Vector, DataType_>::index() const
    {
        return this->_imp->_index;
    }

    template <typename DataType_> struct Implementation<ConstElementIterator<storage::Sparse, container::Vector, DataType_> >
    {
        /// Our parent vector.
        SparseVector<DataType_> _vector;

        /// Our index.
        unsigned long _index;

        /// Our position in the index table.
        unsigned long _pos;

        Implementation(const SparseVector<DataType_> & vector, unsigned long index, unsigned long pos) :
            _vector(vector),
            _index(index),
            _pos(pos)
        {
        }

        Implementation(const Implementation<ElementIterator<storage::Sparse, container::Vector, DataType_> > & other) :
            _vector(other._vector),
            _index(other._index),
            _pos(other._pos)
        {
        }
    };

    template <typename DataType_>
    ConstElementIterator<storage::Sparse, container::Vector, DataType_>::ConstElementIterator(const SparseVector<DataType_> & vector, unsigned long index) :
        PrivateImplementationPattern<ConstElementIterator<storage::Sparse, container::Vector, DataType_>, Single>(
                new Implementation<ConstElementIterator<storage::Sparse, container::Vector, DataType_> >(vector, index, 0))
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::Sparse, container::Vector, DataType_>::ConstElementIterator(const ConstElementIterator & other) :
        PrivateImplementationPattern<ConstElementIterator<storage::Sparse, container::Vector, DataType_>, Single>(
                new Implementation<ConstElementIterator<storage::Sparse, container::Vector, DataType_> >(*other._imp))
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::Sparse, container::Vector, DataType_>::ConstElementIterator(
            const ElementIterator<storage::Sparse, container::Vector, DataType_> & other) :
        PrivateImplementationPattern<ConstElementIterator<storage::Sparse, container::Vector, DataType_>, Single>(
                new Implementation<ConstElementIterator<storage::Sparse, container::Vector, DataType_> >(*other._imp))
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::Sparse, container::Vector, DataType_>::~ConstElementIterator()
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::Sparse, container::Vector, DataType_> &
    ConstElementIterator<storage::Sparse, container::Vector, DataType_>::operator= (
            const ConstElementIterator<storage::Sparse, container::Vector, DataType_> & other)
    {
        if (&other == this)
            return *this;

        this->_imp->_vector = other._imp->_vector;
        this->_imp->_index = other._imp->_index;
        this->_imp->_pos = other._imp->_pos;

        return *this;
    }

    template <typename DataType_>
    ConstElementIterator<storage::Sparse, container::Vector, DataType_> &
    ConstElementIterator<storage::Sparse, container::Vector, DataType_>::operator++ ()
    {
        CONTEXT("When incrementing ConstElementIterator<Sparse, Vector> by one:");

        ++this->_imp->_index;
        while ((this->_imp->_pos < this->_imp->_vector._imp->_used_elements) &&
                (this->_imp->_vector._imp->_indices[this->_imp->_pos] < this->_imp->_index))
            ++this->_imp->_pos;

        return *this;
    }

    template <typename DataType_>
    ConstElementIterator<storage::Sparse, container::Vector, DataType_> &
    ConstElementIterator<storage::Sparse, container::Vector, DataType_>::operator+= (const unsigned long step)
    {
        CONTEXT("When incrementing ConstElementIterator<Sparse, Vector> by '" + stringify(step) + "':");

        this->_imp->_index += step;
        while ((this->_imp->_pos < this->_imp->_vector._imp->_used_elements)
                && (this->_imp->_vector._imp->_indices[this->_imp->_pos] < this->_imp->_index))
            ++this->_imp->_pos;

        return *this;
    }

    template <typename DataType_>
    const DataType_ &
    ConstElementIterator<storage::Sparse, container::Vector, DataType_>::operator* () const
    {
        CONTEXT("When accessing unassignable element at index '" + stringify(this->_imp->_index) + "':");

        if (this->_imp->_vector._imp->_indices[this->_imp->_pos] != this->_imp->_index)
            return this->_imp->_vector._zero_element;
        else
            return this->_imp->_vector._imp->_elements[this->_imp->_pos];
    }


    template <typename DataType_>
    bool
    ConstElementIterator<storage::Sparse, container::Vector, DataType_>::operator< (
            const ConstElementIterator<storage::Sparse, container::Vector, DataType_> & other) const
    {
        return this->_imp->_index < other._imp->_index;
    }

    template <typename DataType_>
    bool
    ConstElementIterator<storage::Sparse, container::Vector, DataType_>::operator== (
            const ConstElementIterator<storage::Sparse, container::Vector, DataType_> & other) const
    {
        return ((this->_imp->_vector._imp == other._imp->_vector._imp) && (this->_imp->_index == other._imp->_index));
    }

    template <typename DataType_>
    bool
    ConstElementIterator<storage::Sparse, container::Vector, DataType_>::operator!= (
            const ConstElementIterator<storage::Sparse, container::Vector, DataType_> & other) const
    {
        return ((this->_imp->_vector._imp != other._imp->_vector._imp) || (this->_imp->_index != other._imp->_index));
    }

    template <typename DataType_>
    unsigned long
    ConstElementIterator<storage::Sparse, container::Vector, DataType_>::index() const
    {
        return this->_imp->_index;
    }


    template <typename DataType_> struct Implementation<ElementIterator<storage::SparseNonZero, container::Vector, DataType_> >
    {
        /// Our parent vector.
        SparseVector<DataType_> _vector;

        /// Our index.
        unsigned long _index;

        /// Our position in the index table.
        unsigned long _pos;

        Implementation(SparseVector<DataType_> & vector, unsigned long index, unsigned long pos) :
            _vector(vector),
            _index(index),
            _pos(pos)
        {
        }
    };

    template <typename DataType_>
    ElementIterator<storage::SparseNonZero, container::Vector, DataType_>::ElementIterator(SparseVector<DataType_> & vector,
            unsigned long pos) :
        PrivateImplementationPattern<ElementIterator<storage::SparseNonZero, container::Vector, DataType_>, Single>(
                new Implementation<ElementIterator<storage::SparseNonZero, container::Vector, DataType_> >
                (vector, vector._imp->_indices[pos], pos))
    {
    }

    template <typename DataType_>
    ElementIterator<storage::SparseNonZero, container::Vector, DataType_>::ElementIterator(
            const ElementIterator<storage::SparseNonZero, container::Vector, DataType_> & other) :
        PrivateImplementationPattern<ElementIterator<storage::SparseNonZero, container::Vector, DataType_>, Single>(
                new Implementation<ElementIterator<storage::SparseNonZero, container::Vector, DataType_> >(*other._imp))
    {
    }

    template <typename DataType_>
    ElementIterator<storage::SparseNonZero, container::Vector, DataType_>::~ElementIterator()
    {
    }

    template <typename DataType_>
    ElementIterator<storage::SparseNonZero, container::Vector, DataType_> &
    ElementIterator<storage::SparseNonZero, container::Vector, DataType_>::operator= (
            const ElementIterator<storage::SparseNonZero, container::Vector, DataType_> & other)
    {
        if (&other == this)
            return *this;

        this->_imp->_vector = other._imp->_vector;
        this->_imp->_index = other._imp->_index;
        this->_imp->_pos = other._imp->_pos;

        return *this;
    }

    template <typename DataType_>
    ElementIterator<storage::SparseNonZero, container::Vector, DataType_> &
    ElementIterator<storage::SparseNonZero, container::Vector, DataType_>::operator++ ()
    {
        CONTEXT("When incrementing ElementIterator<SparseNonZero, Vector> by one:");

        if (this->_imp->_vector._imp->_indices[this->_imp->_pos] == this->_imp->_index)
        {
            ++this->_imp->_pos;
        }
        else if (this->_imp->_vector._imp->_indices[this->_imp->_pos] < this->_imp->_index)
        {
            while (this->_imp->_vector._imp->_indices[this->_imp->_pos] <= this->_imp->_index)
            {
                ++this->_imp->_pos;
            }
        }
        else if (this->_imp->_vector._imp->_indices[this->_imp->_pos] > this->_imp->_index)
        {
            while (this->_imp->_vector._imp->_indices[this->_imp->_pos] > this->_imp->_index)
            {
                --this->_imp->_pos;
            }
            ++this->_imp->_pos;
        }

        this->_imp->_index = this->_imp->_vector._imp->_indices[this->_imp->_pos];
        return *this;
    }

    template <typename DataType_>
    ElementIterator<storage::SparseNonZero, container::Vector, DataType_> &
    ElementIterator<storage::SparseNonZero, container::Vector, DataType_>::operator+= (const unsigned long step)
    {
        CONTEXT("When incrementing ElementIterator<SparseNonZero, Vector> by '" + stringify(step) + "':");

        //Restore position if the vector was modified

        if (this->_imp->_vector._imp->_indices[this->_imp->_pos] < this->_imp->_index)
        {
            while (this->_imp->_vector._imp->_indices[this->_imp->_pos] < this->_imp->_index)
            {
                this->_imp->_pos++;
            }
        }
        else if (this->_imp->_vector._imp->_indices[this->_imp->_pos] > this->_imp->_index)
        {
            while (this->_imp->_vector._imp->_indices[this->_imp->_pos] > this->_imp->_index)
            {
                this->_imp->_pos--;
            }
        }

        this->_imp->_pos += step;
        if (this->_imp->_pos < this->_imp->_vector.capacity())
            this->_imp->_index = this->_imp->_vector._imp->_indices[this->_imp->_pos];
        else
            this->_imp->_index = this->_imp->_vector.size();

        return *this;
    }

    template <typename DataType_>
    DataType_ &
    ElementIterator<storage::SparseNonZero, container::Vector, DataType_>::operator* () const
    {
        CONTEXT("When accessing assignable element at index '" + stringify(this->_imp->_index) + "':");

        //Restore position if the vector was modified
        if (this->_imp->_vector._imp->_indices[this->_imp->_pos] < this->_imp->_index)
        {
            while (this->_imp->_vector._imp->_indices[this->_imp->_pos] < this->_imp->_index)
            {
                ++this->_imp->_pos;
            }
        }
        else if (this->_imp->_vector._imp->_indices[this->_imp->_pos] > this->_imp->_index)
        {
            while (this->_imp->_vector._imp->_indices[this->_imp->_pos] > this->_imp->_index)
            {
                --this->_imp->_pos;
            }
        }
        return this->_imp->_vector._imp->_elements[this->_imp->_pos];
    }

    template <typename DataType_>
    bool
    ElementIterator<storage::SparseNonZero, container::Vector, DataType_>::operator< (
            const ElementIterator<storage::SparseNonZero, container::Vector, DataType_> & other) const
    {
        return this->_imp->_index < other._imp->_index;
    }

    template <typename DataType_>
    bool
    ElementIterator<storage::SparseNonZero, container::Vector, DataType_>::operator== (
            const ElementIterator<storage::SparseNonZero, container::Vector, DataType_> & other) const
    {
        return ((this->_imp->_vector._imp == other._imp->_vector._imp) && (this->_imp->_index == other._imp->_index));
    }

    template <typename DataType_>
    bool
    ElementIterator<storage::SparseNonZero, container::Vector, DataType_>::operator!= (
            const ElementIterator<storage::SparseNonZero, container::Vector, DataType_> & other) const
    {
        return ((this->_imp->_vector._imp != other._imp->_vector._imp) || (this->_imp->_index != other._imp->_index));
    }

    template <typename DataType_>
    unsigned long
    ElementIterator<storage::SparseNonZero, container::Vector, DataType_>::index() const
    {
        return this->_imp->_index;
    }

    template <typename DataType_> struct Implementation<ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_> >
    {
        /// Our parent vector.
        SparseVector<DataType_> _vector;

        /// Our index.
        unsigned long _index;

        /// Our position in the index table.
        unsigned long _pos;

        Implementation(const SparseVector<DataType_> & vector, unsigned long index, unsigned long pos) :
            _vector(vector),
            _index(index),
            _pos(pos)
        {
        }

        Implementation(const Implementation<ElementIterator<storage::SparseNonZero, container::Vector, DataType_> > & other) :
            _vector(other._vector),
            _index(other._index),
            _pos(other._pos)
        {
        }
    };

    template <typename DataType_>
    ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_>::ConstElementIterator(const SparseVector<DataType_> & vector,
            unsigned long pos) :
        PrivateImplementationPattern<ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_>, Single>(
                new Implementation<ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_> >
                (vector, vector._imp->_indices[pos], pos))
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_>::ConstElementIterator(
            const ElementIterator<storage::SparseNonZero, container::Vector, DataType_> & other) :
        PrivateImplementationPattern<ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_>, Single>(
                new Implementation<ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_> >(*other._imp))
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_>::ConstElementIterator(
            const ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_> & other) :
        PrivateImplementationPattern<ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_>, Single>(
                new Implementation<ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_> >(*other._imp))
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_>::~ConstElementIterator()
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_> &
    ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_>::operator= (
            const ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_> & other)
    {
        if (&other == this)
            return *this;

        this->_imp->_vector = other._imp->_vector;
        this->_imp->_index = other._imp->_index;
        this->_imp->_pos = other._imp->_pos;

        return *this;
    }

    template <typename DataType_>
    ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_> &
    ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_>::operator++ ()
    {
        CONTEXT("When incrementing ConstElementIterator<SparseNonZero, Vector> by one:");

        if (this->_imp->_vector._imp->_indices[this->_imp->_pos] == this->_imp->_index)
        {
            ++this->_imp->_pos;
        }
        else if (this->_imp->_vector._imp->_indices[this->_imp->_pos] < this->_imp->_index)
        {
            while (this->_imp->_vector._imp->_indices[this->_imp->_pos] <= this->_imp->_index)
            {
                ++this->_imp->_pos;
            }
        }
        else if (this->_imp->_vector._imp->_indices[this->_imp->_pos] > this->_imp->_index)
        {
            while (this->_imp->_vector._imp->_indices[this->_imp->_pos] > this->_imp->_index)
            {
                --this->_imp->_pos;
            }
            ++this->_imp->_pos;
        }

        this->_imp->_index = this->_imp->_vector._imp->_indices[this->_imp->_pos];
        return *this;
    }

    template <typename DataType_>
    ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_> &
    ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_>::operator+= (const unsigned long step)
    {
        CONTEXT("When incrementing ConstElementIterator<SparseNonZero, Vector> by '" + stringify(step) + "':");

        //Restore position if the vector was modified

        if (this->_imp->_vector._imp->_indices[this->_imp->_pos] < this->_imp->_index)
        {
            while (this->_imp->_vector._imp->_indices[this->_imp->_pos] < this->_imp->_index)
            {
                this->_imp->_pos++;
            }
        }
        else if (this->_imp->_vector._imp->_indices[this->_imp->_pos] > this->_imp->_index)
        {
            while (this->_imp->_vector._imp->_indices[this->_imp->_pos] > this->_imp->_index)
            {
                this->_imp->_pos--;
            }
        }

        this->_imp->_pos += step;
        if (this->_imp->_pos < this->_imp->_vector.capacity())
            this->_imp->_index = this->_imp->_vector._imp->_indices[this->_imp->_pos];
        else
            this->_imp->_index = this->_imp->_vector.size();

        return *this;
    }

    template <typename DataType_>
    const DataType_ &
    ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_>::operator* () const
    {
        CONTEXT("When accessing unassignable element at position '" + stringify(this->_imp->_pos) + "':");

        //Restore position if the vector was modified
        unsigned long temp_pos(this->_imp->_pos);
        if (this->_imp->_vector._imp->_indices[temp_pos] < this->_imp->_index)
        {
            while (this->_imp->_vector._imp->_indices[temp_pos] < this->_imp->_index)
            {
                ++temp_pos;
            }
        }
        else if (this->_imp->_vector._imp->_indices[temp_pos] > this->_imp->_index)
        {
            while (this->_imp->_vector._imp->_indices[temp_pos] > this->_imp->_index)
            {
                --temp_pos;
            }
        }
        return this->_imp->_vector._imp->_elements[temp_pos];
    }

    template <typename DataType_>
    bool
    ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_>::operator< (
            const ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_> & other) const
    {
        return this->_imp->_index < other._imp->_index;
    }

    template <typename DataType_>
    bool
    ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_>::operator== (
            const ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_> & other) const
    {
        return ((this->_imp->_vector._imp == other._imp->_vector._imp) && (this->_imp->_index == other._imp->_index));
    }

    template <typename DataType_>
    bool
    ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_>::operator!= (
            const ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_> & other) const
    {
        return ((this->_imp->_vector._imp != other._imp->_vector._imp) || (this->_imp->_index != other._imp->_index));
    }

    template <typename DataType_>
    unsigned long
    ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_>::index() const
    {
        return this->_imp->_index;
    }

    template <typename DataType_>
    bool
    operator== (const SparseVector<DataType_> & a, const SparseVector<DataType_> & b)
    {
        if (a.size() != b.size())
            throw VectorSizeDoesNotMatch(a.size(), b.size());

        for (typename SparseVector<DataType_>::ConstElementIterator i(a.begin_elements()), i_end(a.end_elements()),
                j(b.begin_elements()) ; i != i_end ; ++i, ++j)
        {
            if (std::fabs(*i - *j) > std::numeric_limits<DataType_>::epsilon())
                return false;
        }

        return true;
    }

    template <typename DataType_>
    std::ostream &
    operator<< (std::ostream & lhs, const SparseVector<DataType_> & b)
    {
        lhs << "[";
        for (typename SparseVector<DataType_>::ConstElementIterator i(b.begin_elements()), i_end(b.end_elements()) ;
                i != i_end ; ++i)
        {
            lhs << "  " << *i;
        }
        lhs << "]" << std::endl;

        return lhs;
    }
}

#endif
