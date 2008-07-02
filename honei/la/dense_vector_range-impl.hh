/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Volker Jung <volker.jung@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_DENSE_VECTOR_RANGE_IMPL_HH
#define LIBLA_GUARD_DENSE_VECTOR_RANGE_IMPL_HH 1

#include <honei/la/dense_vector.hh>
#include <honei/la/dense_vector-impl.hh>
#include <honei/la/dense_vector_range.hh>
#include <honei/la/element_iterator.hh>
#include <honei/util/assertion.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <honei/util/shared_array-impl.hh>
#include <honei/util/stringify.hh>

#include <algorithm>
#include <string>

namespace honei
{
    /**
     * \brief DenseVectorRange::Implementation is the private implementation class for DenseVectorRange.
     *
     * \ingroup grpvector
     */
    template <typename DataType_> struct Implementation<DenseVectorRange<DataType_> >
    {
        /// Our elements.
        SharedArray<DataType_> elements;

        /// Our size.
        const unsigned long size;

        /// Our offset.
        const unsigned long offset;

        /// Constructor.
        Implementation(const SharedArray<DataType_> & e, unsigned long s, unsigned long o) :
            elements(e),
            size(s),
            offset(o)
        {
        }
    };

    template <typename DataType_>
    DenseVectorRange<DataType_>::DenseVectorRange(const DenseVector<DataType_> & source, const unsigned long size,
            const unsigned long offset) :
        PrivateImplementationPattern<DenseVectorRange<DataType_>, Shared>(new Implementation<DenseVectorRange<DataType_> >(
                    source._imp->elements, size, offset))
    {
        CONTEXT("When creating DenseVectorRange:");
        ASSERT(size > 0, "size is zero!");
        ASSERT(size <= source._imp->size, "size of range is bigger than size of source!");
        ASSERT(offset <= source._imp->size, "offset is out of bounds!");
        ASSERT(offset + size <= source._imp->size, "end of range is out of bounds!");
    }

    template <typename DataType_>
    DenseVectorRange<DataType_>::DenseVectorRange(const SharedArray<DataType_> & e, const unsigned long size,
            const unsigned long offset) :
        PrivateImplementationPattern<DenseVectorRange<DataType_>, Shared>(new Implementation<DenseVectorRange<DataType_> >(
                    e, size, offset))
    {
        CONTEXT("When creating DenseVectorRange:");
        ASSERT(size > 0, "size is zero!");
        ASSERT(size <= e.size(), "size of range is bigger than size of source!");
        ASSERT(offset <= e.size(), "offset is out of bounds!");
        ASSERT(offset + size <= e.size(), "end of range is out of bounds!");
    }

    template <typename DataType_>
    DenseVectorRange<DataType_>::DenseVectorRange(const DenseVectorRange<DataType_> & other) :
        PrivateImplementationPattern<DenseVectorRange<DataType_>, Shared>(new Implementation<DenseVectorRange<DataType_> >(
                    other._imp->elements, other._imp->size, other._imp->offset))
    {
    }

    template <typename DataType_>
    DenseVectorRange<DataType_>::DenseVectorRange(const DenseVectorRange<DataType_> & source, const unsigned long size,
            const unsigned long offset) :
        PrivateImplementationPattern<DenseVectorRange<DataType_>, Shared>(new Implementation<DenseVectorRange<DataType_> >(
                    source._imp->elements, size, offset + source._imp->offset))
    {
        CONTEXT("When creating DenseVectorRange:");
        ASSERT(size > 0, "size is zero!");
        ASSERT(size <= source._imp->size, "size of range is bigger than size of source!");
        ASSERT(offset <= source._imp->size, "offset is out of bounds!");
        ASSERT(offset + size <= source._imp->size, "end of range is out of bounds!");
    }

    template <typename DataType_>
    DenseVectorRange<DataType_>::~DenseVectorRange()
    {
    }

    template <typename DataType_>
    typename Vector<DataType_>::ConstElementIterator DenseVectorRange<DataType_>::begin_elements() const
    {
        return ConstElementIterator(new DenseElementIterator(*this, 0));
    }

    template <typename DataType_>
    typename Vector<DataType_>::ConstElementIterator DenseVectorRange<DataType_>::end_elements() const
    {
        return ConstElementIterator(new DenseElementIterator(*this, this->_imp->size));
    }

    template <typename DataType_>
    typename Vector<DataType_>::ConstElementIterator DenseVectorRange<DataType_>::element_at(unsigned long index) const
    {
        return ConstElementIterator(new DenseElementIterator(*this, index));
    }

    template <typename DataType_>
    typename Vector<DataType_>::ElementIterator DenseVectorRange<DataType_>::begin_elements()
    {
        return ElementIterator(new DenseElementIterator(*this, 0));
    }

    template <typename DataType_>
    typename Vector<DataType_>::ElementIterator DenseVectorRange<DataType_>::end_elements()
    {
        return ElementIterator(new DenseElementIterator(*this, this->_imp->size));
    }

    template <typename DataType_>
    typename Vector<DataType_>::ElementIterator DenseVectorRange<DataType_>::element_at(unsigned long index)
    {
        return ElementIterator(new DenseElementIterator(*this, index));
    }

    template <typename DataType_>
    unsigned long DenseVectorRange<DataType_>::size() const
    {
        return this->_imp->size;
    }

    template <typename DataType_>
    const DataType_ & DenseVectorRange<DataType_>::operator[] (unsigned long index) const
    {
        CONTEXT("When retrieving DenseVectorRange element, unassignable:");
        ASSERT(index < this->_imp->size && index >= 0, "index is out of bounds!");
        return this->_imp->elements[index + this->_imp->offset];
    }

    template <typename DataType_>
    DataType_ & DenseVectorRange<DataType_>::operator[] (unsigned long index)
    {
        CONTEXT("When retrieving DenseVectorRange element, assignable:");
        ASSERT(index < this->_imp->size && index >= 0, "index is out of bounds!");
        return this->_imp->elements[index + this->_imp->offset];
    }

    template <typename DataType_>
    unsigned long DenseVectorRange<DataType_>::offset() const
    {
        return this->_imp->offset;
    }

    template <typename DataType_>
    DenseVectorRange<DataType_> DenseVectorRange<DataType_>::range(unsigned long size, unsigned long offset) const
    {
        DenseVectorRange<DataType_> result(*this, size, offset);
        return result;
    }

    template <typename DataType_>
    inline DataType_ * DenseVectorRange<DataType_>::elements() const
    {
        return this->_imp->elements.get() + this->_imp->offset;
    }

    template <typename DataType_>
    inline unsigned long DenseVectorRange<DataType_>::memid() const
    {
        return (unsigned long)this->_imp->elements.get();
    }

    template <typename DataType_>
    inline void * DenseVectorRange<DataType_>::address() const
    {
        return this->_imp->elements.get() + this->_imp->offset;
    }

    template <typename DataType_>
    void * DenseVectorRange<DataType_>::read(tags::TagValue memory) const
    {
        return MemoryArbiter::instance()->read(memory, this->memid(), this->_imp->elements.get() + this->_imp->offset, this->size() * sizeof(DataType_));
    }

    template <typename DataType_>
    void * DenseVectorRange<DataType_>::write(tags::TagValue memory) const
    {
        return MemoryArbiter::instance()->write(memory, this->memid(), this->_imp->elements.get() + this->_imp->offset, this->size() * sizeof(DataType_));
    }

    template <typename DataType_>
    void * DenseVectorRange<DataType_>::write_only(tags::TagValue memory) const
    {
        return MemoryArbiter::instance()->write_only(memory, this->memid(), this->_imp->elements.get() + this->_imp->offset, this->size() * sizeof(DataType_));
    }

    template <typename DataType_>
    void DenseVectorRange<DataType_>::release_read() const
    {
        MemoryArbiter::instance()->release_read(this->memid());
    }

    template <typename DataType_>
    void DenseVectorRange<DataType_>::release_write() const
    {
        MemoryArbiter::instance()->release_write(this->memid());
    }

    template <typename DataType_>
    DenseVector<DataType_> DenseVectorRange<DataType_>::copy() const
    {
        DenseVector<DataType_> result(this->_imp->size);
        DataType_ * source(this->_imp->elements.get());
        DataType_ * target(result.elements());
        for (unsigned long i(0) ; i < this->_imp->size ; i++)
        {
            target[i] = source[i + this->_imp->offset];
        }

        /// \todo: Use TypeTraits<DataType_>::copy()

        return result;
    }


    /**
     * \brief DenseVectorRange::DenseElementIterator is a simple iterator implementation for dense vector ranges.
     *
     * \ingroup grpvector
     */
    template <> template <typename DataType_> class DenseVectorRange<DataType_>::DenseElementIterator :
        public VectorElementIterator
        {
            private:
                /// Our parent range.
                const DenseVectorRange<DataType_> & _vector;

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
                DenseElementIterator(const DenseVectorRange<DataType_> & vector, unsigned long index) :
                    _vector(vector),
                    _index(index)
                {
                }

                /// Copy-constructor.
                DenseElementIterator(DenseElementIterator const & other) :
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
                virtual DenseElementIterator & operator++ ()
                {
                    CONTEXT("When incrementing iterator by one:");

                    ++_index;

                    return *this;
                }

                /// In-place-add operator.
                virtual DenseElementIterator & operator+= (const unsigned long step)
                {
                    CONTEXT("When incrementing iterator by '" + stringify(step) + "':");

                    _index += step;

                    return *this;
                }

                /// Dereference operator that returns an assignable reference.
                virtual DataType_ & operator* ()
                {
                    CONTEXT("When accessing assignable element at index '" + stringify(_index) + "':");

                    return _vector._imp->elements[_index + _vector._imp->offset];
                }

                /// Dereference operator that returns an unassignable reference.
                virtual const DataType_ & operator* () const
                {
                    CONTEXT("When accessing unassignable element at index '" + stringify(_index) + "':");

                    return _vector._imp->elements[_index + _vector._imp->offset];
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
