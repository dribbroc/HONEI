/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_DENSE_VECTOR_IMPL_HH
#define LIBLA_GUARD_DENSE_VECTOR_IMPL_HH 1

#include <honei/la/dense_vector.hh>
#include <honei/la/dense_vector_range.hh>
#include <honei/la/sparse_vector.hh>
#include <honei/la/element_iterator.hh>
#include <honei/util/assertion.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <honei/util/shared_array.hh>
#include <honei/util/stringify.hh>
#include <honei/util/type_traits.hh>
#include <honei/util/memory_arbiter.hh>

#include <algorithm>
#include <string>

namespace honei
{
    // Forward declarations.
    template <typename DataType_> class DenseMatrix;
    template <typename DataType_> class DenseVectorRange;

    /**
     * \brief DenseVector::Implementation is the private implementation class for DenseVector.
     *
     * \ingroup grpvector
     */
    template <typename DataType_> struct Implementation<DenseVector<DataType_> >
    {
        /// Our elements.
        SharedArray<DataType_> elements;

        /// Our size.
        const unsigned long size;

        /// Our offset.
        const unsigned long offset;

        /// Our stepsize.
        const unsigned long stepsize;

        /// Constructor.
        Implementation(unsigned long s, unsigned long o, unsigned long ss) :
            elements(o + ss * (s + 1) - 1),
            size(s),
            offset(o),
            stepsize(ss)
        {
        }

        /// Constructor.
        Implementation(const SharedArray<DataType_> & e, unsigned long s, unsigned long o, unsigned long ss) :
            elements(e),
            size(s),
            offset(o),
            stepsize(ss)
        {
        }

        ~Implementation()
        {
        }
    };

    template <typename DataType_>
    DenseVector<DataType_>::DenseVector(const unsigned long size, const SharedArray<DataType_> & elements,
            unsigned long offset, unsigned stepsize) :
        PrivateImplementationPattern<DenseVector<DataType_>, Shared>(new Implementation<DenseVector<DataType_> >(elements, size, offset, stepsize))
    {
        CONTEXT("When creating DenseVector:");
        ASSERT(size > 0, "size is zero!");
    }

    template <typename DataType_>
    DenseVector<DataType_>::DenseVector(const unsigned long size, const unsigned long offset,
            const unsigned long stepsize) :
        PrivateImplementationPattern<DenseVector<DataType_>, Shared>(new Implementation<DenseVector<DataType_> >(size, offset, stepsize))
    {
        CONTEXT("When creating DenseVector:");
        ASSERT(size > 0, "size is zero!");
    }

    template <typename DataType_>
    DenseVector<DataType_>::DenseVector(DenseVector<DataType_> & src, unsigned long size) :
        PrivateImplementationPattern<DenseVector<DataType_>, Shared>(new Implementation<DenseVector<DataType_> >(src->_imp->elements, size, 0, 1))
    {
        CONTEXT("When creating DenseVector from beginning parts of another DenseVector:");
        ASSERT(size > 0, "size is zero!");
        ASSERT(size <= src.size(), "size is to big!");
    }

    template <typename DataType_>
    DenseVector<DataType_>::DenseVector(const unsigned long size, DataType_ value, unsigned long offset,
            unsigned long stepsize) :
        PrivateImplementationPattern<DenseVector<DataType_>, Shared>(new Implementation<DenseVector<DataType_> >(size, offset, stepsize))
    {
        CONTEXT("When creating DenseVector:");
        ASSERT(size > 0, "size is zero!");

        /// \todo Replace by the following once DenseVector has lost Range and Slice functionality:
        // TypeTraits<DataType_>::fill(_imp->elements.get(), size, value);

        DataType_ *  target(this->_imp->elements.get());
        for (unsigned long i(offset) ; i < (stepsize * size + offset) ; i += stepsize)
            target[i] = value;
    }

    template <typename DataType_>
    DenseVector<DataType_>::DenseVector(const SparseVector<DataType_> & other) :
        PrivateImplementationPattern<DenseVector<DataType_>, Shared>(new Implementation<DenseVector<DataType_> >(other.size(), 0, 1))
    {
        CONTEXT("When creating DenseVector form SparseVector:");
        ASSERT(other.size() > 0, "size is zero!");

        TypeTraits<DataType_>::fill(this->_imp->elements.get(), other.size(), DataType_(0));

        for (typename Vector<DataType_>::ConstElementIterator i(other.begin_non_zero_elements()),
                i_end(other.end_non_zero_elements()) ; i != i_end ; ++i)
        {
            (*this)[i.index()] = *i;
        }
    }

    template <typename DataType_>
    DenseVector<DataType_>::DenseVector(const DenseVector<DataType_> & other) :
        PrivateImplementationPattern<DenseVector<DataType_>, Shared>(other._imp)
    {
    }

    template <typename DataType_>
    DenseVector<DataType_>::~DenseVector()
    {
    }

    template <typename DataType_>
    typename Vector<DataType_>::ConstElementIterator DenseVector<DataType_>::begin_elements() const
    {
        return ConstElementIterator(new DenseElementIterator(*this, 0));
    }

    template <typename DataType_>
    typename Vector<DataType_>::ConstElementIterator DenseVector<DataType_>::end_elements() const
    {
        return ConstElementIterator(new DenseElementIterator(*this, this->_imp->size));
    }

    template <typename DataType_>
    typename Vector<DataType_>::ConstElementIterator DenseVector<DataType_>::element_at(unsigned long index) const
    {
        return ConstElementIterator(new DenseElementIterator(*this, index));
    }

    template <typename DataType_>
    typename Vector<DataType_>::ElementIterator DenseVector<DataType_>::begin_elements()
    {
        return ElementIterator(new DenseElementIterator(*this, 0));
    }

    template <typename DataType_>
    typename Vector<DataType_>::ElementIterator DenseVector<DataType_>::end_elements()
    {
        return ElementIterator(new DenseElementIterator(*this, this->_imp->size));
    }

    template <typename DataType_>
    typename Vector<DataType_>::ElementIterator DenseVector<DataType_>::element_at(unsigned long index)
    {
        return ElementIterator(new DenseElementIterator(*this, index));
    }

    template <typename DataType_>
    unsigned long DenseVector<DataType_>::size() const
    {
        return this->_imp->size;
    }

    template <typename DataType_>
    const DataType_ & DenseVector<DataType_>::operator[] (unsigned long index) const
    {
        CONTEXT("When retrieving DenseVector element, unassignable:");
        ASSERT(index < this->_imp->size && index >= 0, "index is out of bounds!");
        return this->_imp->elements[this->_imp->stepsize * index + this->_imp->offset];
    }

    template <typename DataType_>
    DataType_ & DenseVector<DataType_>::operator[] (unsigned long index)
    {
        CONTEXT("When retrieving DenseVector element, assignable:");
        ASSERT(index < this->_imp->size && index >= 0, "index is out of bounds!");
        return this->_imp->elements[this->_imp->stepsize * index + this->_imp->offset];
    }

    template <typename DataType_>
    unsigned long DenseVector<DataType_>::offset() const
    {
        return this->_imp->offset;
    }

    template <typename DataType_>
    DenseVectorRange<DataType_> DenseVector<DataType_>::range(unsigned long size, unsigned long offset) const
    {
        DenseVectorRange<DataType_> result((*this) , size, offset);
        return result;
    }

    template <typename DataType_>
    inline DataType_ * DenseVector<DataType_>::elements() const
    {
        return this->_imp->elements.get();
    }

    template <typename DataType_>
    unsigned long DenseVector<DataType_>::memid() const
    {
        return (unsigned long)this->_imp->elements.get();
    }

    template <typename DataType_>
    inline void * DenseVector<DataType_>::address() const
    {
        return this->_imp->elements.get();
    }

    template <typename DataType_>
    void * DenseVector<DataType_>::read(tags::TagValue memory) const
    {
        return MemoryArbiter::instance()->read(memory, this->memid(), this->address(), this->size() * sizeof(DataType_));
    }

    template <typename DataType_>
    void * DenseVector<DataType_>::write(tags::TagValue memory) const
    {
        return MemoryArbiter::instance()->write(memory, this->memid(), this->address(), this->size() * sizeof(DataType_));
    }

    template <typename DataType_>
    void * DenseVector<DataType_>::write_only(tags::TagValue memory) const
    {
        return MemoryArbiter::instance()->write_only(memory, this->memid(), this->address(), this->size() * sizeof(DataType_));
    }

    template <typename DataType_>
    void DenseVector<DataType_>::release_read() const
    {
        MemoryArbiter::instance()->release_read(this->memid());
    }

    template <typename DataType_>
    void DenseVector<DataType_>::release_write() const
    {
        MemoryArbiter::instance()->release_write(this->memid());
    }

    template <typename DataType_>
    DenseVector<DataType_> DenseVector<DataType_>::copy() const
    {
        DenseVector result(this->_imp->size);
        this->read();

        TypeTraits<DataType_>::copy(this->_imp->elements.get(), result._imp->elements.get(), this->_imp->size);

        this->release_read();
        return result;
    }

    /**
     * \brief DenseVector::DenseElementIterator is a simple iterator implementation for dense vectors.
     *
     * \ingroup grpvector
     **/
    template <> template <typename DataType_> class DenseVector<DataType_>::DenseElementIterator :
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

                    return _vector._imp->elements[_vector._imp->stepsize * _index + _vector._imp->offset];
                }

                /// Dereference operator that returns an unassignable reference.
                virtual const DataType_ & operator* () const
                {
                    CONTEXT("When accessing unassignable element at index '" + stringify(_index) + "':");

                    return _vector._imp->elements[_vector._imp->stepsize * _index + _vector._imp->offset];
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
