/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Michael Abshoff <michael.abshoff@fsmath.mathematik.uni-dortmund.de>
 * Copyright (c) 2007 Sven Mallach <mallach@honei.org>
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

#pragma once
#ifndef LIBLA_GUARD_DENSE_VECTOR_IMPL_HH
#define LIBLA_GUARD_DENSE_VECTOR_IMPL_HH 1

#include <honei/la/dense_vector.hh>
#include <honei/la/dense_vector_range.hh>
#include <honei/la/vector_error.hh>
#include <honei/la/element_iterator.hh>
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
#include <limits>
#include <cmath>

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
    DenseVector<DataType_>::DenseVector(const unsigned long size, const SharedArray<DataType_> & elements) :
        PrivateImplementationPattern<DenseVector<DataType_>, Shared>(new Implementation<DenseVector<DataType_> >(elements, size, 0, 1))
    {
        CONTEXT("When creating DenseVector:");
        ASSERT(size > 0, "size is zero!");
    }

    template <typename DataType_>
    DenseVector<DataType_>::DenseVector(const unsigned long size) :
        PrivateImplementationPattern<DenseVector<DataType_>, Shared>(new Implementation<DenseVector<DataType_> >(size, 0, 1))
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
    DenseVector<DataType_>::DenseVector(const unsigned long size, DataType_ value) :
        PrivateImplementationPattern<DenseVector<DataType_>, Shared>(new Implementation<DenseVector<DataType_> >(size, 0, 1))
    {
        CONTEXT("When creating DenseVector:");
        ASSERT(size > 0, "size is zero!");

        TypeTraits<DataType_>::fill(this->_imp->elements.get(), size, value);
    }

    template <typename DataType_>
    DenseVector<DataType_>::DenseVector(const SparseVector<DataType_> & other) :
        PrivateImplementationPattern<DenseVector<DataType_>, Shared>(new Implementation<DenseVector<DataType_> >(other.size(), 0, 1))
    {
        CONTEXT("When creating DenseVector form SparseVector:");
        ASSERT(other.size() > 0, "size is zero!");

        TypeTraits<DataType_>::fill(this->_imp->elements.get(), other.size(), DataType_(0));

        for (typename SparseVector<DataType_>::NonZeroConstElementIterator i(other.begin_non_zero_elements()),
                i_end(other.end_non_zero_elements()) ; i != i_end ; ++i)
        {
            (*this)[i.index()] = *i;
        }
    }

    template <typename DataType_>
    DenseVector<DataType_>::DenseVector(const DenseVector<DataType_> & other) :
        DenseVectorContinuousBase<DataType_>(),
        PrivateImplementationPattern<DenseVector<DataType_>, Shared>(other._imp)
    {
    }

    template <typename DataType_>
    DenseVector<DataType_>::~DenseVector()
    {
    }

    template <typename DataType_>
    typename DenseVector<DataType_>::ConstElementIterator
    DenseVector<DataType_>::begin_elements() const
    {
        return ConstElementIterator(this->_imp->elements, 0, 0, 1);
    }

    template <typename DataType_>
    typename DenseVector<DataType_>::ConstElementIterator
    DenseVector<DataType_>::end_elements() const
    {
        return ConstElementIterator(this->_imp->elements, this->_imp->size, 0, 1);
    }

    template <typename DataType_>
    typename DenseVector<DataType_>::ConstElementIterator
    DenseVector<DataType_>::element_at(unsigned long index) const
    {
        return ConstElementIterator(this->_imp->elements, index, 0, 1);
    }

    template <typename DataType_>
    typename DenseVector<DataType_>::ElementIterator
    DenseVector<DataType_>::begin_elements()
    {
        return ElementIterator(this->_imp->elements, 0, 0, 1);
    }

    template <typename DataType_>
    typename DenseVector<DataType_>::ElementIterator
    DenseVector<DataType_>::end_elements()
    {
        return ElementIterator(this->_imp->elements, this->_imp->size, 0, 1);
    }

    template <typename DataType_>
    typename DenseVector<DataType_>::ElementIterator
    DenseVector<DataType_>::element_at(unsigned long index)
    {
        return ElementIterator(this->_imp->elements, index, 0, 1);
    }

    template <typename DataType_>
    unsigned long DenseVector<DataType_>::size() const
    {
        return this->_imp->size;
    }

    template <typename DataType_>
    inline const DataType_ & DenseVector<DataType_>::operator[] (unsigned long index) const
    {
        CONTEXT("When retrieving DenseVector element, unassignable:");
        ASSERT(index < this->_imp->size, "index is out of bounds!");
        return this->_imp->elements[index];
    }

    template <typename DataType_>
    inline DataType_ & DenseVector<DataType_>::operator[] (unsigned long index)
    {
        CONTEXT("When retrieving DenseVector element, assignable:");
        ASSERT(index < this->_imp->size, "index is out of bounds!");
        return this->_imp->elements[index];
    }

    template <typename DataType_>
    unsigned long DenseVector<DataType_>::offset() const
    {
        return 0;
    }

    template <typename DataType_>
    unsigned long DenseVector<DataType_>::stepsize() const
    {
        return 1;
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
    inline SharedArray<DataType_> & DenseVector<DataType_>::array() const
    {
        return this->_imp->elements;
    }

    template <typename DataType_>
    void * DenseVector<DataType_>::memid() const
    {
        return this->_imp->elements.get();
    }

    template <typename DataType_>
    inline void * DenseVector<DataType_>::address() const
    {
        return this->_imp->elements.get();
    }

    template <typename DataType_>
    void * DenseVector<DataType_>::lock(LockMode mode, tags::TagValue memory) const
    {
        return MemoryArbiter::instance()->lock(mode, memory, this->memid(), this->address(), this->size() * sizeof(DataType_));
    }

    template <typename DataType_>
    void DenseVector<DataType_>::unlock(LockMode mode) const
    {
        MemoryArbiter::instance()->unlock(mode, this->memid());
    }

    template <typename DataType_>
    DenseVector<DataType_> DenseVector<DataType_>::copy() const
    {
        DenseVector result(this->_imp->size);
        this->lock(lm_read_only);

        TypeTraits<DataType_>::copy(this->_imp->elements.get(), result._imp->elements.get(), this->_imp->size);

        this->unlock(lm_read_only);
        return result;
    }

    template <typename DataType_> struct Implementation<ConstElementIterator<storage::Dense, container::Vector, DataType_> >
    {
        SharedArray<DataType_> elements;

        unsigned long index;

        unsigned long offset;

        unsigned long stepsize;

        Implementation(const SharedArray<DataType_> & elements, unsigned long index, unsigned long offset, unsigned long stepsize) :
            elements(elements),
            index(index),
            offset(offset),
            stepsize(stepsize)
        {
        }

        Implementation(const Implementation<ElementIterator<storage::Dense, container::Vector, DataType_> > & other) :
            elements(other.elements),
            index(other.index),
            offset(other.offset),
            stepsize(other.stepsize)
        {
        }
    };

    template <typename DataType_>
    ConstElementIterator<storage::Dense, container::Vector, DataType_>::ConstElementIterator(const SharedArray<DataType_> & elements,
            unsigned long index, unsigned long offset, unsigned long stepsize) :
        PrivateImplementationPattern<ConstElementIterator<storage::Dense, container::Vector, DataType_>, Single>(
                new Implementation<ConstElementIterator<storage::Dense, container::Vector, DataType_> >(elements, index, offset, stepsize))
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::Dense, container::Vector, DataType_>::ConstElementIterator(const ConstElementIterator & other) :
        PrivateImplementationPattern<ConstElementIterator<storage::Dense, container::Vector, DataType_>, Single>(
                new Implementation<ConstElementIterator<storage::Dense, container::Vector, DataType_> >(*other._imp))
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::Dense, container::Vector, DataType_>::ConstElementIterator(
            const ElementIterator<storage::Dense, container::Vector, DataType_> & other) :
        PrivateImplementationPattern<ConstElementIterator<storage::Dense, container::Vector, DataType_>, Single>(
                new Implementation<ConstElementIterator<storage::Dense, container::Vector, DataType_> >(*other._imp))
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::Dense, container::Vector, DataType_>::~ConstElementIterator()
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::Dense, container::Vector, DataType_> &
    ConstElementIterator<storage::Dense, container::Vector, DataType_>::operator= (
            const ConstElementIterator<storage::Dense, container::Vector, DataType_> & other)
    {
        if (&other == this)
            return *this;

        this->_imp->elements = other._imp->elements;
        this->_imp->index = other._imp->index;
        this->_imp->offset = other._imp->offset;
        this->_imp->stepsize = other._imp->stepsize;

        return *this;
    }

    template <typename DataType_>
    ConstElementIterator<storage::Dense, container::Vector, DataType_> &
    ConstElementIterator<storage::Dense, container::Vector, DataType_>::operator++ ()
    {
        CONTEXT("When incrementing ConstElementIterator<Dense, Vector> by one:");

        ++this->_imp->index;

        return *this;
    }

    template <typename DataType_>
    ConstElementIterator<storage::Dense, container::Vector, DataType_> &
    ConstElementIterator<storage::Dense, container::Vector, DataType_>::operator+= (const unsigned long step)
    {
        CONTEXT("When incrementing ConstElementIterator<Dense, Vector> by '" + stringify(step) + "':");

        this->_imp->index += step;

        return *this;
    }

    template <typename DataType_>
    const DataType_ &
    ConstElementIterator<storage::Dense, container::Vector, DataType_>::operator* () const
    {
        CONTEXT("When accessing unassignable element at index '" + stringify(this->_imp->index) + "':");

        return this->_imp->elements[this->_imp->stepsize * this->_imp->index + this->_imp->offset];
    }

    template <typename DataType_>
    bool
    ConstElementIterator<storage::Dense, container::Vector, DataType_>::operator< (
            const ConstElementIterator<storage::Dense, container::Vector, DataType_> & other) const
    {
        return this->_imp->index < other._imp->index;
    }

    template <typename DataType_>
    bool
    ConstElementIterator<storage::Dense, container::Vector, DataType_>::operator== (
            const ConstElementIterator<storage::Dense, container::Vector, DataType_> & other) const
    {
        return ((this->_imp->elements.get() == other._imp->elements.get()) && (this->_imp->index == other._imp->index));
    }

    template <typename DataType_>
    bool
    ConstElementIterator<storage::Dense, container::Vector, DataType_>::operator!= (
            const ConstElementIterator<storage::Dense, container::Vector, DataType_> & other) const
    {
        return ((this->_imp->elements.get() != other._imp->elements.get()) || (this->_imp->index != other._imp->index));
    }

    template <typename DataType_>
    unsigned long
    ConstElementIterator<storage::Dense, container::Vector, DataType_>::index() const
    {
        return this->_imp->index;
    }

    template <typename DataType_> struct Implementation<ElementIterator<storage::Dense, container::Vector, DataType_> >
    {
        SharedArray<DataType_> elements;

        unsigned long index;

        unsigned long offset;

        unsigned long stepsize;

        Implementation(const SharedArray<DataType_> & elements, unsigned long index, unsigned long offset, unsigned long stepsize) :
            elements(elements),
            index(index),
            offset(offset),
            stepsize(stepsize)
        {
        }
    };

    template <typename DataType_>
    ElementIterator<storage::Dense, container::Vector, DataType_>::ElementIterator(const SharedArray<DataType_> & elements,
            unsigned long index, unsigned long offset, unsigned long stepsize) :
        PrivateImplementationPattern<ElementIterator<storage::Dense, container::Vector, DataType_>, Single>(
                new Implementation<ElementIterator<storage::Dense, container::Vector, DataType_> >(elements, index, offset, stepsize))
    {
    }

    template <typename DataType_>
    ElementIterator<storage::Dense, container::Vector, DataType_>::ElementIterator(
            const ElementIterator<storage::Dense, container::Vector, DataType_> & other) :
        PrivateImplementationPattern<ElementIterator<storage::Dense, container::Vector, DataType_>, Single>(
                new Implementation<ElementIterator<storage::Dense, container::Vector, DataType_> >(*other._imp))
    {
    }

    template <typename DataType_>
    ElementIterator<storage::Dense, container::Vector, DataType_>::~ElementIterator()
    {
    }

    template <typename DataType_>
    ElementIterator<storage::Dense, container::Vector, DataType_> &
    ElementIterator<storage::Dense, container::Vector, DataType_>::operator= (
            const ElementIterator<storage::Dense, container::Vector, DataType_> & other)
    {
        if (&other == this)
            return *this;

        this->_imp->elements = other._imp->elements;
        this->_imp->index = other._imp->index;
        this->_imp->offset = other._imp->offset;
        this->_imp->stepsize = other._imp->stepsize;

        return *this;
    }

    template <typename DataType_>
    ElementIterator<storage::Dense, container::Vector, DataType_> &
    ElementIterator<storage::Dense, container::Vector, DataType_>::operator++ ()
    {
        CONTEXT("When incrementing ElementIterator<Dense, Vector> by one:");

        ++this->_imp->index;

        return *this;
    }

    template <typename DataType_>
    ElementIterator<storage::Dense, container::Vector, DataType_> &
    ElementIterator<storage::Dense, container::Vector, DataType_>::operator+= (const unsigned long step)
    {
        CONTEXT("When incrementing ElementIterator<Dense, Vector> by '" + stringify(step) + "':");

        this->_imp->index += step;

        return *this;
    }

    template <typename DataType_>
    DataType_ &
    ElementIterator<storage::Dense, container::Vector, DataType_>::operator* () const
    {
        CONTEXT("When accessing assignable element at index '" + stringify(this->_imp->index) + "':");

        return this->_imp->elements[this->_imp->stepsize * this->_imp->index + this->_imp->offset];
    }

    template <typename DataType_>
    bool
    ElementIterator<storage::Dense, container::Vector, DataType_>::operator< (
            const ElementIterator<storage::Dense, container::Vector, DataType_> & other) const
    {
        return this->_imp->index < other._imp->index;
    }

    template <typename DataType_>
    bool
    ElementIterator<storage::Dense, container::Vector, DataType_>::operator== (
            const ElementIterator<storage::Dense, container::Vector, DataType_> & other) const
    {
        return ((this->_imp->elements.get() == other._imp->elements.get()) && (this->_imp->index == other._imp->index));
    }

    template <typename DataType_>
    bool
    ElementIterator<storage::Dense, container::Vector, DataType_>::operator!= (
            const ElementIterator<storage::Dense, container::Vector, DataType_> & other) const
    {
        return ((this->_imp->elements.get() != other._imp->elements.get()) || (this->_imp->index != other._imp->index));
    }

    template <typename DataType_>
    unsigned long
    ElementIterator<storage::Dense, container::Vector, DataType_>::index() const
    {
        return this->_imp->index;
    }

    template <typename DataType_>
    bool
    operator== (const DenseVectorBase<DataType_> & a, const DenseVectorBase<DataType_> & b)
    {
        if (a.size() != b.size())
            throw VectorSizeDoesNotMatch(a.size(), b.size());

        a.lock(lm_read_only);
        b.lock(lm_read_only);
        for (typename DenseVectorBase<DataType_>::ConstElementIterator i(a.begin_elements()), i_end(a.end_elements()),
                j(b.begin_elements()) ; i != i_end ; ++i, ++j)
        {
            if (*i != *i || *j != *j)
            {
                a.unlock(lm_read_only);
                b.unlock(lm_read_only);
                return false;
            }
            if (std::abs(*i - *j) > std::numeric_limits<DataType_>::epsilon())
            {
                a.unlock(lm_read_only);
                b.unlock(lm_read_only);
                return false;
            }
        }
        a.unlock(lm_read_only);
        b.unlock(lm_read_only);

        return true;
    }

    template <>
    bool
    operator== (const DenseVectorBase<unsigned long> & a, const DenseVectorBase<unsigned long> & b)
    {
        if (a.size() != b.size())
            throw VectorSizeDoesNotMatch(a.size(), b.size());

        a.lock(lm_read_only);
        b.lock(lm_read_only);
        for (DenseVectorBase<unsigned long>::ConstElementIterator i(a.begin_elements()), i_end(a.end_elements()),
                j(b.begin_elements()) ; i != i_end ; ++i, ++j)
        {
            if (*i != *i || *j != *j)
            {
                a.unlock(lm_read_only);
                b.unlock(lm_read_only);
                return false;
            }
            if (*i > *j && *i - *j > std::numeric_limits<unsigned long>::epsilon())
            {
                a.unlock(lm_read_only);
                b.unlock(lm_read_only);
                return false;
            }
            if (*j > *i && *j - *i > std::numeric_limits<unsigned long>::epsilon())
            {
                a.unlock(lm_read_only);
                b.unlock(lm_read_only);
                return false;
            }
        }
        a.unlock(lm_read_only);
        b.unlock(lm_read_only);

        return true;
    }

    template <typename DataType_>
    std::ostream &
    operator<< (std::ostream & lhs, const DenseVector<DataType_> & b)
    {
        b.lock(lm_read_only);
        lhs << "[";
        for (typename DenseVector<DataType_>::ConstElementIterator i(b.begin_elements()), i_end(b.end_elements()) ;
                i != i_end ; ++i)
        {
            lhs << "  " << *i;
        }
        lhs << "]" << std::endl;
        b.unlock(lm_read_only);

        return lhs;
    }
}

#endif
