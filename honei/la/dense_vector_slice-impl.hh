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

#ifndef LIBLA_GUARD_DENSE_VECTOR_SLICE_IMPL_HH
#define LIBLA_GUARD_DENSE_VECTOR_SLICE_IMPL_HH 1

#include <honei/la/dense_matrix_tile.hh>
#include <honei/la/dense_vector-impl.hh>
#include <honei/la/dense_vector_slice-impl.hh>
#include <honei/la/element_iterator.hh>
#include <honei/la/vector.hh>
#include <honei/util/assertion.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <honei/util/shared_array-impl.hh>
#include <honei/util/stringify.hh>

#include <iterator>

namespace honei
{
    template <typename DataType_> class DenseMatrixTile;

    /// Implementation for DenseVectorSlice template.
    template <typename DataType_> struct Implementation<DenseVectorSlice<DataType_> >
    {
        /// Pointer to our elements.
        SharedArray<DataType_> elements;

        /// Our size.
        const unsigned long size;

        /// Our offset.
        const unsigned long offset;

        /// Our stepsize.
        const unsigned long stepsize;

        /// Constructor.
        Implementation(const SharedArray<DataType_> & e, const unsigned long s, const unsigned long o, const unsigned long ss) :
            elements(e),
            size(s),
            offset(o),
            stepsize(ss)
        {
            CONTEXT("When creating DenseVectorSlice:");
            ASSERT(stepsize > 0, "stepsize is zero!");
            ASSERT(size > 0, "size is zero!");
            ASSERT(elements.size() > offset + (size - 1) * stepsize, "end of slice is beyond end of source!");
        }
    };

    /// Private Constructor
    template <typename DataType_>
    DenseVectorSlice<DataType_>::DenseVectorSlice(const SharedArray<DataType_> & elements, const unsigned long size,
            const unsigned long offset, const unsigned long stepsize) :
        PrivateImplementationPattern<DenseVectorSlice<DataType_>, Shared>(new Implementation<DenseVectorSlice<DataType_> >(
                    elements, size, offset, stepsize))
    {
    }

    /// Constructor for creation of a new DenseVectorSlice from some DenseVectorBase implementation type
    template <typename DataType_>
    DenseVectorSlice<DataType_>::DenseVectorSlice(const DenseVectorBase<DataType_> & source, const unsigned long size,
            const unsigned long offset, const unsigned long stepsize) :
        PrivateImplementationPattern<DenseVectorSlice<DataType_>, Shared>(new Implementation<DenseVectorSlice<DataType_> >(
                    source.array(), size, offset, stepsize))
    {
    }

    /// Constructor for creation of a new DenseVectorSlice from some DenseVectorContinuousBase implementation type
    template <typename DataType_>
    DenseVectorSlice<DataType_>::DenseVectorSlice(const DenseVectorContinuousBase<DataType_> & source, const unsigned long size,
            const unsigned long offset, const unsigned long stepsize) :
        PrivateImplementationPattern<DenseVectorSlice<DataType_>, Shared>(new Implementation<DenseVectorSlice<DataType_> >(
                    source.array(), size, offset, stepsize))
    {
    }

    /// Copy-Constructor
    template <typename DataType_>
    DenseVectorSlice<DataType_>::DenseVectorSlice(const DenseVectorSlice<DataType_> & other) :
        PrivateImplementationPattern<DenseVectorSlice<DataType_>, Shared>(new Implementation<DenseVectorSlice<DataType_> >(
                    other._imp->elements, other._imp->size, other._imp->offset, other._imp->stepsize))
    {
    }

    template <typename DataType_>
    DenseVectorSlice<DataType_>::~DenseVectorSlice()
    {
    }

    template <typename DataType_>
    typename DenseVectorSlice<DataType_>::ConstElementIterator
    DenseVectorSlice<DataType_>::begin_elements() const
    {
        return ConstElementIterator(this->_imp->elements, 0, this->_imp->offset, this->_imp->stepsize);
    }

    template <typename DataType_>
    typename DenseVectorSlice<DataType_>::ConstElementIterator
    DenseVectorSlice<DataType_>::end_elements() const
    {
        return ConstElementIterator(this->_imp->elements, this->_imp->size, this->_imp->offset, this->_imp->stepsize);
    }

    template <typename DataType_>
    typename DenseVectorSlice<DataType_>::ConstElementIterator
    DenseVectorSlice<DataType_>::element_at(unsigned long index) const
    {
        return ConstElementIterator(this->_imp->elements, index, this->_imp->offset, this->_imp->stepsize);
    }

    template <typename DataType_>
    typename DenseVectorSlice<DataType_>::ElementIterator
    DenseVectorSlice<DataType_>::begin_elements()
    {
        return ElementIterator(this->_imp->elements, 0, this->_imp->offset, this->_imp->stepsize);
    }

    template <typename DataType_>
    typename DenseVectorSlice<DataType_>::ElementIterator
    DenseVectorSlice<DataType_>::end_elements()
    {
        return ElementIterator(this->_imp->elements, this->_imp->size, this->_imp->offset, this->_imp->stepsize);
    }

    template <typename DataType_>
    typename DenseVectorSlice<DataType_>::ElementIterator
    DenseVectorSlice<DataType_>::element_at(unsigned long index)
    {
        return ElementIterator(this->_imp->elements, index, this->_imp->offset, this->_imp->stepsize);
    }

    template <typename DataType_>
    unsigned long DenseVectorSlice<DataType_>::stepsize() const
    {
        return this->_imp->stepsize;
    }

    template <typename DataType_>
    unsigned long DenseVectorSlice<DataType_>::size() const
    {
        return this->_imp->size;
    }

    template <typename DataType_>
    unsigned long DenseVectorSlice<DataType_>::offset() const
    {
        return this->_imp->offset;
    }

    template <typename DataType_>
    const DataType_ & DenseVectorSlice<DataType_>::operator[] (unsigned long index) const
    {
        CONTEXT("When retrieving DenseVectorSlice element, unassignable:");
        ASSERT(index < this->_imp->size && index >= 0, "index is out of bounds!");
        return this->_imp->elements[this->_imp->stepsize * index + this->_imp->offset];
    }

    template <typename DataType_>
    DataType_ & DenseVectorSlice<DataType_>::operator[] (unsigned long index)
    {
        CONTEXT("When retrieving DenseVectorSlice element, assignable:");
        ASSERT(index < this->_imp->size && index >= 0, "index is out of bounds!");
        return this->_imp->elements[this->_imp->stepsize * index + this->_imp->offset];
    }

    template <typename DataType_>
    inline DataType_ * DenseVectorSlice<DataType_>::elements() const
    {
        return this->_imp->elements.get() + this->_imp->offset;
    }

    template <typename DataType_>
    inline SharedArray<DataType_> & DenseVectorSlice<DataType_>::array() const
    {
        return this->_imp->elements;
    }

    template <typename DataType_>
    void * DenseVectorSlice<DataType_>::memid() const
    {
        return this->_imp->elements.get();
    }

    template <typename DataType_>
    inline void * DenseVectorSlice<DataType_>::address() const
    {
        return this->_imp->elements.get();
    }

    template <typename DataType_>
    void * DenseVectorSlice<DataType_>::lock(LockMode mode, tags::TagValue memory) const
    {
        return MemoryArbiter::instance()->lock(mode, memory, this->memid(), this->address(), this->size() * sizeof(DataType_));
    }

    template <typename DataType_>
    void DenseVectorSlice<DataType_>::unlock(LockMode mode) const
    {
        MemoryArbiter::instance()->unlock(mode, this->memid());
    }

    template <typename DataType_>
    DenseVector<DataType_> DenseVectorSlice<DataType_>::copy() const
    {
        DenseVector<DataType_> result(this->_imp->size);
        DataType_ * source(this->_imp->elements.get());
        DataType_ * target(result.elements());
        for (unsigned long i(0) ; i < this->_imp->size ; i++)
        {
            target[i] = source[this->_imp->stepsize * i + this->_imp->offset];
        }

        ///\todo: Use TypeTraits.

        return result;
    }

    template <typename DataType_>
    bool
    operator== (const DenseVectorSlice<DataType_> & a, const DenseVectorSlice<DataType_> & b)
    {
        if (a.size() != b.size())
            throw VectorSizeDoesNotMatch(a.size(), b.size());

        for (typename DenseVectorSlice<DataType_>::ConstElementIterator i(a.begin_elements()), i_end(a.end_elements()),
                j(b.begin_elements()) ; i != i_end ; ++i, ++j)
        {
            if (std::fabs(*i - *j) > std::numeric_limits<DataType_>::epsilon())
                return false;
        }

        return true;
    }

    template <typename DataType_>
    std::ostream &
    operator<< (std::ostream & lhs, const DenseVectorSlice<DataType_> & b)
    {
        lhs << "[";
        for (typename DenseVectorSlice<DataType_>::ConstElementIterator i(b.begin_elements()), i_end(b.end_elements()) ;
                i != i_end ; ++i)
        {
            lhs << "  " << *i;
        }
        lhs << "]" << std::endl;

        return lhs;
    }
}

#endif
