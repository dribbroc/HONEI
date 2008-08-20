/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef LIBLA_GUARD_CONST_VECTOR_IMPL_HH
#define LIBLA_GUARD_CONST_VECTOR_IMPL_HH 1

#include <honei/la/const_vector.hh>
#include <honei/la/vector_error.hh>
#include <honei/util/assertion.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <honei/util/shared_array.hh>
#include <honei/util/stringify.hh>
#include <honei/util/memory_arbiter.hh>

#include <algorithm>
#include <cmath>
#include <string>

namespace honei
{
    // Forward declarations.
    template <typename VectorType_> class VectorImplementation;

    /**
     * \brief ConstVector::Implementation is the private implementation class for ConstVector.
     *
     * \ingroup grpvector
     */
    template <typename DataType_>
    struct Implementation<ConstVector<DataType_> >
    {
        public:
            /// Our size.
            virtual const unsigned long size() const = 0;;

            /// Our offset.
            virtual const unsigned long offset() const = 0;

            /// Return our memory id
            virtual void * memid() const = 0;

            /// Return the address of our data
            virtual void * address() const = 0;

            /// Request a memory access lock for our data.
            virtual void * lock(LockMode mode, tags::TagValue memory = tags::CPU::memory_value) const = 0;

            /// Release a memory access lock for our data.
            virtual void unlock(LockMode mode) const = 0;

            /// Retrieve element at position index, unassignable
            virtual const DataType_ & element(unsigned long index) const = 0;

            virtual ~Implementation()
            {
            }
    };

    template <typename VectorType_ >
    struct VectorImplementation<ConstVector<VectorType_> > :
        public Implementation<ConstVector<typename VectorType_::DataType> >
    {
        private:
            VectorType_ _v;

        public:
            /// Our size.
            virtual const unsigned long size() const
            {
                return _v.size();
            }

            /// Our offset.
            virtual const unsigned long offset() const
            {
                return _v.offset();
            }

            /// Return our memory id
            virtual void * memid() const
            {
                return _v.memid();
            }

            /// Return the address of our data
            virtual void * address() const
            {
                return _v.address();
            }

            /// Request a memory access lock for our data.
            virtual void * lock(LockMode mode, tags::TagValue memory = tags::CPU::memory_value) const
            {
                return _v.lock(mode, memory);
            }

            /// Release a memory access lock for our data.
            virtual void unlock(LockMode mode) const
            {
                _v.unlock(mode);
            }

            /// Constructor
            VectorImplementation(const VectorType_ src) :
                _v(src)
            {
            }

            /// Retrieve element at position index, unassignable
            virtual const typename VectorType_::DataType & element(unsigned long index) const
            {
                return _v[index];
            }
    };

    template <typename DataType_>
    ConstVector<DataType_>::ConstVector(const DenseVector<DataType_> & src) :
        PrivateImplementationPattern<ConstVector<DataType_>, Shared>(new VectorImplementation<ConstVector<DenseVector<DataType_> > >(src))
    {
        CONTEXT("When creating ConstVector from DenseVector:");
    }

    template <typename DataType_>
    ConstVector<DataType_>::ConstVector(const DenseVectorRange<DataType_> & src) :
        PrivateImplementationPattern<ConstVector<DataType_>, Shared>(new VectorImplementation<ConstVector<DenseVectorRange<DataType_> > >(src))
    {
        CONTEXT("When creating ConstVector from DenseVectorRange:");
    }

    template <typename DataType_>
    ConstVector<DataType_>::ConstVector(const DenseVectorSlice<DataType_> & src) :
        PrivateImplementationPattern<ConstVector<DataType_>, Shared>(new VectorImplementation<ConstVector<DenseVectorSlice<DataType_> > >(src))
    {
        CONTEXT("When creating ConstVector from DenseVectorSlice:");
    }

    template <typename DataType_>
    ConstVector<DataType_>::ConstVector(const ConstVector<DataType_> & other) :
        PrivateImplementationPattern<ConstVector<DataType_>, Shared>(other._imp)
    {
    }

    template <typename DataType_>
    ConstVector<DataType_>::~ConstVector()
    {
    }

    template <typename DataType_>
    typename ConstVector<DataType_>::ConstElementIterator
    ConstVector<DataType_>::begin_elements() const
    {
        ConstElementIterator result(*this, 0);
        return result;
    }

    template <typename DataType_>
    typename ConstVector<DataType_>::ConstElementIterator
    ConstVector<DataType_>::end_elements() const
    {
        ConstElementIterator result(*this, this->size());
        return result;
    }

    template <typename DataType_>
    typename ConstVector<DataType_>::ConstElementIterator
    ConstVector<DataType_>::element_at(unsigned long index) const
    {
        ConstElementIterator result(*this, index);
        return result;
    }

    template <typename DataType_>
    unsigned long
    ConstVector<DataType_>::size() const
    {
        return this->_imp->size();
    }

    template <typename DataType_>
    const DataType_ &
    ConstVector<DataType_>::operator[] (unsigned long index) const
    {
        CONTEXT("When retrieving ConstVector element, unassignable:");
        ASSERT(index < this->_imp->size() && index >= 0, "index is out of bounds!");

        return this->_imp->element(index);
    }

    template <typename DataType_>
    unsigned long
    ConstVector<DataType_>::offset() const
    {
        return this->_imp->offset();
    }


    template <typename DataType_>
    void *
    ConstVector<DataType_>::memid() const
    {
        return this->_imp->memid();
    }

    template <typename DataType_>
    void *
    ConstVector<DataType_>::address() const
    {
        return this->_imp->address();
    }

    template <typename DataType_>
    void *
    ConstVector<DataType_>::lock(LockMode mode, tags::TagValue memory) const
    {
        ASSERT(mode == lm_read_only, "Only read_only access is allowed!");
        return this->_imp->lock(mode, memory);
    }

    template <typename DataType_>
    void
    ConstVector<DataType_>::unlock(LockMode mode) const
    {
        ASSERT(mode == lm_read_only, "Only read_only access is allowed!");
        this->_imp->unlock(mode);
    }

    template <typename DataType_>
    bool
    operator== (const ConstVector<DataType_> & a, const ConstVector<DataType_> & b)
    {
        if (a.size() != b.size())
            throw VectorSizeDoesNotMatch(a.size(), b.size());

        for (unsigned long index(0) ; index < a.size() ; ++index)
        {
            if (std::fabs(a[index] - b[index]) <= std::numeric_limits<DataType_>::epsilon())
            {
                continue;
            }
            else return false;
        }
        return true;

    }

    template <typename DataType_>
    std::ostream &
    operator<< (std::ostream & lhs, const ConstVector<DataType_> & vector)
    {
        lhs << "[ ";
        for (unsigned long index(0) ; index < vector.size() ; ++index)
        {
            lhs << vector[index] << " ";
        }
        lhs << "]";

        return lhs;
    }

    template <typename DataType_> struct Implementation<ConstElementIterator<storage::Const, container::Vector, DataType_> >
    {
        const ConstVector<DataType_> & vector;

        unsigned long index;

        Implementation(const ConstVector<DataType_> & vector, unsigned long index) :
            vector(vector),
            index(index)
        {
        }

        Implementation(const Implementation<ElementIterator<storage::Const, container::Vector, DataType_> > & other) :
            vector(other.vector),
            index(other.index)
        {
        }
    };

    template <typename DataType_>
    ConstElementIterator<storage::Const, container::Vector, DataType_>::ConstElementIterator(const ConstVector<DataType_> & vector,
            unsigned long index) :
        PrivateImplementationPattern<ConstElementIterator<storage::Const, container::Vector, DataType_>, Single>(
                new Implementation<ConstElementIterator<storage::Const, container::Vector, DataType_> >(vector, index))
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::Const, container::Vector, DataType_>::ConstElementIterator(const ConstElementIterator & other) :
        PrivateImplementationPattern<ConstElementIterator<storage::Const, container::Vector, DataType_>, Single>(
                new Implementation<ConstElementIterator<storage::Const, container::Vector, DataType_> >(*other._imp))
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::Const, container::Vector, DataType_>::~ConstElementIterator()
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::Const, container::Vector, DataType_> &
    ConstElementIterator<storage::Const, container::Vector, DataType_>::operator++ ()
    {
        CONTEXT("When incrementing ConstElementIterator<Const, Vector> by one:");

        ++this->_imp->index;

        return *this;
    }

    template <typename DataType_>
    ConstElementIterator<storage::Const, container::Vector, DataType_> &
    ConstElementIterator<storage::Const, container::Vector, DataType_>::operator+= (const unsigned long step)
    {
        CONTEXT("When incrementing ConstElementIterator<Const, Vector> by '" + stringify(step) + "':");

        this->_imp->index += step;

        return *this;
    }

    template <typename DataType_>
    const DataType_ &
    ConstElementIterator<storage::Const, container::Vector, DataType_>::operator* () const
    {
        CONTEXT("When accessing unassignable element at index '" + stringify(this->_imp->index) + "':");

        return this->_imp->vector[this->_imp->index];
    }

    template <typename DataType_>
    bool
    ConstElementIterator<storage::Const, container::Vector, DataType_>::operator< (
            const ConstElementIterator<storage::Const, container::Vector, DataType_> & other) const
    {
        return this->_imp->index < other._imp->index;
    }

    template <typename DataType_>
    bool
    ConstElementIterator<storage::Const, container::Vector, DataType_>::operator== (
            const ConstElementIterator<storage::Const, container::Vector, DataType_> & other) const
    {
        return ((&(this->_imp->vector) == &(other._imp->vector)) && (this->_imp->index == other._imp->index));
    }

    template <typename DataType_>
    bool
    ConstElementIterator<storage::Const, container::Vector, DataType_>::operator!= (
            const ConstElementIterator<storage::Const, container::Vector, DataType_> & other) const
    {
        return ((&(this->_imp->vector) != &(other._imp->vector)) || (this->_imp->index != other._imp->index));
    }

    template <typename DataType_>
    unsigned long
    ConstElementIterator<storage::Const, container::Vector, DataType_>::index() const
    {
        return this->_imp->index;
    }
}
#endif
