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
#include <honei/util/assertion.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <honei/util/shared_array.hh>
#include <honei/util/stringify.hh>
#include <honei/util/memory_arbiter.hh>

#include <algorithm>
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
            virtual unsigned long memid() const = 0;

            /// Return the address of our data
            virtual void * address() const = 0;

            /// Request a read lock for our data.
            virtual void * read(tags::TagValue memory) const = 0;

            /// Release a read lock for our data.
            virtual void release_read() const = 0;

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
            virtual unsigned long memid() const
            {
                return _v.memid();
            }

            /// Return the address of our data
            virtual void * address() const
            {
                return _v.address();
            }

            /// Request a read lock for our data.
            virtual void * read(tags::TagValue memory) const
            {
                return _v.read(memory);
            }

            /// Release a read lock for our data.
            virtual void release_read() const
            {
                _v.release_read();
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
    ConstVector<DataType_>::ConstVector(const ConstVector<DataType_> & other) :
        PrivateImplementationPattern<ConstVector<DataType_>, Shared>(other._imp)
    {
    }

    template <typename DataType_>
    ConstVector<DataType_>::~ConstVector()
    {
    }

    template <typename DataType_>
    typename ConstVector<DataType_>::ConstVectorElementIterator ConstVector<DataType_>::begin_elements() const
    {
        ConstVectorElementIterator result(*this, 0);
        return result;
    }

    template <typename DataType_>
    typename ConstVector<DataType_>::ConstVectorElementIterator ConstVector<DataType_>::end_elements() const
    {
        ConstVectorElementIterator result(*this, this->size());
        return result;
    }

    template <typename DataType_>
    typename ConstVector<DataType_>::ConstVectorElementIterator ConstVector<DataType_>::element_at(unsigned long index) const
    {
        ConstVectorElementIterator result(*this, index);
        return result;
    }

    template <typename DataType_>
    unsigned long ConstVector<DataType_>::size() const
    {
        return this->_imp->size();
    }

    template <typename DataType_>
    const DataType_ & ConstVector<DataType_>::operator[] (unsigned long index) const
    {
        CONTEXT("When retrieving ConstVector element, unassignable:");
        ASSERT(index < this->_imp->size() && index >= 0, "index is out of bounds!");

        return this->_imp->element(index);
    }

    template <typename DataType_>
    unsigned long ConstVector<DataType_>::offset() const
    {
        return this->_imp->offset();
    }


    template <typename DataType_>
    unsigned long ConstVector<DataType_>::memid() const
    {
        return this->_imp->memid();
    }

    template <typename DataType_>
    void * ConstVector<DataType_>::address() const
    {
        return this->_imp->address();
    }

    template <typename DataType_>
    void * ConstVector<DataType_>::read(tags::TagValue memory) const
    {
        return this->_imp->read(memory);
    }

    template <typename DataType_>
    void ConstVector<DataType_>::release_read() const
    {
        this->_imp->release_read();
    }

    template <typename DataType_>
    bool
    operator== (const ConstVector<DataType_> & a, const ConstVector<DataType_> & b)
    {
        if (a.size() != b.size())
            throw VectorSizeDoesNotMatch(a.size(), b.size());

        for (unsigned long index(0) ; index < a.size() ; ++index)
        {
            if (fabs(a[index] - b[index]) <= std::numeric_limits<DataType_>::epsilon())
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
}

#endif
