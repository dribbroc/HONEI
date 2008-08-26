/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#ifndef HONEI_GUARD_LA_DENSE_VECTOR_BASE_HH
#define HONEI_GUARD_LA_DENSE_VECTOR_BASE_HH 1

#include <honei/la/element_iterator.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/util/private_implementation_pattern.hh>

namespace honei
{
    /**
     * DenseVectorBase is the abstract baseclass for all dense vector types.
     *
     * \ingroup grpvector
     */
    template <typename DataType_> class DenseVectorBase
    {
        public:
            /// Type of the const iterator over our elements.
            typedef ConstElementIterator<storage::Dense, container::Vector, DataType_> ConstElementIterator;

            /// Returns const iterator pointing to the first element of the vector.
            virtual ConstElementIterator begin_elements() const = 0;

            /// Returns const iterator pointing behind the last element of the vector.
            virtual ConstElementIterator end_elements() const = 0;

            /// Returns const iterator pointing to a given element of the vector.
            virtual ConstElementIterator element_at(unsigned long index) const = 0;

            /// Type of the iterator over our elements.
            typedef ElementIterator<storage::Dense, container::Vector, DataType_> ElementIterator;

            /// Returns iterator pointing to the first element of the vector.
            virtual ElementIterator begin_elements() = 0;

            /// Returns iterator pointing behind the last element of the vector.
            virtual ElementIterator end_elements() = 0;

            /// Returns iterator pointing to a given element of the vector.
            virtual ElementIterator element_at(unsigned long index) = 0;

            /// Returns our size.
            virtual unsigned long size() const = 0;

            /// Return a pointer to our elements.
            virtual DataType_ * elements() const = 0;

            /// Retrieves element by index, zero-based, unassignable.
            virtual const DataType_ & operator[] (unsigned long index) const = 0;

            /// Retrieves element by index, zero-based, assignable
            virtual DataType_ & operator[] (unsigned long index) = 0;

            /// Type of our elements.
            typedef DataType_ DataType;

            /// Destructor.
            virtual ~DenseVectorBase() {}

            /// Return our memory id
            virtual void * memid() const = 0;

            /// Return the address of our data
            virtual void * address() const = 0;

            /// Request a memory access lock for our data.
            virtual void * lock(LockMode mode, tags::TagValue memory = tags::CPU::memory_value) const = 0;

            /// Release a memory access lock for our data.
            virtual void unlock(LockMode mode) const = 0;
    };

    /**
     * DenseVectorContinousBase is the abstract base class for all dense vectors which
     * keep their data continous in memory.
     *
     * \ingroup grpvector
     */
    template <typename DataType_> class DenseVectorContinuousBase :
        public DenseVectorBase<DataType_>
    {
        public:
            /// Return our offset.
            virtual unsigned long offset() const = 0;

            /// Return a range of our DenseVectorContinuousBase.
            virtual DenseVectorRange<DataType_> range(unsigned long size, unsigned long offset) const = 0;
    };
}

#endif
