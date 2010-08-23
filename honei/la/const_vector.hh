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

#pragma once
#ifndef LIBLA_GUARD_CONST_VECTOR_HH
#define LIBLA_GUARD_CONST_VECTOR_HH 1

#include <honei/util/tags.hh>
#include <honei/util/private_implementation_pattern.hh>
#include <honei/la/element_iterator.hh>
#include <honei/util/stringify.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/dense_vector_range.hh>
#include <honei/la/dense_vector_slice.hh>

namespace honei
{
    /**
     * ConstVector is a vector with read only element access.
     *
     * \ingroup grpvector
     */
    template <typename DataType_> class ConstVector :
        public PrivateImplementationPattern<ConstVector<DataType_>, Shared>
    {
        public:
            /// Type of the const iterator over our elements.
            typedef honei::ConstElementIterator<storage::Const, container::Vector, DataType_> ConstElementIterator;

            /// Constructors
            /// \{

            /**
             * Constructor.
             *
             * \param src The DenseVector, whose data should be used
             */
            ConstVector(const DenseVector<DataType_> & src);

            /**
             * Constructor.
             *
             * \param src The DenseVectorRange, whose data should be used
             */
            ConstVector(const DenseVectorRange<DataType_> & src);

            /**
             * Constructor.
             *
             * \param src The source DenseVectorSlice, whose data should be used
             */
            ConstVector(const DenseVectorSlice<DataType_> & src);


            /// Copy-constructor.
            ConstVector(const ConstVector<DataType_> & other);

            /// Destructor.
            virtual ~ConstVector();

            /// \}


            /// Returns const iterator pointing to the first element of the vector.
            ConstElementIterator begin_elements() const;

            /// Returns const iterator pointing behind the last element of the vector.
            ConstElementIterator end_elements() const;

            /// Returns const iterator pointing to a given element of the vector.
            ConstElementIterator element_at(unsigned long index) const;

            /// Returns our size.
            unsigned long size() const;

            /// Retrieves element by index, zero-based, unassignable.
            const DataType_ & operator[] (unsigned long index) const;

            /// Return our offset.
            unsigned long offset() const;

            /// Return our memory id
            virtual void * memid() const;

            /// Return the address of our data
            virtual void * address() const;

            /// Request a memory access lock for our data.
            virtual void * lock(LockMode mode, tags::TagValue memory = tags::CPU::memory_value) const;

            /// Release a memory access lock for our data.
            virtual void unlock(LockMode mode) const;

            /// \}
    };

    /**
     * Equality operator for ConstVector.
     *
     * Compares if corresponding elements of two const vectors are equal
     * within machine precision.
     */
    template <typename DataType_> bool operator== (const ConstVector<DataType_> & a, const ConstVector<DataType_> & b);

    /**
     * Output operator for ConstVector.
     *
     * Outputs a ConstVector to an output stream.
     */
    template <typename DataType_> std::ostream & operator<< (std::ostream & lhs, const ConstVector<DataType_> & vector);

    extern template class ConstVector<float>;

    extern template class ConstVector<double>;

    extern template class ConstVector<int>;

    extern template class ConstVector<long>;

    extern template class ConstVector<unsigned long>;
}
#endif
