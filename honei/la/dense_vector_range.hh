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

#pragma once
#ifndef LIBLA_GUARD_DENSE_VECTOR_RANGE_HH
#define LIBLA_GUARD_DENSE_VECTOR_RANGE_HH 1

#include <honei/util/tags.hh>
#include <honei/la/dense_vector.hh>
#include <honei/util/private_implementation_pattern.hh>

namespace honei
{
    // Forward declarations.
    template <typename DataType_> class DenseMatrixTile;
    template <typename DataType_> class SparseVector;

    /**
     * \brief DenseVectorRange is a vector which provides access to a contiguous part
     * \brief of a DenseVector and allows to manipulate the accessed part.
     *
     * \ingroup grpvector
     */
    template <typename DataType_> class DenseVectorRange :
        public DenseVectorContinuousBase<DataType_>,
        public PrivateImplementationPattern<DenseVectorRange<DataType_>, Shared>
    {
        private:
            /**
             * Constructor.
             *
             * For use by DenseMatrix and DenseMatrixTile.
             *
             * \param source The shared array the range provides access to.
             * \param size Size of the new dense vector range.
             * \param offset Offset of the vector's data inside the shared array.
             */
            DenseVectorRange(const SharedArray<DataType_> & source, const unsigned long size, const unsigned long offset);

        public:
            /// \name Friends of DenseVector
            /// \{

            friend class DenseMatrix<DataType_>;
            friend class DenseMatrixTile<DataType_>;

            /// \}

            /// Type of the const iterator over our elements.
            typedef honei::ConstElementIterator<storage::Dense, container::Vector, DataType_> ConstElementIterator;

            /// Type of the iterator over our elements.
            typedef honei::ElementIterator<storage::Dense, container::Vector, DataType_> ElementIterator;

            /// Constructors
            /// \{

            /**
             * Constructor.
             *
             * \param source The DenseVector that is accessed by the range.
             * \param size Size of the new dense vector range.
             * \param offset Offset of the vector's data inside the DenseVector's shared array.
             */
            DenseVectorRange(const DenseVector<DataType_> & source, const unsigned long size, const unsigned long offset);

            /// Copy-constructor.
            DenseVectorRange(const DenseVectorRange<DataType_> & other);

            /**
             * Constructor.
             *
             * \param source The DenseVectorRange our range is build of.
             * \param size Size of the new dense vector range.
             * \param offset Offset in the data array.
             */
            DenseVectorRange(const DenseVectorRange<DataType_> & source, const unsigned long size, const unsigned long offset);

            /// Destructor.
            ~DenseVectorRange();

            /// \}

            /**
             * \name Functions inherited by Vector
             * \{
             */

            /// Returns const iterator pointing to the first element of the range.
            virtual ConstElementIterator begin_elements() const;

            /// Returns const iterator pointing behind the last element of the range.
            virtual ConstElementIterator end_elements() const;

            /// Returns const iterator pointing to a given element of the range.
            virtual ConstElementIterator element_at(unsigned long index) const;

            /// Returns iterator pointing to the first element of the range.
            virtual ElementIterator begin_elements();

            /// Returns iterator pointing behind the last element of the range.
            virtual ElementIterator end_elements();

            /// Returns iterator pointing to a given element of the range.
            virtual ElementIterator element_at(unsigned long index);

            /// Returns our size.
            virtual unsigned long size() const;

            /// Retrieves element by index, zero-based, unassignable.
            virtual const DataType_ & operator[] (unsigned long index) const;

            /// Retrieves element by index, zero-based, assignable.
            virtual DataType_ & operator[] (unsigned long index);

            /// \}

            /**
             * \name Functions inherited by DenseVectorContinousBase
             * \{
             */

            /// Return our offset.
            virtual unsigned long offset() const;

            /// Return our stepsize.
            virtual unsigned long stepsize() const;

            /// Return a range of our DenseVectorRange.
            virtual DenseVectorRange<DataType_> range(unsigned long size, unsigned long offset) const;

            /// Return a pointer to our elements.
            virtual DataType_ * elements() const;

            /// Return a reference of the elements array.
            virtual SharedArray<DataType_> & array() const;

            /// Return our memory id
            virtual void * memid() const;

            /// Return the address of our data
            virtual void * address() const;

            /// Request a memory access lock for our data.
            virtual void * lock(LockMode mode, tags::TagValue memory = tags::CPU::memory_value) const;

            /// Release a memory access lock for our data.
            virtual void unlock(LockMode mode) const;

            /// \}

            /// Return a copy of the Vector.
            DenseVector<DataType_> copy() const;
    };

    /**
     * Equality operator for DenseVectorRange.
     *
     * Compares if corresponding elements of two dense vector ranges are equal
     * within machine precision.
     */
    template <typename DataType_> bool operator== (const DenseVectorRange<DataType_> & a, const DenseVectorRange<DataType_> & b);

    /**
     * Output operator for DenseVectorRange.
     *
     * Outputs a dense vector to an output stream.
     */
    template <typename DataType_> std::ostream & operator<< (std::ostream & lhs, const DenseVectorRange<DataType_> & vector);

    extern template class DenseVectorRange<float>;

    extern template bool operator== (const DenseVectorRange<float> & a, const DenseVectorRange<float> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const DenseVectorRange<float> & vector);

    extern template class DenseVectorRange<double>;

    extern template bool operator== (const DenseVectorRange<double> & a, const DenseVectorRange<double> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const DenseVectorRange<double> & vector);

    extern template class DenseVectorRange<int>;

    extern template bool operator== (const DenseVectorRange<int> & a, const DenseVectorRange<int> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const DenseVectorRange<int> & vector);

    extern template class DenseVectorRange<long>;

    extern template bool operator== (const DenseVectorRange<long> & a, const DenseVectorRange<long> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const DenseVectorRange<long> & vector);

    extern template class DenseVectorRange<unsigned long>;

    extern template bool operator== (const DenseVectorRange<unsigned long> & a, const DenseVectorRange<unsigned long> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const DenseVectorRange<unsigned long> & vector);

    extern template class DenseVectorRange<bool>;

    extern template bool operator== (const DenseVectorRange<bool> & a, const DenseVectorRange<bool> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const DenseVectorRange<bool> & vector);

}

#endif
