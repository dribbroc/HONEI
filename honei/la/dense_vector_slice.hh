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
#ifndef LIBLA_GUARD_DENSE_VECTOR_SLICE_HH
#define LIBLA_GUARD_DENSE_VECTOR_SLICE_HH 1

#include <honei/la/dense_vector.hh>
#include <honei/util/private_implementation_pattern.hh>
#ifdef HONEI_GMP
#include <gmpxx.h>
#endif

namespace honei
{
    template <typename DataType_> class DenseMatrixTile;

    /**
     * \brief DenseVectorSlice is a vector with O(size) non-zero elements which
     * \brief provides access to a part of a DenseVector for more efficient handling.
     *
     * \ingroup grpvector
     */
    template <typename DataType_> class DenseVectorSlice :
        public DenseVectorBase<DataType_>,
        public PrivateImplementationPattern<DenseVectorSlice<DataType_>, Shared>
    {
        private:
            /**
             * Constructor.
             *
             * For use by DenseMatrixTile.
             *
             * \param size Size of the new slice.
             * \param elements The shared array this slice provides access to.
             * \param offset Offset of the slice's data inside the shared array.
             * \param stepsize Stepsize between two of the vector's elements inside the shared array.
             */
            DenseVectorSlice(const SharedArray<DataType_> & elements, const unsigned long size, const unsigned long offset,
                    const unsigned long stepsize);

        public:
            /// \name Friends of DenseVectorSlice
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

            /// Copy-constructor.
            DenseVectorSlice(const DenseVectorSlice<DataType_> & other);

            /**
             * Constructor.
             *
             * \param source The vector this slice provides access to.
             * \param size Size of the new slice.
             * \param offset Offset of the slice's data inside the source vector's elements.
             * \param stepsize Stepsize between two of the slice's elements inside the vector's elements.
             */
            DenseVectorSlice(const DenseVectorBase<DataType_> & source, const unsigned long size,
                                const unsigned long offset, const unsigned long stepsize);

            DenseVectorSlice(const DenseVectorContinuousBase<DataType_> & source, const unsigned long size,
                                const unsigned long offset, const unsigned long stepsize);


            /// Destructor.
            ~DenseVectorSlice();

            /// \}

            /**
             * \name Functions inherited by Vector
             * \{
             */

            /// Returns const iterator pointing to the first element of the slice.
            virtual ConstElementIterator begin_elements() const;

            /// Returns const iterator pointing behind the last element of the slice.
            virtual ConstElementIterator end_elements() const;

            /// Returns const iterator pointing to a given element of the slice.
            virtual ConstElementIterator element_at(unsigned long index) const;

            /// Returns iterator pointing to the first element of the slice.
            virtual ElementIterator begin_elements();

            /// Returns iterator pointing behind the last element of the slice.
            virtual ElementIterator end_elements();

            /// Returns iterator pointing to a given element of the slice.
            virtual ElementIterator element_at(unsigned long index);

            /// Retrieves element by index, zero-based, assignable.
            virtual const DataType_ & operator[] (unsigned long index) const;

            /// Retrieves element by index, zero-based, assignable.
            virtual DataType_ & operator[] (unsigned long index);

            /// Returns our size.
            virtual unsigned long size() const;

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

            /// Return a copy of the slice.
            DenseVector<DataType_> copy() const;

            /// Returns our offset.
            unsigned long offset() const;

            /// Returns our stepsize.
            unsigned long stepsize() const;
    };

    /**
     * Equality operator for DenseVectorSlice.
     *
     * Compares if corresponding elements of two dense vector slices are equal
     * within machine precision.
     */
    template <typename DataType_> bool operator== (const DenseVectorSlice<DataType_> & a, const DenseVectorSlice<DataType_> & b);

    /**
     * Output operator for DenseVectorSlice.
     *
     * Outputs a dense vector to an output stream.
     */
    template <typename DataType_> std::ostream & operator<< (std::ostream & lhs, const DenseVectorSlice<DataType_> & vector);

    extern template class DenseVectorSlice<float>;

    extern template bool operator== (const DenseVectorSlice<float> & a, const DenseVectorSlice<float> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const DenseVectorSlice<float> & vector);

    extern template class DenseVectorSlice<double>;

    extern template bool operator== (const DenseVectorSlice<double> & a, const DenseVectorSlice<double> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const DenseVectorSlice<double> & vector);

    extern template class DenseVectorSlice<long>;

    extern template bool operator== (const DenseVectorSlice<long> & a, const DenseVectorSlice<long> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const DenseVectorSlice<long> & vector);

    extern template class DenseVectorSlice<unsigned long>;

    //extern template bool operator== (const DenseVectorSlice<unsigned long> & a, const DenseVectorSlice<unsigned long> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const DenseVectorSlice<unsigned long> & vector);

    extern template class DenseVectorSlice<bool>;

    extern template bool operator== (const DenseVectorSlice<bool> & a, const DenseVectorSlice<bool> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const DenseVectorSlice<bool> & vector);

#ifdef HONEI_GMP
    extern template class DenseVectorSlice<mpf_class>;

    extern template bool operator== (const DenseVectorSlice<mpf_class> & a, const DenseVectorSlice<mpf_class> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const DenseVectorSlice<mpf_class> & vector);
#endif
}

#endif
