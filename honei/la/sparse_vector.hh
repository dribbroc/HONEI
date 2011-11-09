/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Michael Abshoff <michael.abshoff@fsmath.mathematik.uni-dortmund.de>
 * Copyright (c) 2007 Sven Mallach <sven.mallach@honei.org>
 * Copyright (c) 2009 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
#ifndef LIBLA_GUARD_SPARSE_VECTOR_HH
#define LIBLA_GUARD_SPARSE_VECTOR_HH 1

#include <honei/la/element_iterator.hh>
#include <honei/la/sparse_matrix-fwd.hh>
#include <honei/util/assertion.hh>
#include <honei/util/exception.hh>
#include <honei/util/log.hh>
#include <honei/util/private_implementation_pattern.hh>
#include <honei/util/shared_array.hh>
#include <honei/util/stringify.hh>


namespace honei
{
    /**
     * \brief SparseVector is a vector with O(1) non-zero elements which keeps its data
     * \brief non-sequential.
     *
     * \ingroup grpvector
     */
    template <typename DataType_> class SparseVector :
        public PrivateImplementationPattern<SparseVector<DataType_>, Shared>
    {
        private:
            /// Out zero element.
            static const DataType_ _zero_element;

            /**
             * Insert an empty element into the vector and reallocate additional space if necessary.
             *
             * \param position The position at which the new element shall be inserted.
             * \param index The index of the new element.
             */
            void _insert_element(unsigned long position, unsigned long index) const;

        public:
            /// \name Friends of SparseVector
            /// \{

            friend class SparseMatrix<DataType_>;
            friend class honei::ConstElementIterator<storage::Sparse, container::Vector, DataType_>;
            friend class honei::ElementIterator<storage::Sparse, container::Vector, DataType_>;
            friend class honei::ElementIterator<storage::SparseNonZero, container::Vector, DataType_>;
            friend class honei::ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_>;

            /// \}

            /// Type of the const iterator over our elements.
            typedef honei::ConstElementIterator<storage::Sparse, container::Vector, DataType_> ConstElementIterator;

            /// Type of the iterator over our elements.
            typedef honei::ElementIterator<storage::Sparse, container::Vector, DataType_> ElementIterator;

            /// Type of the const iterator over our non zero elements.
            typedef honei::ConstElementIterator<storage::SparseNonZero, container::Vector, DataType_> NonZeroConstElementIterator;

            /// Type of the iterator over our non zero elements.
            typedef honei::ElementIterator<storage::SparseNonZero, container::Vector, DataType_> NonZeroElementIterator;


            /// \name Constructors and Destructors
            /// \{

            /**
             * Constructor.
             *
             * \param size Size of the new sparse vector.
             * \param capacity Number of elements for which the vector shall reserve space.
             */
            SparseVector(unsigned long size, unsigned long capacity = 1);

            /// Copy-constructor.
            SparseVector(const SparseVector<DataType_> & other);

            /// Destructor.
            ~SparseVector();

            /// \}

            /// Returns const iterator pointing to the first element of the vector.
            ConstElementIterator begin_elements() const;

            /// Returns iterator pointing behind the last element of the vector.
            ConstElementIterator end_elements() const;

            /// Returns iterator pointing to the first element of the vector.
            ElementIterator begin_elements();

            /// Returns iterator pointing behind the last element of the vector.
            ElementIterator end_elements();

            /// Returns const iterator pointing to the first non-zero element of the vector.
            NonZeroConstElementIterator begin_non_zero_elements() const;

            /// Returns const iterator pointing behind the last element of the vector.
            NonZeroConstElementIterator end_non_zero_elements() const;

            /// Returns const iterator pointing to a given non zero element of the vector.
            /// Please note that pos != index. pos is the index into the elements array.
            NonZeroConstElementIterator non_zero_element_at(unsigned long pos) const;

            /// Returns iterator pointing to the first non-zero element of the vector.
            NonZeroElementIterator begin_non_zero_elements();

            /// Returns iterator pointing behind the last element of the vector.
            NonZeroElementIterator end_non_zero_elements();

            /// Returns iterator pointing to a given non zero element of the vector.
            /// Please note that pos != index. pos is the index into the elements array.
            NonZeroElementIterator non_zero_element_at(unsigned long pos);

            /// Returns out element capacity.
            unsigned long capacity() const;

            /// Returns our used element number.
            unsigned long used_elements() const;

            /// Returns our size.
            unsigned long size() const;

            /// Return a pointer to our elements.
            DataType_ * elements() const;

            /// Return a pointer to our indices.
            unsigned long * indices() const;

            /// Retrieves element by index, zero-based, unassignable.
            const DataType_ & operator[] (unsigned long index) const;

            /// Retrieves (and inserts empty) element by index, zero-based, assignable.
            DataType_ & operator[] (unsigned long index);

            /// Returns a copy of the vector.
            SparseVector copy() const;
    };

    /**
     * Equality operator for SparseVector.
     *
     * Compares if corresponding elements of two sparse vectors are equal
     * within machine precision.
     */
    template <typename DataType_> bool operator== (const SparseVector<DataType_> & a, const SparseVector<DataType_> & b);

    /**
     * Output operator for SparseVector.
     *
     * Outputs a sparse vector to an output stream.
     */
    template <typename DataType_> std::ostream & operator<< (std::ostream & lhs, const SparseVector<DataType_> & vector);

    extern template class SparseVector<float>;

    extern template bool operator== (const SparseVector<float> & a, const SparseVector<float> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const SparseVector<float> & vector);

    extern template class SparseVector<double>;

    extern template bool operator== (const SparseVector<double> & a, const SparseVector<double> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const SparseVector<double> & vector);

    extern template class SparseVector<long>;

    extern template bool operator== (const SparseVector<long> & a, const SparseVector<long> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const SparseVector<long> & vector);

    extern template class SparseVector<unsigned long>;

    extern template bool operator== (const SparseVector<unsigned long> & a, const SparseVector<unsigned long> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const SparseVector<unsigned long> & vector);

}

#endif
