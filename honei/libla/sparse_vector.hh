/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Michael Abshoff <michael.abshoff@fsmath.mathematik.uni-dortmund.de>
 * Copyright (c) 2007 Sven Mallach <sven.mallach@honei.org>
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

#ifndef LIBLA_GUARD_SPARSE_VECTOR_HH
#define LIBLA_GUARD_SPARSE_VECTOR_HH 1

#include <honei/libla/element_iterator.hh>
#include <honei/libla/vector.hh>
#include <honei/libutil/assertion.hh>
#include <honei/libutil/exception.hh>
#include <honei/libutil/log.hh>
#include <honei/libutil/shared_array-impl.hh>
#include <honei/libutil/stringify.hh>

#include <iterator>
#include <ostream>
#include <string>

namespace honei
{
    /**
     * \brief SparseVector is a vector with O(1) non-zero elements which keeps its data
     * \brief non-sequential.
     *
     * \ingroup grpvector
     */
    template <typename DataType_> class SparseVector :
        public Vector<DataType_>
    {
        private:
            /// Our implementation type.
            class Implementation;

            /// Pointer to our implementation.
            std::tr1::shared_ptr<Implementation> _imp;

            /// Out zero element.
            static const DataType_ _zero_element;

            /**
             * Insert an empty element into the vector and reallocate additional space if necessary.
             *
             * \param position The position at which the new element shall be inserted.
             * \param index The index of the new element.
             */
            void _insert_element(unsigned long position, unsigned long index) const;

            /// Our normal implementation of ElementIteratorBase.
            template <typename ElementType_> class SparseElementIterator;

            /// Our smart implementation of ElementIteratorBase.
            template <typename ElementType_> class NonZeroElementIterator;

            typedef typename Vector<DataType_>::VectorElementIterator VectorElementIterator;

        public:
            friend class SparseElementIterator<DataType_>;
            friend class NonZeroElementIterator<DataType_>;

            /// Type of the const iterator over our elements.
            typedef typename Vector<DataType_>::ConstElementIterator ConstElementIterator;

            /// Type of the iterator over our elements.
            typedef typename Vector<DataType_>::ElementIterator ElementIterator;

            /// \name Constructors and Destructors
            /// \{

            /**
             * Constructor.
             *
             * \param size Size of the new sparse vector.
             * \param capacity Number of elements for which the vector shall reserve space.
             */
            SparseVector(unsigned long size, unsigned long capacity);

            /// Destructor.
            ~SparseVector();

            /// \}

            /// Returns const iterator pointing to the first element of the vector.
            virtual ConstElementIterator begin_elements() const;

            /// Returns iterator pointing behind the last element of the vector.
            virtual ConstElementIterator end_elements() const;

            /// Returns iterator pointing to the first element of the vector.
            virtual ElementIterator begin_elements();

            /// Returns iterator pointing behind the last element of the vector.
            virtual ElementIterator end_elements();

            /// Returns const iterator pointing to the first non-zero element of the vector.
            virtual ConstElementIterator begin_non_zero_elements() const;

            /// Returns const iterator pointing behind the last element of the vector.
            virtual ConstElementIterator end_non_zero_elements() const;

            /// Returns const iterator pointing to a given non zero element of the vector.
            /// Please note that pos != index. pos is the index into the elements array.
            virtual ConstElementIterator non_zero_element_at(unsigned long pos) const;

            /// Returns iterator pointing to the first non-zero element of the vector.
            virtual ElementIterator begin_non_zero_elements();

            /// Returns iterator pointing behind the last element of the vector.
            virtual ElementIterator end_non_zero_elements();

            /// Returns iterator pointing to a given non zero element of the vector.
            /// Please note that pos != index. pos is the index into the elements array.
            virtual ElementIterator non_zero_element_at(unsigned long pos);

            /// Returns out element capacity.
            virtual unsigned long capacity() const;

            /// Returns our used element number.
            virtual unsigned long used_elements() const;

            /// Returns our size.
            virtual unsigned long size() const;

            /// Return a pointer to our elements.
            DataType_ * elements() const;

            /// Return a pointer to our indices.
            unsigned long * indices() const;

            /// Retrieves element by index, zero-based, unassignable.
            virtual const DataType_ & operator[] (unsigned long index) const;

            /// Retrieves (and inserts empty) element by index, zero-based, assignable.
            virtual DataType_ & operator[] (unsigned long index);

            /// Returns a copy of the vector.
            virtual SparseVector copy() const;
    };

    extern template class SparseVector<float>;

    extern template class SparseVector<double>;

    extern template class SparseVector<int>;

    extern template class SparseVector<long>;
}

#endif
