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

#ifndef LIBLA_GUARD_DENSE_VECTOR_SLICE_HH
#define LIBLA_GUARD_DENSE_VECTOR_SLICE_HH 1

#include <libla/dense_vector.hh>

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
        public DenseVectorBase<DataType_>
    {
        private:
            struct Implementation;

            /// Our implementation.
            std::tr1::shared_ptr<Implementation> _imp;

            /// Our implementation of ElementIteratorBase.
            template <typename ElementType_> class DenseElementIterator;

            typedef typename Vector<DataType_>::VectorElementIterator VectorElementIterator;

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
            friend class DenseElementIterator<DataType_>;
            friend class DenseMatrixTile<DataType_>;

            /// Type of the const iterator over our elements.
            typedef typename Vector<DataType_>::ConstElementIterator ConstElementIterator;

            /// Type of the iterator over our elements.
            typedef typename Vector<DataType_>::ElementIterator ElementIterator;

            /// Constructors
            /// \{

            /**
             * Constructor.
             *
             * \param source The vector this slice provides access to.
             * \param size Size of the new slice.
             * \param offset Offset of the slice's data inside the source vector's elements.
             * \param stepsize Stepsize between two of the slice's elements inside the vector's elements.
             */
            DenseVectorSlice(const DenseVector<DataType_> & source, const unsigned long size,
                                const unsigned long offset, const unsigned long stepsize);

            /// Copy-constructor.
            DenseVectorSlice(const DenseVectorSlice<DataType_> & other);

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

            /// \}

            /**
             * \name Functions inherited by DenseVectorContinousBase
             * \{
             */

            /// Return a pointer to our elements.
            virtual DataType_ * elements() const;

            /// \}

            /// Return a copy of the slice.
            DenseVector<DataType_> copy() const;

            /// Returns our offset.
            unsigned long offset() const;

    };

    extern template class DenseVectorSlice<float>;

    extern template class DenseVectorSlice<double>;

    extern template class DenseVectorSlice<int>;

    extern template class DenseVectorSlice<long>;
}

#endif
