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

#ifndef LIBLA_GUARD_DENSE_VECTOR_RANGE_HH
#define LIBLA_GUARD_DENSE_VECTOR_RANGE_HH 1

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
            /// Our implementation of ElementIteratorBase.
            class DenseElementIterator;

            typedef typename Vector<DataType_>::VectorElementIterator VectorElementIterator;

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
            friend class DenseElementIterator;
            friend class DenseMatrix<DataType_>;
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

            /// Return a range of our DenseVectorRange.
            virtual DenseVectorRange<DataType_> range(unsigned long size, unsigned long offset) const;

            /// Return a pointer to our elements.
            virtual DataType_ * elements() const;

            /// Return our memory id
            virtual unsigned long memid() const;

            /// \}

            /// Return a copy of the Vector.
            DenseVector<DataType_> copy() const;
    };

    extern template class DenseVectorRange<float>;

    extern template class DenseVectorRange<double>;

    extern template class DenseVectorRange<int>;

    extern template class DenseVectorRange<long>;

    extern template class DenseVectorRange<unsigned long>;
}

#endif
