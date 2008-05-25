/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Michael Abshoff <michael.abshoff@fsmath.mathematik.uni-dortmund.de>
 * Copyright (c) 2007 Sven Mallach <sven.mallach@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_DENSE_VECTOR_HH
#define LIBLA_GUARD_DENSE_VECTOR_HH 1

#include <honei/la/vector.hh>
#include <honei/util/private_implementation_pattern.hh>

namespace honei
{
    // Forward declarations.
    template <typename DataType_> class DenseMatrix;
    template <typename DataType_> class DenseVectorRange;
    template <typename DataType_> class DenseVectorSlice;
    template <typename DataType_> class SparseVector;

    /**
     * DenseVector is a vector with O(size) non-zero elements which keeps its data
     * aligned and sequential.
     *
     * \ingroup grpvector
     */
    template <typename DataType_> class DenseVector :
        public DenseVectorContinuousBase<DataType_>,
        public PrivateImplementationPattern<DenseVector<DataType_>, Shared>
    {
        private:
            /// Our implementation of ElementIteratorBase.
            class DenseElementIterator;

            typedef typename Vector<DataType_>::VectorElementIterator VectorElementIterator;

            /**
             * Constructor.
             *
             * For use by DenseMatrix.
             *
             * \param size Size of the new dense vector.
             * \param elements SharedArray of the vector's elements.
             * \param offset Offset of the vector's data inside the shared array.
             * \param stepsize Stepsize between two of the vector's elements inside the shared array.
             */
            DenseVector(const unsigned long size, const SharedArray<DataType_> & elements, unsigned long offset = 0,
                    unsigned stepsize = 1);

        public:
            friend class DenseElementIterator;
            friend class DenseMatrix<DataType_>;
            friend class DenseVectorRange<DataType_>;
            friend class DenseVectorSlice<DataType_>;

            /// Type of the const iterator over our elements.
            typedef typename Vector<DataType_>::ConstElementIterator ConstElementIterator;

            /// Type of the iterator over our elements.
            typedef typename Vector<DataType_>::ElementIterator ElementIterator;

            /// Constructors
            /// \{

            /**
             * Constructor.
             *
             * \param size Size of the new dense vector.
             * \param offset Offset of the vector's data inside the shared array.
             * \param stepsize Stepsize between two of the vector's elements inside the shared array.
             */
            DenseVector(const unsigned long size, const unsigned long offset = 0, const unsigned long stepsize = 1);

            /**
             * Constructor.
             *
             * \param size Size of the new dense vector.
             * \param value Default value for all of the vector's elements.
             * \param offset Offset of the vector's data inside the shared array.
             * \param stepsize Stepsize between two of the vector's elements inside the shared array.
             */
            DenseVector(const unsigned long size, DataType_ value, unsigned long offset = 0,
                    unsigned long stepsize = 1);

            /**
             * Constructor.
             *
             * \param other The SparseVector to densify.
             */
            DenseVector(const SparseVector<DataType_> & other);

            /// Copy-constructor.
            DenseVector(const DenseVector<DataType_> & other);

            /// Destructor.
            ~DenseVector();

            /// \}

            /**
             * \name Functions inherited by Vector
             * \{
             */

            /// Returns const iterator pointing to the first element of the vector.
            virtual ConstElementIterator begin_elements() const;

            /// Returns const iterator pointing behind the last element of the vector.
            virtual ConstElementIterator end_elements() const;

            /// Returns const iterator pointing to a given element of the vector.
            virtual ConstElementIterator element_at(unsigned long index) const;

            /// Returns iterator pointing to the first element of the vector.
            virtual ElementIterator begin_elements();

            /// Returns iterator pointing behind the last element of the vector.
            virtual ElementIterator end_elements();

            /// Returns iterator pointing to a given element of the vector.
            virtual ElementIterator element_at(unsigned long index);

            /// Returns our size.
            virtual unsigned long size() const;

            /// Retrieves element by index, zero-based, assignable.
            virtual const DataType_ & operator[] (unsigned long index) const;

            /// Retrieves element by index, zero-based, assignable.
            virtual DataType_ & operator[] (unsigned long index);

            /// \}

            /// \{

            /**
             * Arrow operator
             *
             * Used as hack for BandIterator.
             */

            DenseVector<DataType_> * operator-> ()
            {
                return this;
            }

            const DenseVector<DataType_> * const operator-> () const
            {
                return this;
            }

            /// \}

            /**
             * \name Functions inherited by DenseVectorContinousBase
             * \{
             */

            /// Return our offset.
            virtual unsigned long offset() const;

            /// Return a range of our DenseVector.
            virtual DenseVectorRange<DataType_> range(unsigned long size, unsigned long offset) const;

            /// Return a pointer to our elements.
            virtual DataType_ * elements() const;

            /// \}

            /// Return a copy of the Vector.
            DenseVector<DataType_> copy() const;
    };

    extern template class DenseVector<float>;

    extern template class DenseVector<double>;

    extern template class DenseVector<int>;

    extern template class DenseVector<long>;

    extern template class DenseVector<unsigned long>;
}

#endif
