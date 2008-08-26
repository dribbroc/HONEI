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

#include <honei/la/dense_matrix-fwd.hh>
#include <honei/la/dense_vector_base.hh>
#include <honei/la/dense_vector_range-fwd.hh>
#include <honei/la/dense_vector_slice-fwd.hh>
#include <honei/util/tags.hh>
#include <honei/util/private_implementation_pattern.hh>

namespace honei
{
    // Forward declarations
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
            /// \name Friends of DenseVector
            /// \{

            friend class DenseMatrix<DataType_>;
            friend class DenseVectorRange<DataType_>;
            friend class DenseVectorSlice<DataType_>;

            /// \}

            /// Type of the const iterator over our elements.
            typedef ConstElementIterator<storage::Dense, container::Vector, DataType_> ConstElementIterator;

            /// Type of the iterator over our elements.
            typedef ElementIterator<storage::Dense, container::Vector, DataType_> ElementIterator;

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
             * \param src The source vector, whose data should be used
             * \param size Size of the new dense vector.
             */
            DenseVector(DenseVector<DataType_> & src, unsigned long size);

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

            /// Retrieves element by index, zero-based, unassignable.
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
     * Equality operator for DenseVector.
     *
     * Compares if corresponding elements of two banded matrices are equal
     * within machine precision.
     */
    template <typename DataType_> bool operator== (const DenseVectorBase<DataType_> & a, const DenseVectorBase<DataType_> & b);

    /**
     * Output operator for DenseVector.
     *
     * Outputs a dense vector to an output stream.
     */
    template <typename DataType_> std::ostream & operator<< (std::ostream & lhs, const DenseVector<DataType_> & vector);

    extern template class DenseVector<float>;

    extern template bool operator== (const DenseVectorBase<float> & a, const DenseVectorBase<float> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const DenseVector<float> & vector);

    extern template class DenseVector<double>;

    extern template bool operator== (const DenseVectorBase<double> & a, const DenseVectorBase<double> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const DenseVector<double> & vector);

    extern template class DenseVector<int>;

    extern template bool operator== (const DenseVectorBase<int> & a, const DenseVectorBase<int> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const DenseVector<int> & vector);

    extern template class DenseVector<long>;

    extern template bool operator== (const DenseVectorBase<long> & a, const DenseVectorBase<long> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const DenseVector<long> & vector);

    extern template class DenseVector<unsigned long>;

    extern template bool operator== (const DenseVectorBase<unsigned long> & a, const DenseVectorBase<unsigned long> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const DenseVector<unsigned long> & vector);
}

#endif
