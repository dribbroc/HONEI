/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2011 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
#ifndef LIBLA_GUARD_DENSE_VECTOR_MPI_HH
#define LIBLA_GUARD_DENSE_VECTOR_MPI_HH 1

#include <honei/util/tags.hh>
#include <honei/la/dense_vector.hh>

#include <iostream>

namespace honei
{
    /**
     * DenseVectorMPI is a vector with O(size) non-zero elements which keeps its data
     * aligned and sequential.
     *
     * \ingroup grpvector
     */
    template <typename DataType_> class DenseVectorMPI
    {
        private:
            shared_ptr<DenseVector<DataType_> > _vector;
            unsigned long _size;
            unsigned long _offset;
            unsigned long _rank;
            unsigned long _com_size;

        public:

            /// Constructors
            /// \{

            /**
             * Constructor.
             *
             * \param size Size of the new dense vector.
             */
            DenseVectorMPI(const DenseVector<DataType_> & src, const unsigned long rank, const unsigned long com_size) :
                _rank(rank),
                _com_size(com_size)
            {
                unsigned long part_size(src.size() / com_size);
                unsigned long rest(src.size() - (part_size * com_size));
                if (rank < rest)
                    ++part_size;
                _size = part_size;

                unsigned local_offset(0);
                for (unsigned long i(0) ; i < _rank ; ++i)
                {
                    local_offset += src.size() / com_size;
                    if (i < rest)
                        ++local_offset;
                }
                _offset = local_offset;

                _vector.reset(new DenseVector<DataType_>(_size));
                for (unsigned long i(0) ; i < _size ; ++i)
                {
                    (*_vector)[i] = src[i + _offset];
                }
            }


            /// Copy-constructor.
            DenseVectorMPI(const DenseVectorMPI<DataType_> & other) :
                _size(other._size),
                _offset(other._offset),
                _rank(other._rank),
                _com_size(other._com_size)
            {
                _vector.reset(new DenseVector<DataType_> (*other._vector));
            }

            /// Destructor.
            virtual ~DenseVectorMPI()
            {
            }

            /// \}

            /// Returns our size.
            virtual unsigned long size() const
            {
                return _size;
            }

            /// Returns our offset into the origin vector.
            virtual unsigned long offset() const
            {
                return _offset;
            }

            /// Retrieves element by index, zero-based, unassignable.
            virtual const DataType_ & operator[] (unsigned long index) const
            {
                return (*_vector)[index];
            }

            /// Retrieves element by index, zero-based, assignable.
            virtual DataType_ & operator[] (unsigned long index)
            {
                return (*_vector)[index];
            }

            /// \{


            /// Return a pointer to our elements.
            virtual DataType_ * elements() const
            {
                return _vector->elements();
            }

            /// Return a reference of the elements array.
            virtual SharedArray<DataType_> & array() const
            {
                return _vector->array();
            }

            /*/// Return our memory id
            virtual void * memid() const;

            /// Return the address of our data
            virtual void * address() const;

            /// Request a memory access lock for our data.
            virtual void * lock(LockMode mode, tags::TagValue memory = tags::CPU::memory_value) const;

            /// Release a memory access lock for our data.
            virtual void unlock(LockMode mode) const;*/

            /// \}

            /// Return a copy of the Vector.
            DenseVectorMPI<DataType_> copy() const
            {
                DenseVectorMPI<DataType_> result(*this);
                result._vector.reset(new DenseVector<DataType_>(this->_vector->copy()));
                return result;
            }
    };

    /**
     * Equality operator for DenseVectorMPI.
     *
     * Compares if corresponding elements of two dense vectors are equal
     * within machine precision.
     */
    template <typename DataType_> bool operator== (const DenseVectorMPI<DataType_> & a, const DenseVectorMPI<DataType_> & b);

    /**
     * Output operator for DenseVectorMPI.
     *
     * Outputs a dense vector to an output stream.
     */
    template <typename DataType_> std::ostream & operator<< (std::ostream & lhs, const DenseVectorMPI<DataType_> & vector);

    /*extern template class DenseVectorMPI<float>;

    extern template bool operator== (const DenseVectorMPI<float> & a, const DenseVectorMPI<float> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const DenseVectorMPI<float> & vector);

    extern template class DenseVectorMPI<double>;

    extern template bool operator== (const DenseVectorMPI<double> & a, const DenseVectorMPI<double> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const DenseVectorMPI<double> & vector);

    extern template class DenseVectorMPI<long>;

    extern template bool operator== (const DenseVectorMPI<long> & a, const DenseVectorMPI<long> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const DenseVectorMPI<long> & vector);

    extern template class DenseVectorMPI<unsigned long>;

    extern template bool operator== (const DenseVectorMPI<unsigned long> & a, const DenseVectorMPI<unsigned long> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const DenseVectorMPI<unsigned long> & vector);*/
}

#endif
