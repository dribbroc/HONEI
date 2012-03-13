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
#ifndef MPI_GUARD_DENSE_VECTOR_MPI_HH
#define MPI_GUARD_DENSE_VECTOR_MPI_HH 1

#include <honei/util/tags.hh>
#include <honei/util/configuration.hh>
#include <honei/la/dense_vector.hh>
#include <honei/backends/mpi/operations.hh>

namespace honei
{
    template <typename DT_> class DenseVectorMPI
    {
        private:
            shared_ptr<DenseVector<DT_> > _vector;
            unsigned long _orig_size;
            unsigned long _offset;
            unsigned long _rank;
            unsigned long _com_size;

        public:

            static unsigned long calc_size(unsigned long _orig_size, MPI_Comm com = MPI_COMM_WORLD)
            {
                unsigned long size;
                unsigned long _rank(mpi::mpi_comm_rank(com));
                unsigned long _com_size(mpi::mpi_comm_size(com));

                unsigned long min_part_size(Configuration::instance()->get_value("mpi::min_part_size", 1));
                unsigned long pre_part_size(_orig_size / _com_size);
                if (pre_part_size >= min_part_size)
                {
                    unsigned long part_size(pre_part_size);
                    unsigned long rest(_orig_size - (part_size * _com_size));
                    if (_rank < rest)
                        ++part_size;
                    return part_size;
                }
                else
                {
                    unsigned long count(_orig_size / min_part_size);
                    if (count == 0)
                        count = 1;
                    unsigned long part_size(_orig_size / count);
                    unsigned long rest(_orig_size - (part_size * count));
                    if (_rank < rest)
                        ++part_size;
                    size = part_size;
                    if (_rank >= count)
                        size = 0;
                    return size;
                }
            }

            static unsigned long calc_offset(unsigned long _orig_size, MPI_Comm com = MPI_COMM_WORLD)
            {
                unsigned long size;
                unsigned long _rank(mpi::mpi_comm_rank(com));
                unsigned long _com_size(mpi::mpi_comm_size(com));

                unsigned long min_part_size(Configuration::instance()->get_value("mpi::min_part_size", 1));
                unsigned long pre_part_size(_orig_size / _com_size);
                if (pre_part_size >= min_part_size)
                {
                    unsigned long part_size(pre_part_size);
                    unsigned long rest(_orig_size - (part_size * _com_size));
                    if (_rank < rest)
                        ++part_size;
                    size = part_size;

                    unsigned local_offset(0);
                    for (unsigned long i(0) ; i < _rank ; ++i)
                    {
                        local_offset += _orig_size / _com_size;
                        if (i < rest)
                            ++local_offset;
                    }
                    return local_offset;
                }
                else
                {
                    unsigned long count(_orig_size / min_part_size);
                    if (count == 0)
                        count = 1;
                    unsigned long part_size(_orig_size / count);
                    unsigned long rest(_orig_size - (part_size * count));
                    if (_rank < rest)
                        ++part_size;
                    size = part_size;
                    if (_rank >= count)
                        size = 0;

                    unsigned local_offset(0);
                    for (unsigned long i(0) ; i < _rank ; ++i)
                    {
                        local_offset += _orig_size / count;
                        if (i < rest)
                            ++local_offset;
                    }
                    return local_offset;
                }
            }

            /// Constructors
            /// \{

            /**
             * Constructor.
             *
             * \param src The src for the new dense vector.
             */
            DenseVectorMPI(const DenseVector<DT_> & src,  MPI_Comm com = MPI_COMM_WORLD) :
                _orig_size(src.size()),
                _rank(mpi::mpi_comm_rank(com)),
                _com_size(mpi::mpi_comm_size(com))
            {
                unsigned long size = calc_size(_orig_size, com);
                _offset = calc_offset(_orig_size, com);

                _vector.reset(new DenseVector<DT_>(size));
                for (unsigned long i(0) ; i < size ; ++i)
                {
                    (*_vector)[i] = src[i + _offset];
                }
            }

            /**
             * Constructor.
             */
            DenseVectorMPI(const DenseVector<DT_> & src,  unsigned long global_size, MPI_Comm com = MPI_COMM_WORLD) :
                _orig_size(global_size),
                _rank(mpi::mpi_comm_rank(com)),
                _com_size(mpi::mpi_comm_size(com))
            {
                _offset = calc_offset(_orig_size, com);

                _vector.reset(new DenseVector<DT_>(src.copy()));
            }

            /**
             * Constructor.
             *
             * \param size Size of the new dense vector.
             */
            DenseVectorMPI(unsigned long src_size,  MPI_Comm com = MPI_COMM_WORLD) :
                _orig_size(src_size),
                _rank(mpi::mpi_comm_rank(com)),
                _com_size(mpi::mpi_comm_size(com))
        {
                unsigned long size = calc_size(_orig_size, com);
                _offset = calc_offset(_orig_size, com);

            _vector.reset(new DenseVector<DT_>(size));
        }

            /**
             * Constructor.
             *
             * \param size Size of the new dense vector.
             * \param value Value the vector will be filled with.
             */
            DenseVectorMPI(unsigned long src_size,  DT_ value, MPI_Comm com = MPI_COMM_WORLD) :
                _orig_size(src_size),
                _rank(mpi::mpi_comm_rank(com)),
                _com_size(mpi::mpi_comm_size(com))
        {
                unsigned long size = calc_size(_orig_size, com);
                _offset = calc_offset(_orig_size, com);

            _vector.reset(new DenseVector<DT_>(size, value));
        }


            /// Copy-constructor.
            DenseVectorMPI(const DenseVectorMPI<DT_> & other) :
                _orig_size(other._orig_size),
                _offset(other._offset),
                _rank(other._rank),
                _com_size(other._com_size)
        {
            _vector.reset(new DenseVector<DT_> (*other._vector));
            //_vector = std::tr1::shared_ptr<DenseVector<DT_> >(new DenseVector<DT_>(*other._vector));
        }

            /// Destructor.
            virtual ~DenseVectorMPI()
            {
            }

            /// \}

            /// Returns our size.
            virtual unsigned long local_size() const
            {
                return _vector->size();
            }

            /// Returns our original size.
            virtual unsigned long size() const
            {
                return _orig_size;
            }

            /// Returns our offset into the origin vector.
            virtual unsigned long offset() const
            {
                return _offset;
            }

            /// Retrieves element by index, zero-based, unassignable.
            virtual const DT_ & operator[] (unsigned long index) const
            {
                return (*_vector)[index];
            }

            /// Retrieves element by index, zero-based, assignable.
            virtual DT_ & operator[] (unsigned long index)
            {
                return (*_vector)[index];
            }

            const DenseVector<DT_> & vector() const
            {
                return *_vector;
            }

            DenseVector<DT_> & vector()
            {
                return *_vector;
            }

            /// \{


            /// Return a pointer to our elements.
            virtual DT_ * elements() const
            {
                return _vector->elements();
            }

            /// Return a reference of the elements array.
            virtual SharedArray<DT_> & array() const
            {
                return _vector->array();
            }

            /*/// Return our memory id
              virtual void * memid() const;

            /// Return the address of our data
            virtual void * address() const;*/

            /// Request a memory access lock for our data.
            virtual void * lock(LockMode mode, tags::TagValue memory = tags::CPU::memory_value) const
            {
                return _vector->lock(mode, memory);
            }

            /// Release a memory access lock for our data.
            virtual void unlock(LockMode mode) const
            {
                return _vector->unlock(mode);
            }

            /// \}

            /// Return a copy of the Vector.
            DenseVectorMPI<DT_> copy() const
            {
                DenseVectorMPI<DT_> result(*this);
                result._vector.reset(new DenseVector<DT_>(this->_vector->copy()));
                //result._vector = std::tr1::shared_ptr<DenseVector<DT_> >(new DenseVector<DT_>(this->_vector->copy()));
                return result;
            }
    };

    /**
     * Equality operator for DenseVectorMPI.
     *
     * Compares if corresponding elements of two dense vectors are equal
     * within machine precision.
     */
    template <typename DT_> bool operator== (const DenseVectorMPI<DT_> & a, const DenseVectorMPI<DT_> & b)
    {
            if (a.size() != b.size()) return false;
            if (a.offset() != b.offset()) return false;
            return a.vector()==b.vector();
    }

    /**
     * Output operator for DenseVectorMPI.
     *
     * Outputs a dense vector to an output stream.
     */
    template <typename DT_> std::ostream & operator<< (std::ostream & lhs, const DenseVectorMPI<DT_> & vector)
    {
        lhs << vector;
        return lhs;
    }

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
