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
#ifndef MPI_GUARD_SPARSE_MATRIX_ELL_MPI_HH
#define MPI_GUARD_SPARSE_MATRIX_ELL_MPI_HH 1

#include <honei/util/tags.hh>
#include <honei/la/sparse_matrix_ell.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/la/dense_vector.hh>

namespace honei
{
    template <typename DT_> class SparseMatrixELLMPI
    {
        private:
            shared_ptr<SparseMatrixELL<DT_> > _matrix;
            unsigned long _offset;
            unsigned long _rank;
            unsigned long _com_size;

        public:

            /// Constructors
            /// \{

            /**
             * Constructor.
             */
            SparseMatrixELLMPI(const SparseMatrix<DT_> & src, const unsigned long rank, const unsigned long com_size) :
                _rank(rank),
                _com_size(com_size)
            {
                unsigned long part_size(src.rows() / com_size);
                unsigned long rest(src.rows() - (part_size * com_size));
                if (rank < rest)
                    ++part_size;
                unsigned long rows = part_size;

                unsigned local_offset(0);
                for (unsigned long i(0) ; i < _rank ; ++i)
                {
                    local_offset += src.rows() / com_size;
                    if (i < rest)
                        ++local_offset;
                }
                _offset = local_offset;

                SparseMatrix<DT_> src_part(rows, src.columns());
                for (unsigned long row(0) ; row < rows ; ++row)
                {
                    for (unsigned long i(0) ; i < src[row + _offset].used_elements() ; ++i)
                    {
                        src_part(row, (src[row + _offset].indices())[i], (src[row + _offset].elements())[i]);
                    }
                }
                _matrix.reset(new SparseMatrixELL<DT_>(src_part));
            }


            /// Copy-constructor.
            SparseMatrixELLMPI(const SparseMatrixELLMPI<DT_> & other) :
                _offset(other._offset),
                _rank(other._rank),
                _com_size(other._com_size)
            {
                _matrix.reset(new SparseMatrixELL<DT_> (*other._matrix));
            }

            /// Destructor.
            virtual ~SparseMatrixELLMPI()
            {
            }

            /// \}

            /// Returns our size.
            virtual unsigned long size() const
            {
                return _matrix->size();
            }

            /// Returns our row count.
            virtual unsigned long rows() const
            {
                return _matrix->rows();
            }

            /// Returns our column count.
            virtual unsigned long columns() const
            {
                return _matrix->columns();
            }

            /// Returns our offset into the origin vector.
            virtual unsigned long offset() const
            {
                return _offset;
            }

            const DT_ operator()(unsigned long i, unsigned long j) const
            {
                return (*_matrix)(i, j);
            }

            const SparseMatrixELL<DT_> & matrix() const
            {
                return *_matrix;
            }

            SparseMatrixELL<DT_> & matrix()
            {
                return *_matrix;
            }

            /// \{


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
            SparseMatrixELLMPI<DT_> copy() const
            {
                SparseMatrixELLMPI<DT_> result(*this);
                result._matrix.reset(new SparseMatrixELL<DT_>(this->_matrix->copy()));
                return result;
            }
    };

    /**
     * Equality operator for SparseMatrixELLMPI.
     *
     * Compares if corresponding elements of two dense vectors are equal
     * within machine precision.
     */
    template <typename DT_> bool operator== (const SparseMatrixELLMPI<DT_> & a, const SparseMatrixELLMPI<DT_> & b);

    /**
     * Output operator for SparseMatrixELLMPI.
     *
     * Outputs a dense vector to an output stream.
     */
    template <typename DT_> std::ostream & operator<< (std::ostream & lhs, const SparseMatrixELLMPI<DT_> & vector);

    /*extern template class SparseMatrixELLMPI<float>;

    extern template bool operator== (const SparseMatrixELLMPI<float> & a, const SparseMatrixELLMPI<float> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const SparseMatrixELLMPI<float> & vector);

    extern template class SparseMatrixELLMPI<double>;

    extern template bool operator== (const SparseMatrixELLMPI<double> & a, const SparseMatrixELLMPI<double> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const SparseMatrixELLMPI<double> & vector);

    extern template class SparseMatrixELLMPI<long>;

    extern template bool operator== (const SparseMatrixELLMPI<long> & a, const SparseMatrixELLMPI<long> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const SparseMatrixELLMPI<long> & vector);

    extern template class SparseMatrixELLMPI<unsigned long>;

    extern template bool operator== (const SparseMatrixELLMPI<unsigned long> & a, const SparseMatrixELLMPI<unsigned long> & b);

    extern template std::ostream & operator<< (std::ostream & lhs, const SparseMatrixELLMPI<unsigned long> & vector);*/
}

#endif
