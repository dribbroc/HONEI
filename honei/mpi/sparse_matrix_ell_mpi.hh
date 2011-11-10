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
#include <honei/backends/mpi/operations.hh>

#include <vector>
#include <set>

namespace honei
{
    template <typename DT_> class SparseMatrixELLMPI
    {
        private:
            shared_ptr<SparseMatrixELL<DT_> > _inner;
            shared_ptr<SparseMatrixELL<DT_> > _outer;
            std::set<unsigned long> _missing_indices;
            SharedArray<std::set<unsigned long> > _alien_indices;
            unsigned long _rows;
            unsigned long _columns;
            unsigned long _offset;
            unsigned long _x_offset;
            unsigned long _col_part_size;
            unsigned long _rank;
            unsigned long _com_size;

        public:

            /// Constructors
            /// \{

            /**
             * Constructor.
             */
            SparseMatrixELLMPI(const SparseMatrix<DT_> & src, MPI_Comm com = MPI_COMM_WORLD) :
                _alien_indices(mpi::mpi_comm_size(com))
            {
                int irank;
                mpi::mpi_comm_rank(&irank, com);
                _rank = irank;
                int icom_size;
                mpi::mpi_comm_size(&icom_size, com);
                _com_size = icom_size;

                unsigned long part_size(src.rows() / _com_size);
                unsigned long rest(src.rows() - (part_size * _com_size));
                if (_rank < rest)
                    ++part_size;
                _rows = part_size;
                _columns = src.columns();

                unsigned local_offset(0);
                for (unsigned long i(0) ; i < _rank ; ++i)
                {
                    local_offset += src.rows() / _com_size;
                    if (i < rest)
                        ++local_offset;
                }
                _offset = local_offset;

                unsigned long col_part_size(src.columns() / _com_size);
                unsigned long col_rest(src.columns() - (col_part_size * _com_size));
                if (_rank < col_rest)
                    ++col_part_size;
                _col_part_size = col_part_size;

                unsigned local_x_offset(0);
                for (unsigned long i(0) ; i < _rank ; ++i)
                {
                    local_x_offset += src.columns() / _com_size;
                    if (i < col_rest)
                        ++local_x_offset;
                }
                _x_offset = local_x_offset;

                // matrix fenster ausschneiden und in src_part speichern
                SparseMatrix<DT_> src_part(_rows, _columns);
                for (unsigned long row(0) ; row < _rows ; ++row)
                {
                    for (unsigned long i(0) ; i < src[row + _offset].used_elements() ; ++i)
                    {
                        src_part(row, (src[row + _offset].indices())[i], (src[row + _offset].elements())[i]);
                    }
                }


                // inneren teil ohne abhaengigkeiten in inner speichern
                SparseMatrix<DT_> inner(_rows, _col_part_size);
                for (unsigned long col(0) ; col < _col_part_size ; ++col)
                {
                    for (unsigned long i(0) ; i < src_part.column(col + _x_offset).used_elements() ; ++i)
                        inner((src_part.column(col + _x_offset).indices())[i], col, (src_part.column(col + _x_offset).elements())[i]);
                }


                // extern abhaengige spalten berechnen
                for (unsigned long col(0) ; col < _x_offset ; ++col)
                {
                    if (src_part.column(col).used_elements() != 0)
                            _missing_indices.insert(col);
                }
                for (unsigned long col(_x_offset + _col_part_size) ; col < _columns  ; ++col)
                {
                    if (src_part.column(col).used_elements() != 0)
                            _missing_indices.insert(col);
                }

                // outer matrix mit allen abhaengigkeiten nach draussen erstellen
                SparseMatrix<DT_> outer_comp(_rows, _missing_indices.size());
                {
                    unsigned long cix(0);
                    for (std::set<unsigned long>::iterator ci(_missing_indices.begin()) ; ci != _missing_indices.end() ; ++ci, ++cix)
                    {
                        for (unsigned long i(0) ; i < src_part.column(*ci).used_elements() ; ++i)
                        {
                            outer_comp((src_part.column(*ci).indices())[i], cix, (src_part.column(*ci).elements())[i]);
                        }
                    }
                }

                _inner.reset(new SparseMatrixELL<DT_>(inner));
                _outer.reset(new SparseMatrixELL<DT_>(outer_comp));

                // liste an alle anderen prozesse schicken
                for (unsigned long rank(0) ; rank < _com_size ; ++rank)
                {
                    if (rank == _rank)
                    {
                        unsigned long count(_missing_indices.size());
                        mpi::mpi_bcast(&count, 1, rank);

                        for (std::set<unsigned long>::iterator i(_missing_indices.begin()) ; i != _missing_indices.end() ; ++i)
                        {
                            unsigned long index(*i);
                            mpi::mpi_bcast(&index, 1, rank);
                        }
                    }
                    else
                    {
                        unsigned long count;
                        mpi::mpi_bcast(&count, 1, rank);
                        for (unsigned long i(0) ; i < count ; ++i)
                        {
                            unsigned long index;
                            mpi::mpi_bcast(&index, 1, rank);
                            if (index >= _x_offset && index < _x_offset + col_part_size)
                                _alien_indices[rank].insert(index);
                        }
                    }
                }
            }


            /// Copy-constructor.
            SparseMatrixELLMPI(const SparseMatrixELLMPI<DT_> & other) :
                _alien_indices(other._alien_indices),
                _rows(other._rows),
                _columns(other._columns),
                _offset(other._offset),
                _x_offset(other._x_offset),
                _col_part_size(other._col_part_size),
                _rank(other._rank),
                _com_size(other._com_size)
            {
                // \TODO create a real copy of _alien_indices
                _inner.reset(new SparseMatrixELL<DT_> (*other._inner));
                _outer.reset(new SparseMatrixELL<DT_> (*other._outer));
            }

            /// Destructor.
            virtual ~SparseMatrixELLMPI()
            {
            }

            /// \}

            /// Returns our size.
            virtual unsigned long size() const
            {
                return _rows * _columns;
            }

            /// Returns our row count.
            virtual unsigned long rows() const
            {
                return _rows;
            }

            /// Returns our column count.
            virtual unsigned long columns() const
            {
                return _columns;
            }

            /// Returns our offset into the origin vector.
            virtual unsigned long offset() const
            {
                return _offset;
            }

            /// Returns our offset into the origin vector.
            virtual unsigned long x_offset() const
            {
                return _x_offset;
            }

            const DT_ operator()(unsigned long i, unsigned long j) const
            {
                if (j >= _x_offset && j < _x_offset + _col_part_size)
                {
                    return (*_inner)(i, j - _x_offset);
                }
                else
                {
                    unsigned long cix(0);
                    for (std::set<unsigned long>::iterator ci(_missing_indices.begin()) ; ci != _missing_indices.end() ; ++ci, ++cix)
                    {
                        if (*ci == j)
                            return (*_outer)(i, cix);
                    }
                    return DT_(0);
                }
            }

            const SparseMatrixELL<DT_> & inner_matrix() const
            {
                return *_inner;
            }

            /*SparseMatrixELL<DT_> & inner_matrix()
            {
                return *_inner;
            }*/

            const SparseMatrixELL<DT_> & outer_matrix() const
            {
                return *_outer;
            }

            /*SparseMatrixELL<DT_> & outer_matrix()
            {
                return *_outer;
            }*/

            const std::set<unsigned long> & missing_indices() const
            {
                return _missing_indices;
            }

            const std::set<unsigned long> & alien_indices(unsigned long i) const
            {
                return _alien_indices[i];
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
                result._inner.reset(new SparseMatrixELL<DT_>(this->_inner->copy()));
                result._outer.reset(new SparseMatrixELL<DT_>(this->_outer->copy()));
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
