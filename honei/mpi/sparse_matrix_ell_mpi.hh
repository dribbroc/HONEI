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
            shared_ptr<SparseMatrix<DT_> > _outer;
            std::set<unsigned long> _missing_indices;
            SharedArray<std::set<unsigned long> > _alien_indices;
            unsigned long _rows;
            unsigned long _columns;
            unsigned long _offset;
            unsigned long _rank;
            unsigned long _com_size;
            std::vector<unsigned long> _before_cols;
            std::vector<unsigned long> _middle_cols;
            std::vector<unsigned long> _behind_cols;

        public:

            /// Constructors
            /// \{

            /**
             * Constructor.
             */
            SparseMatrixELLMPI(const SparseMatrix<DT_> & src, const unsigned long rank, const unsigned long com_size) :
                _alien_indices(com_size),
                _rank(rank),
                _com_size(com_size)
            {
                unsigned long part_size(src.rows() / com_size);
                unsigned long rest(src.rows() - (part_size * com_size));
                if (rank < rest)
                    ++part_size;
                _rows = part_size;
                _columns = src.columns();

                unsigned local_offset(0);
                for (unsigned long i(0) ; i < _rank ; ++i)
                {
                    local_offset += src.rows() / com_size;
                    if (i < rest)
                        ++local_offset;
                }
                _offset = local_offset;

                SparseMatrix<DT_> src_part(_rows, src.columns());
                for (unsigned long row(0) ; row < _rows ; ++row)
                {
                    for (unsigned long i(0) ; i < src[row + _offset].used_elements() ; ++i)
                    {
                        src_part(row, (src[row + _offset].indices())[i], (src[row + _offset].elements())[i]);
                    }
                }


                std::vector<unsigned long> inner_row_index;
                std::vector<unsigned long> outer_row_index;

                // matrix fenster ausschneiden und in src_part speichern
                for (unsigned long row(0) ; row < _rows ; ++row)
                {
                    bool outter(false);
                    for (unsigned long i(0) ; i < src_part[row].used_elements() ; ++i)
                    {
                        if ((src_part[row].indices())[i] < _offset || (src_part[row].indices())[i] >= _offset + _rows)
                        {
                            outter = true;
                            break;
                        }
                    }
                    if (outter)
                        outer_row_index.push_back(row);
                    else
                        inner_row_index.push_back(row);
                }

                // inneren teil ohne abhaengigkeiten in inner speichern
                SparseMatrix<DT_> inner(src_part.copy());
                SparseVector<DT_> zeros(_columns, 0);
                for (unsigned long row(0) ; row < outer_row_index.size() ; ++row)
                {
                    inner[outer_row_index.at(row)] = zeros;
                }
                // alle spalten in inner nach links ruecken, damit man mit 0 in einen rechte seite vektor einsteigen kann
                for (unsigned long row(0) ; row < _rows ; ++row)
                {
                    for (unsigned long i(0) ; i < inner[row].used_elements() ; ++i)
                    {
                        (inner[row].indices())[i] -= _offset;
                    }
                }

                // aeusseren teil mit abhaengigkeit nach draussen speichern
                SparseMatrix<DT_> outer(_rows, _columns);
                for (unsigned long row(0) ; row < outer_row_index.size() ; ++row)
                {
                    outer[outer_row_index.at(row)] = src_part[outer_row_index.at(row)];
                }

                _inner.reset(new SparseMatrixELL<DT_>(inner));
                _outer.reset(new SparseMatrix<DT_>(outer));
                _outer->_synch_column_vectors();

                // liste der fehlenden x eintraege erstellen
                for (unsigned long row(0) ; row < _rows ; ++row)
                {
                    for (unsigned long i(0) ; i < outer[row].used_elements() ; ++i)
                    {
                        if ((outer[row].indices())[i] < _offset || (outer[row].indices())[i] >= _offset + _rows)
                            _missing_indices.insert((outer[row].indices())[i]);
                    }
                }

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
                            if (index >= _offset && index < _offset + _rows)
                                _alien_indices[rank].insert(index);
                        }
                    }
                }

                //cols der outer matrix, die nicht leer sind
                for (unsigned long col(0) ; col < _columns && col < _offset; ++col)
                {
                    if (_outer->column(col).used_elements() > 0)
                        _before_cols.push_back(col);
                }
                for (unsigned long col(_offset) ; col < _columns && (col < _offset + _rows)  ; ++col)
                {
                    if (_outer->column(col).used_elements() > 0)
                        _middle_cols.push_back(col);
                }
                for (unsigned long col(_offset + _rows) ; col < _columns; ++col)
                {
                    if (_outer->column(col).used_elements() > 0)
                        _behind_cols.push_back(col);
                }
            }


            /// Copy-constructor.
            SparseMatrixELLMPI(const SparseMatrixELLMPI<DT_> & other) :
                _alien_indices(other._alien_indices),
                _rows(other._rows),
                _columns(other._columns),
                _offset(other._offset),
                _rank(other._rank),
                _com_size(other._com_size)
            {
                // \TODO create a real copy of _alien_indices
                _inner.reset(new SparseMatrixELL<DT_> (*other._inner));
                _outer.reset(new SparseMatrix<DT_> (*other._outer));
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

            const DT_ operator()(unsigned long i, unsigned long j) const
            {
                DT_ result;
                if (j >= _offset)
                    result = (unsigned long)(*_inner)(i, j - _offset) | (unsigned long)(*_outer)(i, j);
                else
                    result = (unsigned long)(*_outer)(i, j);

                return result;
            }

            const SparseMatrixELL<DT_> & inner_matrix() const
            {
                return *_inner;
            }

            SparseMatrixELL<DT_> & inner_matrix()
            {
                return *_inner;
            }

            const SparseMatrix<DT_> & outer_matrix() const
            {
                return *_outer;
            }

            SparseMatrix<DT_> & outer_matrix()
            {
                return *_outer;
            }

            const std::set<unsigned long> & missing_indices() const
            {
                return _missing_indices;
            }

            const std::set<unsigned long> & alien_indices(unsigned long i) const
            {
                return _alien_indices[i];
            }

            const std::vector<unsigned long> & before_cols() const
            {
                return _before_cols;
            }

            const std::vector<unsigned long> & middle_cols() const
            {
                return _middle_cols;
            }

            const std::vector<unsigned long> & behind_cols() const
            {
                return _behind_cols;
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
                result._outer.reset(new SparseMatrix<DT_>(this->_outer->copy()));
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
