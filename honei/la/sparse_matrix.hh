/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Sven Mallach <sven.mallach@honei.org>
 * Copyright (c) 2009, 2010 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
#ifndef LIBLA_GUARD_SPARSE_MATRIX_HH
#define LIBLA_GUARD_SPARSE_MATRIX_HH 1

#include <honei/la/element_iterator.hh>
#include <honei/la/sparse_vector.hh>
#include <honei/la/sparse_matrix_ell.hh>
#include <honei/la/sparse_matrix_csr.hh>
#include <honei/la/banded_matrix_qx.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/matrix_error.hh>
#include <honei/util/shared_array-impl.hh>
#include <limits>

namespace honei
{
    /**
     * \brief SparseMatrix is a matrix with O(rows) sparse row-vectors.
     *
     * \ingroup grpmatrix
     **/
    template <typename DataType_> class SparseMatrix
    {
        private:
            /// Our row-vectors' initial capacity.
            unsigned long _capacity;

            /// Our columns.
            unsigned long _columns;

            /// Our rows.
            unsigned long _rows;

            /// Our row-vectors.
            SharedArray<shared_ptr<SparseVector<DataType_> > > _row_vectors;

            /// Our column-vectors.
            SharedArray<shared_ptr<SparseVector<DataType_> > > _column_vectors;

            /// Our zero-vector.
            SparseVector<DataType_> _zero_vector;

            /// Our zero-column-vector.
            SparseVector<DataType_> _zero_column_vector;

        public:
            friend class honei::ConstElementIterator<storage::Sparse, container::Matrix, DataType_>;
            friend class honei::ElementIterator<storage::Sparse, container::Matrix, DataType_>;
            friend class honei::ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_>;
            friend class honei::ElementIterator<storage::SparseNonZero, container::Matrix, DataType_>;
            friend struct honei::Implementation<honei::ConstElementIterator<storage::Sparse, container::Matrix, DataType_> >;
            friend struct honei::Implementation<honei::ElementIterator<storage::Sparse, container::Matrix, DataType_> >;
            friend struct honei::Implementation<honei::ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_> >;
            friend struct honei::Implementation<honei::ElementIterator<storage::SparseNonZero, container::Matrix, DataType_> >;


            /// Type of the const iterator over our elements.
            typedef honei::ConstElementIterator<storage::Sparse, container::Matrix, DataType_> ConstElementIterator;

            /// Type of the iterator over our elements.
            typedef honei::ElementIterator<storage::Sparse, container::Matrix, DataType_> ElementIterator;

            /// Type of the const iterator over our non zero elements.
            typedef honei::ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_> NonZeroConstElementIterator;

            /// Type of the iterator over our non zero elements.
            typedef honei::ElementIterator<storage::SparseNonZero, container::Matrix, DataType_> NonZeroElementIterator;

            void _synch_column_vectors();


            /// \name Constructors
            /// \{

            /**
             * Constructor.
             *
             * \param columns Number of columns of the new matrix.
             * \param rows Number of rows of the new matrix.
             * \param capacity Initial capacity for non-zero elements.
             **/
            SparseMatrix(unsigned long rows, unsigned long columns, unsigned long capacity = 1) :
                _capacity(capacity),
                _columns(columns),
                _rows(rows),
                _row_vectors(rows + 1),
                _column_vectors(columns + 1),
                _zero_vector(columns, 1),
                _zero_column_vector(rows, 1)
            {
                CONTEXT("When creating SparseMatrix:");
                ASSERT(rows > 0, "number of rows is zero!");
                ASSERT(columns > 0, "number of columns is zero!");
                ASSERT(columns >= capacity, "capacity '" + stringify(capacity) + "' exceeds row-vector size '" +
                        stringify(columns) + "'!");

                _row_vectors[rows].reset(new SparseVector<DataType_>(columns, 1));
                _column_vectors[columns].reset(new SparseVector<DataType_>(rows, 1));
            }

            explicit SparseMatrix(DenseMatrix<DataType_> & src, unsigned long capacity = 1) :
                _capacity(capacity),
                _columns(src.columns()),
                _rows(src.rows()),
                _row_vectors(src.rows() + 1),
                _column_vectors(src.columns() + 1),
                _zero_vector(src.columns(), 1),
                _zero_column_vector(src.rows(), 1)
            {
                CONTEXT("When creating SparseMatrix:");
                ASSERT(src.columns() >= capacity, "capacity '" + stringify(capacity) + "' exceeds row-vector size '" +
                        stringify(src.columns()) + "'!");

                _row_vectors[src.rows()].reset(new SparseVector<DataType_>(src.columns(), 1));
                _column_vectors[src.columns()].reset(new SparseVector<DataType_>(src.rows(), 1));

                typename DenseMatrix<DataType_>::ElementIterator i(src.begin_elements());
                typename SparseMatrix<DataType_>::ElementIterator si(this->begin_elements());
                while (i < src.end_elements())
                {
                    /// \todo Removed hardcoded zero element.
                    if (*i != DataType_(0))
                        *si = *i;
                    ++i;
                    ++si;
                }

                _synch_column_vectors();
            }

            explicit SparseMatrix(const SparseMatrixELL<DataType_> & src, unsigned long capacity = 1) :
                _capacity(capacity),
                _columns(src.columns()),
                _rows(src.rows()),
                _row_vectors(src.rows() + 1),
                _column_vectors(src.columns() + 1),
                _zero_vector(src.columns(), 1),
                _zero_column_vector(src.rows(), 1)
            {
                CONTEXT("When creating SparseMatrix from SparseMatrixELL:");
                ASSERT(src.columns() >= capacity, "capacity '" + stringify(capacity) + "' exceeds row-vector size '" +
                        stringify(src.columns()) + "'!");

                _row_vectors[src.rows()].reset(new SparseVector<DataType_>(src.columns(), 1));
                _column_vectors[src.columns()].reset(new SparseVector<DataType_>(src.rows(), 1));

                for (unsigned long row(0) ; row < src.rows() ; ++row)
                {
                    for (unsigned long i(row*src.threads()) ; i < src.Aj().size() ; i += src.stride())
                    {
                        for (unsigned long thread(0) ; thread < src.threads() ; ++thread)
                        {
                            if (src.Ax()[i + thread] != 0)
                                (*this)(row, src.Aj()[i + thread]) = src.Ax()[i + thread];
                        }
                    }
                }

                _synch_column_vectors();
            }

            explicit SparseMatrix(const SparseMatrixCSR<DataType_> & src, unsigned long capacity = 1) :
                _capacity(capacity),
                _columns(src.columns()),
                _rows(src.rows()),
                _row_vectors(src.rows() + 1),
                _column_vectors(src.columns() + 1),
                _zero_vector(src.columns(), 1),
                _zero_column_vector(src.rows(), 1)
            {
                CONTEXT("When creating SparseMatrix from SparseMatrixCSR:");
                ASSERT(src.columns() >= capacity, "capacity '" + stringify(capacity) + "' exceeds row-vector size '" +
                        stringify(src.columns()) + "'!");

                _row_vectors[src.rows()].reset(new SparseVector<DataType_>(src.columns(), 1));
                _column_vectors[src.columns()].reset(new SparseVector<DataType_>(src.rows(), 1));

                for (unsigned long row(0) ; row < src.rows() ; ++row)
                {
                    for (unsigned long i(src.Ar()[row]) ; i < src.Ar()[row + 1] ; ++i)
                    {
                        (*this)(row, src.Aj()[i]) = src.Ax()[i * src.blocksize()];
                        for (unsigned long blocki(1) ; blocki < src.blocksize() ; ++blocki)
                        {
                            if(src.Ax()[i + blocki] != DataType_(0))
                            {
                                (*this)(row, src.Aj()[i]+blocki) = src.Ax()[i * src.blocksize() + blocki];
                            }
                            else
                            {
                                break;
                            }
                        }
                    }
                }

                _synch_column_vectors();
            }

            explicit SparseMatrix(const BandedMatrixQx<Q1Type, DataType_> & src, unsigned long capacity = 1) :
                _capacity(capacity),
                _columns(src.columns()),
                _rows(src.rows()),
                _row_vectors(src.rows() + 1),
                _column_vectors(src.columns() + 1),
                _zero_vector(src.columns(), 1),
                _zero_column_vector(src.rows(), 1)
        {
            CONTEXT("When creating SparseMatrix from BandedMatrixQ1:");
            ASSERT(src.columns() >= capacity, "capacity '" + stringify(capacity) + "' exceeds row-vector size '" +
                    stringify(src.columns()) + "'!");

            _row_vectors[src.rows()].reset(new SparseVector<DataType_>(src.columns(), 1));
            _column_vectors[src.columns()].reset(new SparseVector<DataType_>(src.rows(), 1));

            /*for (unsigned long row(0) ; row < src.rows() ; ++row)
              {
              for (unsigned long column(0) ; column < src.columns() ; ++column)
              {
              if (src(row, column) != DataType_(0))
              (*this)(row, column) = src(row, column);
              }
              }*/
            long root(src.root());

            DenseVector<DataType_> band_ll(src.band(LL));
            for (unsigned long i(root + 1) ; i < band_ll.size() ; ++i)
            {
                if (band_ll.elements()[i] != DataType_(0))
                    (*this)(i, -root - 1 + i) = band_ll.elements()[i];
            }
            DenseVector<DataType_> band_ld(src.band(LD));
            for (unsigned long i(root) ; i < band_ld.size() ; ++i)
            {
                if (band_ld.elements()[i] != DataType_(0))
                    (*this)(i, -root + i) = band_ld.elements()[i];
            }
            DenseVector<DataType_> band_lu(src.band(LU));
            for (unsigned long i(root -1) ; i < band_lu.size() ; ++i)
            {
                if (band_lu.elements()[i] != DataType_(0))
                    (*this)(i, -root + 1 + i) = band_lu.elements()[i];
            }
            DenseVector<DataType_> band_dl(src.band(DL));
            for (unsigned long i(1) ; i < band_dl.size() ; ++i)
            {
                if (band_dl.elements()[i] != DataType_(0))
                    (*this)(i, - 1 + i) = band_dl.elements()[i];
            }
            DenseVector<DataType_> band_dd(src.band(DD));
            for (unsigned long i(0) ; i < band_dd.size() ; ++i)
            {
                if (band_dd.elements()[i] != DataType_(0))
                    (*this)(i, i) = band_dd.elements()[i];
            }
            DenseVector<DataType_> band_du(src.band(DU));
            for (unsigned long i(0) ; i < band_du.size() - 1; ++i)
            {
                if (band_du.elements()[i] != DataType_(0))
                    (*this)(i, 1 + i) = band_du.elements()[i];
            }
            DenseVector<DataType_> band_ul(src.band(UL));
            for (unsigned long i(0) ; i < band_ul.size() - root + 1; ++i)
            {
                if (band_ul.elements()[i] != DataType_(0))
                    (*this)(i, root - 1 + i) = band_ul.elements()[i];
            }
            DenseVector<DataType_> band_ud(src.band(UD));
            for (unsigned long i(0) ; i < band_ud.size() - root ; ++i)
            {
                if (band_ud.elements()[i] != DataType_(0))
                    (*this)(i, root + i) = band_ud.elements()[i];
            }
            DenseVector<DataType_> band_uu(src.band(UU));
            for (unsigned long i(0) ; i < band_uu.size() - root + -1; ++i)
            {
                if (band_uu.elements()[i] != DataType_(0))
                    (*this)(i, root + 1 + i) = band_uu.elements()[i];
            }

            _synch_column_vectors();
        }

            explicit SparseMatrix(unsigned long rows, unsigned long columns, unsigned long * row_indices,
                    unsigned long * column_indices, DataType_ * data, unsigned long nnz) :
                _capacity(1),
                _columns(columns),
                _rows(rows),
                _row_vectors(rows + 1),
                _column_vectors(columns + 1),
                _zero_vector(columns, 1),
                _zero_column_vector(rows, 1)
            {
                CONTEXT("When creating SparseMatrix:");
                ASSERT(rows > 0, "number of rows is zero!");
                ASSERT(columns > 0, "number of columns is zero!");
                ASSERT(columns >= _capacity, "capacity '" + stringify(_capacity) + "' exceeds row-vector size '" +
                        stringify(columns) + "'!");

                _row_vectors[rows].reset(new SparseVector<DataType_>(columns, 1));
                _column_vectors[columns].reset(new SparseVector<DataType_>(rows, 1));

                for (unsigned long i(0) ; i < nnz ; ++i)
                {
                    (*this)(row_indices[i], column_indices[i]) = data[i];
                }

                _synch_column_vectors();
            }

            ~SparseMatrix()
            {
                MemoryPool<tags::CPU>::instance()->release_free();
            }

            /// \}

            /// Returns iterator pointing to the first element of the matrix.
            ConstElementIterator begin_elements() const
            {
                return ConstElementIterator(*this, 0);
            }

            /// Returns iterator pointing behind the last element of the matrix.
            ConstElementIterator end_elements() const
            {
                return ConstElementIterator(*this, _rows * _columns);
            }

            /// Returns iterator pointing to the first element of the matrix.
            ElementIterator begin_elements()
            {
                return ElementIterator(*this, 0);
            }

            /// Returns iterator pointing behind the last element of the matrix.
            ElementIterator end_elements()
            {
                return ElementIterator(*this, _rows * _columns);
            }

            /// Returns const iterator pointing to the first non-zero element of the matrix.
            NonZeroConstElementIterator begin_non_zero_elements() const
            {
                return NonZeroConstElementIterator(*this, 0);
            }

            /// Returns const iterator pointing behind the last element of the matrix.
            NonZeroConstElementIterator end_non_zero_elements() const
            {
                return NonZeroConstElementIterator(*this, 1 /* Dummy */);
            }

            /// Returns iterator pointing to the first non-zero element of the matrix.
            NonZeroElementIterator begin_non_zero_elements()
            {
                return NonZeroElementIterator(*this, 0);
            }

            /// Returns iterator pointing behind the last element of the matrix.
            NonZeroElementIterator end_non_zero_elements()
            {
                return NonZeroElementIterator(*this, 1 /* Dummy */);
            }

            /// Returns the number of our columns.
            unsigned long columns() const
            {
                return _columns;
            }

            /// Returns the number of our rows.
            unsigned long rows() const
            {
                return _rows;
            }

            /// Returns our size = rows * columns.
            unsigned long size() const
            {
                return _rows * _columns;
            }

            /// Return number of non zero elements.
            unsigned long used_elements() const
            {
                unsigned long ue(0);
                for (unsigned long i(0) ; i < rows() ; ++i)
                    if (_row_vectors[i])
                        ue += _row_vectors[i]->used_elements();
                return ue;
            }

            /// Retrieves row vector by index, zero-based, unassignable.
            const SparseVector<DataType_> & operator[] (unsigned long row) const
            {
                if (! _row_vectors[row])
                    return _zero_vector;

                return *_row_vectors[row];
            }

            /// Retrieves row vector by index, zero-based, assignable.
            SparseVector<DataType_> & operator[] (unsigned long row)
            {
                if (! _row_vectors[row])
                    _row_vectors[row].reset(new SparseVector<DataType_>(_columns, _capacity));

                return *_row_vectors[row];
            }

            /// Retrieves column vector by index, zero-based, unassignable.
            /// Only usable with operator(row, column, value) !
            const SparseVector<DataType_> & column(unsigned long column) const
            {
                if (! _column_vectors[column])
                    return _zero_column_vector;

                return *_column_vectors[column];
            }

            /// Retrieves column vector by index, zero-based, assignable.
            /// Only usable with operator(row, column, value) !
            SparseVector<DataType_> & column(unsigned long column)
            {
                if (! _column_vectors[column])
                    _column_vectors[column].reset(new SparseVector<DataType_>(_rows, _capacity));

                return *_column_vectors[column];
            }

            /// Retrieves element at (row, column), unassignable.
            const DataType_ & operator() (unsigned long row, unsigned long column) const
            {
                if (! _row_vectors[row])
                    /// \todo access element directly
                    return _zero_vector[0];

                /// \todo access element directly
                return ((const SparseVector<DataType_>)*_row_vectors[row])[column];
            }

            /// Retrieves element at (row, column), assignable.
            DataType_ & operator() (unsigned long row, unsigned long column)
            {
                if (! _row_vectors[row])
                    _row_vectors[row].reset(new SparseVector<DataType_>(_columns, _capacity));

                /// \todo access element directly
                return (*_row_vectors[row])[column];
            }

            /// Writes element at (row, column).
            void operator() (unsigned long row, unsigned long column, DataType_ value)
            {
                if (! _row_vectors[row])
                    _row_vectors[row].reset(new SparseVector<DataType_>(_columns, _capacity));

                if (! _column_vectors[column])
                    _column_vectors[column].reset(new SparseVector<DataType_>(_rows, _capacity));

                /// \todo access element directly
                (*_row_vectors[row])[column] = value;
                (*_column_vectors[column])[row] = value;
            }

            /// Returns a copy of the matrix.
            SparseMatrix copy() const
            {
                CONTEXT("When creating a copy:");
                SparseMatrix result(_rows, _columns, _capacity);

                for (unsigned long i(0) ; i < _rows ; ++i)
                {
                    if (_row_vectors[i])
                        result._row_vectors[i].reset(new SparseVector<DataType_>(_row_vectors[i]->copy()));
                }

                for (unsigned long i(0) ; i < _columns ; ++i)
                {
                    if (_column_vectors[i])
                        result._column_vectors[i].reset(new SparseVector<DataType_>(_column_vectors[i]->copy()));
                }

                return result;
            }


    };

    template <typename DataType_>
    void SparseMatrix<DataType_>::_synch_column_vectors()
    {
        for (typename SparseMatrix<DataType_>::NonZeroConstElementIterator i(this->begin_non_zero_elements()) ;
                i != this->end_non_zero_elements() ; ++i)
        {
            this->column(i.column())[i.row()] = *i;
        }
    }

    template <typename DataType_> struct Implementation<ConstElementIterator<storage::Sparse, container::Matrix, DataType_> >
    {
        /// Our parent matrix.
        SparseMatrix<DataType_> _matrix;

        /// Our index.
        unsigned long _index;

        /// Our row-vector's iterator.
        typename SparseVector<DataType_>::ConstElementIterator _iter;

        /// Our row index.
        unsigned long _row;

        static typename SparseVector<DataType_>::ConstElementIterator _get_iterator(const SparseMatrix<DataType_> & matrix,
                unsigned long index)
        {
            if (! matrix._row_vectors[index / matrix._columns])
                matrix._row_vectors[index / matrix._columns].reset(new SparseVector<DataType_>(matrix._columns,
                            matrix._capacity));

            return matrix._row_vectors[index / matrix._columns]->begin_elements();
        }

        Implementation(const SparseMatrix<DataType_> & matrix, unsigned long index) :
            _matrix(matrix),
            _index(index),
            _iter(_get_iterator(matrix, index)),
            _row(index / matrix._columns)
        {
        }

        Implementation(const Implementation<ElementIterator<storage::Sparse, container::Matrix, DataType_> > & other) :
            _matrix(other._matrix),
            _index(other._index),
            _iter(other._iter),
            _row(other._row)
        {
        }
    };

    template <typename DataType_>
    ConstElementIterator<storage::Sparse, container::Matrix, DataType_>::ConstElementIterator(const SparseMatrix<DataType_> & matrix, unsigned long index) :
        PrivateImplementationPattern<ConstElementIterator<storage::Sparse, container::Matrix, DataType_>, Single>(
                new Implementation<ConstElementIterator<storage::Sparse, container::Matrix, DataType_> >(matrix, index))
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::Sparse, container::Matrix, DataType_>::ConstElementIterator(const ConstElementIterator & other) :
        PrivateImplementationPattern<ConstElementIterator<storage::Sparse, container::Matrix, DataType_>, Single>(
                new Implementation<ConstElementIterator<storage::Sparse, container::Matrix, DataType_> >(*other._imp))
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::Sparse, container::Matrix, DataType_>::ConstElementIterator(
            const ElementIterator<storage::Sparse, container::Matrix, DataType_> & other) :
        PrivateImplementationPattern<ConstElementIterator<storage::Sparse, container::Matrix, DataType_>, Single>(
                new Implementation<ConstElementIterator<storage::Sparse, container::Matrix, DataType_> >(*other._imp))
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::Sparse, container::Matrix, DataType_>::~ConstElementIterator()
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::Sparse, container::Matrix, DataType_> &
    ConstElementIterator<storage::Sparse, container::Matrix, DataType_>::operator= (
            const ConstElementIterator<storage::Sparse, container::Matrix, DataType_> & other)
    {
        if (&other == this)
            return *this;

        this->_imp->_matrix = other._imp->_matrix;
        this->_imp->_index = other._imp->_index;
        this->_imp->_iter = other._imp->_iter;
        this->_imp->_row = other._imp->_row;

        return *this;
    }

    template <typename DataType_>
    ConstElementIterator<storage::Sparse, container::Matrix, DataType_> &
    ConstElementIterator<storage::Sparse, container::Matrix, DataType_>::operator++ ()
    {
        CONTEXT("When incrementing ConstElementIterator<Sparse, Matrix> by one:");

        ++this->_imp->_index;
        unsigned long row(this->_imp->_index / this->_imp->_matrix._columns);
        if (row != this->_imp->_row)
        {
            if (! this->_imp->_matrix._row_vectors[row])
                this->_imp->_matrix._row_vectors[row].reset(new SparseVector<DataType_>(this->_imp->_matrix._columns,
                            this->_imp->_matrix._capacity));

            this->_imp->_iter = this->_imp->_matrix._row_vectors[row]->begin_elements();
            this->_imp->_row = row;
        }
        else
        {
            ++this->_imp->_iter;
        }

        return *this;
    }

    template <typename DataType_>
    ConstElementIterator<storage::Sparse, container::Matrix, DataType_> &
    ConstElementIterator<storage::Sparse, container::Matrix, DataType_>::operator+= (const unsigned long step)
    {
        CONTEXT("When incrementing ConstElementIterator<Sparse, Matrix> by '" + stringify(step) + "':");

        this->_imp->_index += step;
        if (this->_imp->_matrix._rows * this->_imp->_matrix._columns >= this->_imp->_index)
            this->_imp->_index = this->_imp->_matrix._rows + this->_imp->_matrix._columns;

        unsigned long row(this->_imp->_index / this->_imp->_matrix._columns);
        if (row != this->_imp->_row)
        {
            if (! this->_imp->_matrix._row_vectors[row])
                this->_imp->_matrix._row_vectors[row].reset(new SparseVector<DataType_>(this->_imp->_matrix._columns,
                            this->_imp->_matrix._capacity));

            this->_imp->_iter = this->_imp->_matrix._row_vectors[row]->begin_elements();
            this->_imp->_row = row;
        }
        else
        {
            ++this->_imp->_iter;
        }

        return *this;
    }


    template <typename DataType_>
    const DataType_ &
    ConstElementIterator<storage::Sparse, container::Matrix, DataType_>::operator* () const
    {
        CONTEXT("When accessing unassignable element at index '" + stringify(this->_imp->_index) + "':");

        return *(this->_imp->_iter);
    }

    template <typename DataType_>
    bool
    ConstElementIterator<storage::Sparse, container::Matrix, DataType_>::operator< (
            const ConstElementIterator<storage::Sparse, container::Matrix, DataType_> & other) const
    {
        return this->_imp->_index < other._imp->_index;
    }

    /// \todo Check if both iterators belong to the same matrix
    template <typename DataType_>
    bool
    ConstElementIterator<storage::Sparse, container::Matrix, DataType_>::operator== (
            const ConstElementIterator<storage::Sparse, container::Matrix, DataType_> & other) const
    {
        return (/*(this->_imp->_matrix == other._imp->_matrix) && */(this->_imp->_index == other._imp->_index));
    }

    template <typename DataType_>
    bool
    ConstElementIterator<storage::Sparse, container::Matrix, DataType_>::operator!= (
            const ConstElementIterator<storage::Sparse, container::Matrix, DataType_> & other) const
    {
        return (/*(!(this->_imp->_matrix == other._imp->_matrix)) || */(this->_imp->_index != other._imp->_index));
    }

    template <typename DataType_>
    unsigned long
    ConstElementIterator<storage::Sparse, container::Matrix, DataType_>::index() const
    {
        return this->_imp->_index;
    }

    template <typename DataType_>
    unsigned long
    ConstElementIterator<storage::Sparse, container::Matrix, DataType_>::column() const
    {
        return this->_imp->_iter.index();
    }

    template <typename DataType_>
    unsigned long
    ConstElementIterator<storage::Sparse, container::Matrix, DataType_>::row() const
    {
        return this->_imp->_row;
    }

    template <typename DataType_> struct Implementation<ElementIterator<storage::Sparse, container::Matrix, DataType_> >
    {
        /// Our parent matrix.
        SparseMatrix<DataType_> & _matrix;

        /// Our index.
        unsigned long _index;

        /// Our row-vector's iterator.
        typename SparseVector<DataType_>::ElementIterator _iter;

        /// Our row index.
        unsigned long _row;

        typename SparseVector<DataType_>::ElementIterator _get_iterator(SparseMatrix<DataType_> & matrix,
                unsigned long index)
        {
            if (! matrix._row_vectors[index / matrix._columns])
                matrix._row_vectors[index / matrix._columns].reset(new SparseVector<DataType_>(matrix._columns,
                            matrix._capacity));

            return matrix._row_vectors[index / matrix._columns]->begin_elements();
        }

        Implementation(SparseMatrix<DataType_> & matrix, unsigned long index) :
            _matrix(matrix),
            _index(index),
            _iter(_get_iterator(matrix, index)),
            _row(index / matrix._columns)
        {
        }
    };

    template <typename DataType_>
    ElementIterator<storage::Sparse, container::Matrix, DataType_>::ElementIterator(SparseMatrix<DataType_> & matrix, unsigned long index) :
        PrivateImplementationPattern<ElementIterator<storage::Sparse, container::Matrix, DataType_>, Single>(
                new Implementation<ElementIterator<storage::Sparse, container::Matrix, DataType_> >(matrix, index))
    {
    }

    template <typename DataType_>
    ElementIterator<storage::Sparse, container::Matrix, DataType_>::ElementIterator(const ElementIterator & other) :
        PrivateImplementationPattern<ElementIterator<storage::Sparse, container::Matrix, DataType_>, Single>(
                new Implementation<ElementIterator<storage::Sparse, container::Matrix, DataType_> >(*other._imp))
    {
    }

    template <typename DataType_>
    ElementIterator<storage::Sparse, container::Matrix, DataType_>::~ElementIterator()
    {
    }

    template <typename DataType_>
    ElementIterator<storage::Sparse, container::Matrix, DataType_> &
    ElementIterator<storage::Sparse, container::Matrix, DataType_>::operator= (
            const ElementIterator<storage::Sparse, container::Matrix, DataType_> & other)
    {
        if (&other == this)
            return *this;

        this->_imp->_matrix = other._imp->_matrix;
        this->_imp->_index = other._imp->_index;
        this->_imp->_iter = other._imp->_iter;
        this->_imp->_row = other._imp->_row;

        return *this;
    }

    template <typename DataType_>
    ElementIterator<storage::Sparse, container::Matrix, DataType_> &
    ElementIterator<storage::Sparse, container::Matrix, DataType_>::operator++ ()
    {
        CONTEXT("When incrementing ElementIterator<Sparse, Matrix> by one:");

        ++this->_imp->_index;
        unsigned long row(this->_imp->_index / this->_imp->_matrix._columns);
        if (row != this->_imp->_row)
        {
            if (! this->_imp->_matrix._row_vectors[row])
                this->_imp->_matrix._row_vectors[row].reset(new SparseVector<DataType_>(this->_imp->_matrix._columns,
                            this->_imp->_matrix._capacity));

            this->_imp->_iter = this->_imp->_matrix._row_vectors[row]->begin_elements();
            this->_imp->_row = row;
        }
        else
        {
            ++this->_imp->_iter;
        }

        return *this;
    }

    template <typename DataType_>
    ElementIterator<storage::Sparse, container::Matrix, DataType_> &
    ElementIterator<storage::Sparse, container::Matrix, DataType_>::operator+= (const unsigned long step)
    {
        CONTEXT("When incrementing ElementIterator<Sparse, Matrix> by '" + stringify(step) + "':");

        this->_imp->_index += step;
        if (this->_imp->_matrix._rows * this->_imp->_matrix._columns >= this->_imp->_index)
            this->_imp->_index = this->_imp->_matrix._rows + this->_imp->_matrix._columns;

        unsigned long row(this->_imp->_index / this->_imp->_matrix._columns);
        if (row != this->_imp->_row)
        {
            if (! this->_imp->_matrix._row_vectors[row])
                this->_imp->_matrix._row_vectors[row].reset(new SparseVector<DataType_>(this->_imp->_matrix._columns,
                            this->_imp->_matrix._capacity));

            this->_imp->_iter = this->_imp->_matrix._row_vectors[row]->begin_elements();
            this->_imp->_row = row;
        }
        else
        {
            ++this->_imp->_iter;
        }

        return *this;
    }


    template <typename DataType_>
    DataType_ &
    ElementIterator<storage::Sparse, container::Matrix, DataType_>::operator* () const
    {
        CONTEXT("When accessing assignable element at index '" + stringify(this->_imp->_index) + "':");

        return *(this->_imp->_iter);
    }

    template <typename DataType_>
    bool
    ElementIterator<storage::Sparse, container::Matrix, DataType_>::operator< (
            const ElementIterator<storage::Sparse, container::Matrix, DataType_> & other) const
    {
        return this->_imp->_index < other._imp->_index;
    }

    /// \todo Check if both iterators belong to the same matrix
    template <typename DataType_>
    bool
    ElementIterator<storage::Sparse, container::Matrix, DataType_>::operator== (
            const ElementIterator<storage::Sparse, container::Matrix, DataType_> & other) const
    {
        return (/*(this->_imp->_matrix == other._imp->_matrix) && */(this->_imp->_index == other._imp->_index));
    }

    template <typename DataType_>
    bool
    ElementIterator<storage::Sparse, container::Matrix, DataType_>::operator!= (
            const ElementIterator<storage::Sparse, container::Matrix, DataType_> & other) const
    {
        return (/*(!(this->_imp->_matrix == other._imp->_matrix)) || */(this->_imp->_index != other._imp->_index));
    }

    template <typename DataType_>
    unsigned long
    ElementIterator<storage::Sparse, container::Matrix, DataType_>::index() const
    {
        return this->_imp->_index;
    }

    template <typename DataType_>
    unsigned long
    ElementIterator<storage::Sparse, container::Matrix, DataType_>::column() const
    {
        return this->_imp->_iter.index();
    }

    template <typename DataType_>
    unsigned long
    ElementIterator<storage::Sparse, container::Matrix, DataType_>::row() const
    {
        return this->_imp->_row;
    }

    template <typename DataType_> struct Implementation<ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_> >
    {
        /// Our parent matrix.
        SparseMatrix<DataType_> _matrix;

        /// Our index.
        unsigned long _index;

        /// Our row-vector's iterator for the current position.
        typename SparseVector<DataType_>::NonZeroConstElementIterator _iter;

        /// Our row-vector's iterator for the end.
        typename SparseVector<DataType_>::NonZeroConstElementIterator _end;

        /// Our column index.
        unsigned long _column;

        /// Our row index.
        unsigned long _row;

        /// Find the pair of iterators for the next existing row.
        void _find_next_row()
        {
            for ( ; _row < _matrix._rows + 1; ++_row)
            {
                if (_row == _matrix._rows)
                {
                    _iter = _matrix._row_vectors[_row]->begin_non_zero_elements();
                    _end = _matrix._row_vectors[_row]->end_non_zero_elements();
                    _index = _row * _matrix._columns;
                    return;
                }

                if (! _matrix._row_vectors[_row])
                    continue;

                _iter = _matrix._row_vectors[_row]->begin_non_zero_elements();
                _end = _matrix._row_vectors[_row]->end_non_zero_elements();


                if (_iter == _end)
                    continue;

                break;
            }

            _column = _iter.index();
            _index = _row * _matrix._columns + _column;
        }

        Implementation(const SparseMatrix<DataType_> & matrix, unsigned long index) :
                _matrix(matrix),
                _index(index),
                _iter(matrix._row_vectors[matrix._rows]->begin_non_zero_elements()), // Dummy
                _end(matrix._row_vectors[matrix._rows]->end_non_zero_elements()), // Dummy
                _column(0),
                _row(0)
        {
            if (index != 0)
            {
                _index = matrix._rows * matrix._columns;
                _row = matrix._rows;
                _find_next_row();
            }
            else
                _find_next_row();
        }

        Implementation(const Implementation<ElementIterator<storage::SparseNonZero, container::Matrix, DataType_> > & other) :
            _matrix(other._matrix),
            _index(other._index),
            _iter(other._iter),
            _end(other._end),
            _column(other._column),
            _row(other._row)
        {
        }
    };

    template <typename DataType_>
    ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_>::ConstElementIterator(const SparseMatrix<DataType_> & matrix, unsigned long index) :
        PrivateImplementationPattern<ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_>, Single>(
                new Implementation<ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_> >(matrix, index))
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_>::ConstElementIterator(const ConstElementIterator & other) :
        PrivateImplementationPattern<ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_>, Single>(
                new Implementation<ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_> >(*other._imp))
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_>::ConstElementIterator(
            const ElementIterator<storage::SparseNonZero, container::Matrix, DataType_> & other) :
        PrivateImplementationPattern<ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_>, Single>(
                new Implementation<ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_> >(*other._imp))
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_>::~ConstElementIterator()
    {
    }


    template <typename DataType_>
    ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_> &
    ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_>::operator= (
            const ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_> & other)
    {
        if (&other == this)
            return *this;

        this->_imp->_matrix = other._imp->_matrix;
        this->_imp->_index = other._imp->_index;
        this->_imp->_iter = other._imp->_iter;
        this->_imp->_end = other._imp->_end;
        this->_imp->_column = other._imp->_column;
        this->_imp->_row = other._imp->_row;

        return *this;
    }

    template <typename DataType_>
    ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_> &
    ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_>::operator++ ()
    {
        CONTEXT("When incrementing ConstElementIterator<SparseNonZero, Matrix> by one:");

        ++this->_imp->_iter;
        this->_imp->_column= this->_imp->_iter.index();

        if (this->_imp->_iter == this->_imp->_end)
        {
            ++this->_imp->_row;
            this->_imp->_column = 0;
            this->_imp->_find_next_row();
        }

        this->_imp->_index = this->_imp->_row * this->_imp->_matrix._columns + this->_imp->_column;

        return *this;
    }

    template <typename DataType_>
    ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_> &
    ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_>::operator+= (const unsigned long /*step*/)
    {
        throw InternalError("Not implemented yet!");
    }

    template <typename DataType_>
    const DataType_ &
    ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_>::operator* () const
    {
        CONTEXT("When accessing unassignable element at index '" + stringify(this->_imp->_index) + "':");

        return *(this->_imp->_iter);
    }

    template <typename DataType_>
    bool
    ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_>::operator< (
            const ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_> & other) const
    {
        return this->_imp->_index < other._imp->_index;
    }

    /// \todo Check if both iterators belong to the same matrix
    template <typename DataType_>
    bool
    ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_>::operator== (
            const ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_> & other) const
    {
        return (/*(this->_imp->_matrix == other._imp->_matrix) && */(this->_imp->_index == other._imp->_index));
    }

    template <typename DataType_>
    bool
    ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_>::operator!= (
            const ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_> & other) const
    {
        return (/*(!(this->_imp->_matrix == other._imp->_matrix)) || */(this->_imp->_index != other._imp->_index));
    }

    template <typename DataType_>
    unsigned long
    ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_>::index() const
    {
        return this->_imp->_index;
    }

    template <typename DataType_>
    unsigned long
    ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_>::column() const
    {
        return this->_imp->_column;
    }

    template <typename DataType_>
    unsigned long
    ConstElementIterator<storage::SparseNonZero, container::Matrix, DataType_>::row() const
    {
        return this->_imp->_row;
    }

    template <typename DataType_> struct Implementation<ElementIterator<storage::SparseNonZero, container::Matrix, DataType_> >
    {
        /// Our parent matrix.
        SparseMatrix<DataType_> & _matrix;

        /// Our index.
        unsigned long _index;

        /// Our row-vector's iterator for the current position.
        typename SparseVector<DataType_>::NonZeroElementIterator _iter;

        /// Our row-vector's iterator for the end.
        typename SparseVector<DataType_>::NonZeroElementIterator _end;

        /// Our column index.
        unsigned long _column;

        /// Our row index.
        unsigned long _row;

        /// Find the pair of iterators for the next existing row.
        void _find_next_row()
        {
            for ( ; _row < _matrix._rows + 1; ++_row)
            {
                if (_row == _matrix._rows)
                {
                    _iter = _matrix._row_vectors[_row]->begin_non_zero_elements();
                    _end = _matrix._row_vectors[_row]->end_non_zero_elements();
                    _index = _row * _matrix._columns;
                    return;
                }

                if (! _matrix._row_vectors[_row])
                    continue;

                _iter = _matrix._row_vectors[_row]->begin_non_zero_elements();
                _end = _matrix._row_vectors[_row]->end_non_zero_elements();


                if (_iter == _end)
                    continue;

                break;
            }

            _column = _iter.index();
            _index = _row * _matrix._columns + _column;
        }

        Implementation(SparseMatrix<DataType_> & matrix, unsigned long index) :
                _matrix(matrix),
                _index(index),
                _iter(matrix._row_vectors[matrix._rows]->begin_non_zero_elements()), // Dummy
                _end(matrix._row_vectors[matrix._rows]->end_non_zero_elements()), // Dummy
                _column(0),
                _row(0)
        {
            if (index != 0)
            {
                _index = matrix._rows * matrix._columns;
                _row = matrix._rows;
                _find_next_row();
            }
            else
                _find_next_row();
        }
    };

    template <typename DataType_>
    ElementIterator<storage::SparseNonZero, container::Matrix, DataType_>::ElementIterator(SparseMatrix<DataType_> & matrix, unsigned long index) :
        PrivateImplementationPattern<ElementIterator<storage::SparseNonZero, container::Matrix, DataType_>, Single>(
                new Implementation<ElementIterator<storage::SparseNonZero, container::Matrix, DataType_> >(matrix, index))
    {
    }

    template <typename DataType_>
    ElementIterator<storage::SparseNonZero, container::Matrix, DataType_>::ElementIterator(const ElementIterator & other) :
        PrivateImplementationPattern<ElementIterator<storage::SparseNonZero, container::Matrix, DataType_>, Single>(
                new Implementation<ElementIterator<storage::SparseNonZero, container::Matrix, DataType_> >(*other._imp))
    {
    }

    template <typename DataType_>
    ElementIterator<storage::SparseNonZero, container::Matrix, DataType_>::~ElementIterator()
    {
    }


    template <typename DataType_>
    ElementIterator<storage::SparseNonZero, container::Matrix, DataType_> &
    ElementIterator<storage::SparseNonZero, container::Matrix, DataType_>::operator= (
            const ElementIterator<storage::SparseNonZero, container::Matrix, DataType_> & other)
    {
        if (&other == this)
            return *this;

        this->_imp->_matrix = other._imp->_matrix;
        this->_imp->_index = other._imp->_index;
        this->_imp->_iter = other._imp->_iter;
        this->_imp->_end = other._imp->_end;
        this->_imp->_column = other._imp->_column;
        this->_imp->_row = other._imp->_row;

        return *this;
    }

    template <typename DataType_>
    ElementIterator<storage::SparseNonZero, container::Matrix, DataType_> &
    ElementIterator<storage::SparseNonZero, container::Matrix, DataType_>::operator++ ()
    {
        CONTEXT("When incrementing ElementIterator<SparseNonZero, Matrix> by one:");

        ++this->_imp->_iter;
        this->_imp->_column = this->_imp->_iter.index();

        if (this->_imp->_iter == this->_imp->_end)
        {
            ++this->_imp->_row;
            this->_imp->_column = 0;
            this->_imp->_find_next_row();
        }

        this->_imp->_index = this->_imp->_row * this->_imp->_matrix._columns + this->_imp->_column;

        return *this;
    }

    template <typename DataType_>
    ElementIterator<storage::SparseNonZero, container::Matrix, DataType_> &
    ElementIterator<storage::SparseNonZero, container::Matrix, DataType_>::operator+= (const unsigned long /*step*/)
    {
        throw InternalError("Not implemented yet!");
    }

    template <typename DataType_>
    DataType_ &
    ElementIterator<storage::SparseNonZero, container::Matrix, DataType_>::operator* () const
    {
        CONTEXT("When accessing assignable element at index '" + stringify(this->_imp->_index) + "':");

        return *(this->_imp->_iter);
    }

    template <typename DataType_>
    bool
    ElementIterator<storage::SparseNonZero, container::Matrix, DataType_>::operator< (
            const ElementIterator<storage::SparseNonZero, container::Matrix, DataType_> & other) const
    {
        return this->_imp->_index < other._imp->_index;
    }

    /// \todo Check if both iterators belong to the same matrix
    template <typename DataType_>
    bool
    ElementIterator<storage::SparseNonZero, container::Matrix, DataType_>::operator== (
            const ElementIterator<storage::SparseNonZero, container::Matrix, DataType_> & other) const
    {
        return (/*(this->_imp->_matrix == other._imp->_matrix) && */(this->_imp->_index == other._imp->_index));
    }

    template <typename DataType_>
    bool
    ElementIterator<storage::SparseNonZero, container::Matrix, DataType_>::operator!= (
            const ElementIterator<storage::SparseNonZero, container::Matrix, DataType_> & other) const
    {
        return (/*(!(this->_imp->_matrix == other._imp->_matrix)) || */(this->_imp->_index != other._imp->_index));
    }

    template <typename DataType_>
    unsigned long
    ElementIterator<storage::SparseNonZero, container::Matrix, DataType_>::index() const
    {
        return this->_imp->_index;
    }

    template <typename DataType_>
    unsigned long
    ElementIterator<storage::SparseNonZero, container::Matrix, DataType_>::column() const
    {
        return this->_imp->_column;
    }

    template <typename DataType_>
    unsigned long
    ElementIterator<storage::SparseNonZero, container::Matrix, DataType_>::row() const
    {
        return this->_imp->_row;
    }




    /**
     * \brief SparseMatrix::NonZeroRowIterator is a smart iterator implementation that iterates over non-zero
     * \brief rows of sparse matrices.
     *
     * \ingroup grpmatrix
     **/
    /*template <> template <typename DataType_> class SparseMatrix<DataType_>::NonZeroRowIterator :
        public VectorIteratorBase<DataType_, SparseVector<DataType_> >
    {
        private:
            /// Our parent matrix.
            const SparseMatrix<DataType_> & _matrix;

            /// Our index.
            unsigned long _index;

        public:
            /// \name Constructors and destructor
            /// \{

            NonZeroRowIterator(const SparseMatrix<DataType_> & matrix, const unsigned long & index) :
                _matrix(matrix),
                _index(index)
            {
            }

            /// Copy-constructor.
            NonZeroRowIterator(NonZeroRowIterator const & other) :
                _matrix(other._matrix),
                _index(other._index)
            {
            }

            /// \}

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            virtual NonZeroRowIterator & operator++ ()
            {
                while (_index < _matrix._rows)
                {
                    _index++;
                    if (_matrix._row_vectors[_index])
                        break;
                }
                return *this;
            }

            /// Dereference operator that returns assignable reference.
            virtual SparseVector<DataType_> & operator* ()
            {
                CONTEXT("When accessing assignable non-zero-row at index '" + stringify(_index) + "':");

                return *_matrix._row_vectors[_index];
            }

            /// Dereference operator that returns umassignable reference.
            virtual const SparseVector<DataType_> & operator* () const
            {
                CONTEXT("When accessing unassignable non-zero-row at index '" + stringify(_index) + "':");

                return *_matrix._row_vectors[_index];
            }

            /// Comparison operator for less-than.
            virtual bool operator< (const VectorIteratorBase<DataType_, SparseVector<DataType_> > & other) const
            {
                return _index < other.index();
            }

            /// Comparison operator for equality.
            virtual bool operator== (const VectorIteratorBase<DataType_, SparseVector<DataType_> > & other) const
            {
                return ((&_matrix == other.parent()) && (_index == other.index()));
            }

            /// Comparison operator for inequality.
            virtual bool operator!= (const VectorIteratorBase<DataType_, SparseVector<DataType_> > & other) const
            {
                return ((&_matrix != other.parent()) || (_index != other.index()));
            }

            /// \}


            /// \name IteratorTraits interface
            /// \{

            /// Returns true if the referenced vector already exists.
            virtual bool exists() const
            {
                return _matrix._row_vectors[_index];
            }

            /// Returns our index.
            virtual const unsigned long index() const
            {
                return _index;
            }

            /// Returns a pointer to our parent container.
            virtual const Matrix<DataType_> * parent() const
            {
                return &_matrix;
            }

            /// \}
    };*/

    template <typename DataType_>
    std::ostream &
    operator<< (std::ostream & lhs, const SparseMatrix<DataType_> & b)
    {
        lhs << "[" << std::endl;
        for (unsigned long row(0) ; row < b.rows() ; ++row)
        {
            lhs << b[row];
        }
        lhs << "]" << std::endl;

        return lhs;
    }


    template <typename DataType_>
    bool
    operator== (const SparseMatrix<DataType_> & a, const SparseMatrix<DataType_> & b)
    {
        CONTEXT("When comparing two sparse matrices:");

        if (a.columns() != b.columns())
        {
            throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
        }

        if (a.rows() != b.rows())
        {
            throw MatrixRowsDoNotMatch(b.rows(), a.rows());
        }

        for (unsigned long row(0) ; row < a.rows() ; ++row)
        {
            if (! (a[row] == b[row]))
                return false;
        }

        return true;
    }
}

#endif
