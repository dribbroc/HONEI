/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_DENSE_MATRIX_HH
#define LIBLA_GUARD_DENSE_MATRIX_HH 1

#include <libla/element_iterator.hh>
#include <libla/matrix.hh>
#include <libla/dense_matrix_tile.hh>
#include <libla/dense_vector.hh>
#include <libla/sparse_matrix.hh>
#include <libla/matrix_error.hh>
#include <libutil/shared_array-impl.hh>
#include <libutil/stringify.hh>
#include <libutil/type_traits.hh>

#include <algorithm>
#include <cstring>
#include <iterator>

/**
 * \file
 *
 * Implementation of DenseMatrix and related classes.
 *
 * \ingroup grpmatrix
 */
namespace honei
{
    template <typename DataType_> class DenseMatrixTile;

    /**
     * \brief DenseMatrix is a matrix with O(column * row) non-zero elements which keeps its data
     * \brief continuous.
     *
     * \ingroup grpmatrix
     */
    template <typename DataType_> class DenseMatrix :
        public RowAccessMatrix<DataType_>,
        public MutableMatrix<DataType_>
    {
        private:
            /// Pointer to our elements.
            SharedArray<DataType_> _elements;

            /// Our columns.
            unsigned long _columns;

            /// Our rows.
            unsigned long _rows;

            /// Our row-vectors.
            SharedArray<std::tr1::shared_ptr<DenseVector<DataType_> > > _row_vectors;

            /// Our column-vectors.
            SharedArray<std::tr1::shared_ptr<DenseVector<DataType_> > > _column_vectors;

            /// Our implementation of ElementIteratorBase.
            template <typename ElementType_> class DenseElementIterator;

            typedef typename Matrix<DataType_>::MatrixElementIterator MatrixElementIterator;

        public:
            friend class DenseElementIterator<DataType_>;
            friend class DenseMatrixTile<DataType_>;

            /// Type of the const iterator over our elements.
            typedef typename Matrix<DataType_>::ConstElementIterator ConstElementIterator;

            /// Type of the iterator over our elements.
            typedef typename MutableMatrix<DataType_>::ElementIterator ElementIterator;

            /// Type of the const iterator over our vectors.
            typedef VectorIteratorWrapper<DataType_, const DataType_> ConstVectorIterator;

            /// Type of the iterator over our vectors.
            typedef VectorIteratorWrapper<DataType_, DataType_> VectorIterator;

            /// \name Constructors
            /// \{

            /**
             * Constructor.
             *
             * \param columns Number of columns of the new matrix.
             * \param rows Number of rows of the new matrix.
             */
            DenseMatrix(unsigned long rows, unsigned long columns) :
                _elements(rows * columns),
                _columns(columns),
                _column_vectors(columns),
                _rows(rows),
                _row_vectors(rows)
            {
                CONTEXT("When creating DenseMatrix:");
            }

            /**
             * Constructor.
             *
             * \param columns Number of columns of the new matrix.
             * \param rows Number of rows of the new matrix.
             * \param value Default value of each of the new matrice's elements.
             */
            DenseMatrix(unsigned long rows, unsigned long columns, DataType_ value) :
                _elements(rows * columns),
                _columns(columns),
                _column_vectors(columns),
                _rows(rows),
                _row_vectors(rows)
            {
                CONTEXT("When creating DenseMatrix:");
                DataType_ *  target(_elements.get());
                for (unsigned long i(0) ; i < (rows * columns) ; i++)
                    target[i] = value;
            }

            /**
             * Constructor.
             *
             * \param other The SparseMatrix to densify.
             */
            DenseMatrix(const SparseMatrix<DataType_> & other) :
                _elements(other.rows() * other.columns()),
                _columns(other.columns()),
                _column_vectors(other.columns()),
                _rows(other.rows()),
                _row_vectors(other.rows())
            {
                CONTEXT("When creating DenseMatrix form SparseMatrix:");

                /// \todo Use TypeTraits::zero()
                DataType_ *  target(_elements.get());
                DataType_ value(0);
                for (unsigned long i(0) ; i < (other.rows() * other.columns()) ; i++)
                    target[i] = value;

                for (typename Matrix<DataType_>::ConstElementIterator i(other.begin_non_zero_elements()),
                        i_end(other.end_non_zero_elements()) ; i != i_end ; ++i)
                {
                    (*this)(i.row(),i.column()) = *i;
                }
            }

            /**
             * Constructor
             *
             * Create a submatrix from a given source matrix.
             * \param source The source matrix.
             * \param column_offset The source matrix column offset.
             * \param columns Number of columns of the new matrix.
             * \param row_offset The source matrix row offset.
             * \param rows Number of rows of the new matrix.
             */
            DenseMatrix(const DenseMatrix<DataType_> & source, unsigned long column_offset, unsigned long columns,
                    unsigned long row_offset, unsigned long rows) :
                    _elements(columns * rows),
                    _columns(columns),
                    _column_vectors(columns),
                    _rows(rows),
                    _row_vectors(rows)
            {
                if (column_offset + columns > source.columns())
                {
                    throw MatrixColumnsDoNotMatch(column_offset + columns, source.columns());
                }

                if (row_offset + rows > source.rows())
                {
                    throw MatrixRowsDoNotMatch(row_offset + rows, source.rows());
                }

                for (unsigned long i = 0 ; i < rows ; ++i)
                {
                    for (unsigned long j = 0; j < columns ; ++j)
                    {
                        _elements[j + columns * i] = source._elements[j + column_offset  +
                            ((i + row_offset) * source.columns())];
                    }
                }
            }
            /// \}

            /// Returns iterator pointing to the first element of the matrix.
            virtual ConstElementIterator begin_elements() const
            {
                return ConstElementIterator(new DenseElementIterator<DataType_>(*this, 0));
            }

            /// Returns iterator pointing to a given element of the matrix.
            virtual ConstElementIterator element_at(unsigned long index) const
            {
                 return ConstElementIterator(new DenseElementIterator<DataType_>(*this, index));
            }

            /// Returns iterator pointing behind the last element of the matrix.
            virtual ConstElementIterator end_elements() const
            {
                return ConstElementIterator(new DenseElementIterator<DataType_>(*this, _rows * _columns));
            }

            /// Returns iterator pointing to the first element of the matrix.
            virtual ElementIterator begin_elements()
            {
                return ElementIterator(new DenseElementIterator<DataType_>(*this, 0));
            }

            /// Returns iterator pointing to a given element of the matrix.
            virtual ElementIterator element_at(unsigned long index)
            {
                 return ElementIterator(new DenseElementIterator<DataType_>(*this, index));
            }

           /// Returns iterator pointing behind the last element of the matrix.
            virtual ElementIterator end_elements()
            {
                return ElementIterator(new DenseElementIterator<DataType_>(*this, _rows * _columns));
            }

            /// Returns the number of our columns.
            virtual unsigned long columns() const
            {
                return _columns;
            }

            /// Returns the number of our rows.
            virtual unsigned long rows() const
            {
                return _rows;
            }

            /// Retrieves row vector by index, zero-based, unassignable.
            virtual const DenseVector<DataType_> & operator[] (unsigned long row) const
            {
                if (! _row_vectors[row])
                    _row_vectors[row].reset(new DenseVector<DataType_>(_columns, _elements, row * _columns, 1));

                return *_row_vectors[row];
            }

            /// Retrieves row vector by index, zero-based, assignable.
            virtual DenseVector<DataType_> & operator[] (unsigned long row)
            {
                if (! _row_vectors[row])
                    _row_vectors[row].reset(new DenseVector<DataType_>(_columns, _elements, row * _columns, 1));

                return *_row_vectors[row];
            }

            /// Retrieves element at (row, column), unassignable.
            virtual const DataType_ & operator() (unsigned long row, unsigned long column) const
            {
                return _elements.get()[column + row * _columns]; 
            }

            /// Retrieves element at (row, column), assignable.
            virtual DataType_ & operator() (unsigned long row, unsigned long column)
            {
                return _elements.get()[column + row * _columns]; 
            }

            /// Retrieves column vector by index, zero-based, unassignable.
            virtual const DenseVector<DataType_> & column(unsigned long column) const
            {
                if (! _column_vectors[column])
                    _column_vectors[column].reset(new DenseVector<DataType_>(_rows, _elements, column, _columns));

                return *_column_vectors[column];
            }

            /// Retrieves column vector by index, zero-based, assignable.
            virtual DenseVector<DataType_> & column(unsigned long column)
            {
                if (! _column_vectors[column])
                    _column_vectors[column].reset(new DenseVector<DataType_>(_rows, _elements, column, _columns));

                return *_column_vectors[column];
            }

            /// Returns a pointer to our data array.
            inline DataType_ * elements() const
            {
                return _elements.get();
            }

            /// Returns a copy of the matrix.
            DenseMatrix * copy() const
            {
                DenseMatrix * result(new DenseMatrix(_columns, _rows));

                TypeTraits<DataType_>::copy(_elements.get(), result->_elements.get(), _columns * _rows);

                return result;
            }
    };

    /**
     * \brief DenseMatrix::DenseElementIterator is a simple iterator implementation for dense matrices.
     *
     * \ingroup grpmatrix
     */
    template <> template <typename DataType_> class DenseMatrix<DataType_>::DenseElementIterator<DataType_> :
        public MatrixElementIterator
    {
        private:
            /// Our parent matrix.
            const DenseMatrix<DataType_> & _matrix;

            /// Our index.
            unsigned long _index;

        public:
            /// \name Constructors
            /// \{

            /**
             * Constructor.
             *
             * \param matrix The parent matrix that is referenced by the iterator.
             * \param index The index into the matrix.
             */
            DenseElementIterator(const DenseMatrix<DataType_> & matrix, unsigned long index) :
                _matrix(matrix),
                _index(index)
            {
            }

            /// Copy-constructor.
            DenseElementIterator(DenseElementIterator<DataType_> const & other) :
                _matrix(other._matrix),
                _index(other._index)
            {
            }

            /// \}

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            virtual DenseElementIterator<DataType_> & operator++ ()
            {
                CONTEXT("When incrementing iterator by one:");

                ++_index;

                return *this;
            }

            /// In-place-add operator.
            virtual DenseElementIterator<DataType_> & operator+= (const unsigned long step)
            {
                CONTEXT("When incrementing iterator by '" + stringify(step) + "':");

                _index += step;

                return *this;
            }

            /// Dereference operator that returns assignable reference.
            virtual DataType_ & operator* ()
            {
                CONTEXT("When accessing assignable element at index '" + stringify(_index) + "':");

                return _matrix._elements[_index];
            }

            /// Dereference operator that returns umassignable reference.
            virtual const DataType_ & operator* () const
            {
                CONTEXT("When accessing unassignable element at index '" + stringify(_index) + "':");

                return _matrix._elements[_index];
            }

            /// Comparison operator for less-than.
            virtual bool operator< (const IteratorBase<DataType_, Matrix<DataType_> > & other) const
            {
                return _index < other.index();
            }

            /// Comparison operator for equality.
            virtual bool operator== (const IteratorBase<DataType_, Matrix<DataType_> > & other) const
            {
                return ((&_matrix == other.parent()) && (_index == other.index()));
            }

            /// Comparison operator for inequality.
            virtual bool operator!= (const IteratorBase<DataType_, Matrix<DataType_> > & other) const
            {
                return ((&_matrix != other.parent()) || (_index != other.index()));
            }

            /// \}

            /// \name IteratorTraits interface
            /// \{

            /// Returns our index.
            virtual unsigned long index() const
            {
                return _index;
            }

            /// Returns our column index.
            virtual unsigned long column() const
            {
                return _index % _matrix._columns;
            }

            /// Returns our row index.
            virtual unsigned long row() const
            {
                return _index / _matrix._columns;
            }

            /// Returns a pointer to our parent container.
            virtual const Matrix<DataType_> * parent() const
            {
                return &_matrix;
            }

            /// \}
    };
}

#endif
