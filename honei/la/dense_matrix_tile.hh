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

#ifndef LIBLA_GUARD_DENSE_MATRIX_TILE_HH
#define LIBLA_GUARD_DENSE_MATRIX_TILE_HH 1

#include <honei/la/element_iterator.hh>
#include <honei/la/matrix.hh>
#include <honei/la/dense_vector_range.hh>
#include <honei/la/dense_vector_slice.hh>
#include <honei/la/matrix_error.hh>
#include <honei/util/shared_array-impl.hh>
#include <honei/util/stringify.hh>

#include <iterator>

namespace honei
{
    /**
     * \brief DenseMatrixTile is a matrix with O(column * row) non-zero elements.
     * \brief Its purpose is to provide access to a part of a DenseMatrix for more
     * \brief efficient handling.
     *
     * \ingroup grpmatrix
     */
    template <typename DataType_> class DenseMatrixTile :
        public MutableMatrix<DataType_> ,
        public RowAccessMatrix<DataType_>
    {
        private:
            /// Pointer to our elements.
            SharedArray<DataType_> _elements;

            /// Our columns.
            unsigned long _columns;

            /// Our rows.
            unsigned long _rows;

            /// The number of columns of our source matrix.
            unsigned long _source_columns;

            /// Our row-offset.
            unsigned long _row_offset;

            /// Our column-offset.
            unsigned long _column_offset;

            /// Our row-vectors.
            SharedArray<std::tr1::shared_ptr<DenseVectorRange<DataType_> > > _row_vectors;

            /// Our column-vectors.
            SharedArray<std::tr1::shared_ptr<DenseVectorSlice<DataType_> > > _column_vectors;

            /// Our implementation of ElementIteratorBase.
            class DenseElementIterator;

            typedef typename Matrix<DataType_>::MatrixElementIterator MatrixElementIterator;

        public:
            friend class DenseElementIterator;

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
             * \param source The matrix our tile provides access to.
             * \param rows Number of rows of the new matrix.
             * \param columns Number of columns of the new matrix.
             * \param row_offset The row-offset inside the source matrix.
             * \param column_offset The column-offset inside the source matrix.
             */
            DenseMatrixTile(const DenseMatrix<DataType_> & source, const unsigned long rows, const unsigned long columns,
                                const unsigned long row_offset, const unsigned long column_offset) :
                _elements(source._elements),
                _columns(columns),
                _column_vectors(columns),
                _rows(rows),
                _row_vectors(rows),
                _row_offset(row_offset),
                _column_offset(column_offset),
                _source_columns(source._columns)
            {
                CONTEXT("When creating DenseMatrixTile from DenseMatrix: ");
                ASSERT(rows > 0, "number of rows is zero!");
                ASSERT(columns > 0, "number of columns is zero!");
            }

            DenseMatrixTile(const DenseMatrixTile<DataType_> & source, const unsigned long rows, const unsigned long columns,
                                const unsigned long row_offset, const unsigned long column_offset) :
                _elements(source._elements),
                _columns(columns),
                _column_vectors(columns),
                _rows(rows),
                _row_vectors(rows),
                _row_offset(source._row_offset + row_offset),
                _column_offset(source._column_offset + column_offset),
                _source_columns(source._source_columns)
            {
                CONTEXT("When creating DenseMatrixTile from DenseMatrixTile: ");
                ASSERT(rows > 0, "number of rows is zero!");
                ASSERT(columns > 0, "number of columns is zero!");
            }

            /// Returns iterator pointing to the first element of the matrix.
            virtual ConstElementIterator begin_elements() const
            {
                return ConstElementIterator(new DenseElementIterator(*this, 0));
            }

            /// Returns iterator pointing behind the last element of the matrix.
            virtual ConstElementIterator end_elements() const
            {
                return ConstElementIterator(
                        new DenseElementIterator(*this, _rows * _columns));
            }

            /// Returns iterator pointing to the first element of the matrix.
            virtual ElementIterator begin_elements()
            {
                return ElementIterator(new DenseElementIterator(*this, 0));
            }

            /// Returns iterator pointing behind the last element of the matrix.
            virtual ElementIterator end_elements()
            {
                return ElementIterator(
                        new DenseElementIterator(*this, _rows * _columns));
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
            virtual const DenseVectorRange<DataType_> & operator[] (unsigned long row) const
            {
                if (! _row_vectors[row])
                    _row_vectors[row].reset(
                            new DenseVectorRange<DataType_>(_elements, _columns, (_row_offset + row) * _source_columns + _column_offset) );

                return *_row_vectors[row];
            }

            /// Retrieves row vector by index, zero-based, assignable.
            virtual DenseVectorRange<DataType_> & operator[] (unsigned long row)
            {
                if (! _row_vectors[row])
                    _row_vectors[row].reset(
                            new DenseVectorRange<DataType_>(_elements, _columns, (_row_offset + row) * _source_columns + _column_offset) );

                return *_row_vectors[row];
            }

            /// Retrieves element at (row, column), unassignable.
            virtual const DataType_ & operator() (unsigned long row, unsigned long column) const
            {
                return _elements[column + _column_offset + (row + _row_offset) * _source_columns];
            }

            /// Retrieves element at (row, column), assignable.
            virtual DataType_ & operator() (unsigned long row, unsigned long column)
            {
                return _elements[column + _column_offset + (row + _row_offset) * _source_columns];
            }

            /// Retrieves column vector by index, zero-based, unassignable.
            virtual const DenseVectorSlice<DataType_> & column(unsigned long column) const
            {
                if (! _column_vectors[column])
                    _column_vectors[column].reset(
                            new DenseVectorSlice<DataType_>(_elements, _rows,
                                _row_offset * _source_columns + _column_offset + column, _source_columns));

                return *_column_vectors[column];
            }

            /// Retrieves column vector by index, zero-based, assignable.
            virtual DenseVectorSlice<DataType_> & column(unsigned long column)
            {
                if (! _column_vectors[column])
                    _column_vectors[column].reset(
                            new DenseVectorSlice<DataType_>(_elements, _rows,
                                _row_offset * _source_columns + _column_offset + column, _source_columns));

                return *_column_vectors[column];
            }

            /// Returns a pointer to our data array.
            inline DataType_ * elements() const
            {
               return _elements.get() + _row_offset * _source_columns + _column_offset;
            }

            /// Returns a copy of the matrix.
            virtual DenseMatrix<DataType_> copy() const
            {
                DenseMatrix<DataType_> result(_rows, _columns);
                DataType_ * source(_elements.get());
                DataType_ * target(result.elements());

                for (unsigned long i(0) ; i < _rows ; i++)
                {
                    for (unsigned long j(0) ; j < _columns ; ++j)
                    {
                        target[j + i * _rows] = source[(_row_offset + i) * _source_columns + _column_offset + j];
                    }
                }
                ///\todo: Use TypeTraits.

                return result;
            }
    };

    /**
     * \brief DenseMatrixTile::DenseElementIterator is a simple iterator implementation for dense matrix tiles.
     *
     * \ingroup grpmatrix
     */
    template <> template <typename DataType_> class DenseMatrixTile<DataType_>::DenseElementIterator :
        public MatrixElementIterator
    {
        private:
            /// Our parent matrix.
            const DenseMatrixTile<DataType_> & _tile;

            /// Our index.
            unsigned long _index;

            /// Our position inside the matrix's elements.
            unsigned long _pos;

            /// Calculates our position from a given index.
            inline unsigned long _calc_pos(const unsigned long index)
            {
                return (_tile._row_offset + index / _tile._columns) * _tile._source_columns + _tile._column_offset + index % _tile._columns;
            }

        public:
            /// \name Constructors
            /// \{

            /**
             * Constructor.
             *
             * \param matrix The parent matrix tile that is referenced by the iterator.
             * \param index The index into the matrix.
             */
            DenseElementIterator(const DenseMatrixTile<DataType_> & tile, unsigned long index) :
                _tile(tile),
                _index(index),
                _pos(_calc_pos(index))
            {
            }

            /// Copy-constructor.
            DenseElementIterator(DenseElementIterator const & other) :
                _tile(other._tile),
                _index(other._index),
                _pos(other._pos)
            {
            }

            /// \}

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            virtual DenseElementIterator & operator++ ()
            {
                CONTEXT("When incrementing iterator by one:");

                ++_index;
                if (_index % _tile._columns)
                {
                    ++_pos;
                }
                else
                {
                    _pos += _tile._source_columns - _tile._columns + 1;
                }

                return *this;
            }

            /// In-place-add operator.
            virtual DenseElementIterator & operator+= (const unsigned long step)
            {
                CONTEXT("When incrementing iterator by '" + stringify(step) + "':");

                _index += step;
                _pos = _calc_pos(_index);

                return *this;
            }

            /// Dereference operator that returns assignable reference.
            virtual DataType_ & operator* ()
            {
                CONTEXT("When accessing assignable element at index '" + stringify(_index) + "':");

                return _tile._elements[_pos];
            }

            /// Dereference operator that returns umassignable reference.
            virtual const DataType_ & operator* () const
            {
                CONTEXT("When accessing unassignable element at index '" + stringify(_index) + "':");

                return _tile._elements[_pos];
            }

            /// Comparison operator for less-than.
            virtual bool operator< (const IteratorBase<DataType_, Matrix<DataType_> > & other) const
            {
                return _index < other.index();
            }

            /// Comparison operator for equality.
            virtual bool operator== (const IteratorBase<DataType_, Matrix<DataType_> > & other) const
            {
                return ((&_tile == other.parent()) && (_index == other.index()));
            }

            /// Comparison operator for inequality.
            virtual bool operator!= (const IteratorBase<DataType_, Matrix<DataType_> > & other) const
            {
                return ((&_tile != other.parent()) || (_index != other.index()));
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
                return _index % _tile._columns;
            }

            /// Returns our row index.
            virtual unsigned long row() const
            {
                return _index / _tile._columns;
            }

            /// Returns a pointer to our parent container.
            virtual const Matrix<DataType_> * parent() const
            {
                return &_tile;
            }

            /// \}
    };
}

#endif
