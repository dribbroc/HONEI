/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Volker Jung <volker.jung@uni-dortmund.de>
 * Copyright (c) 2009 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
#include <honei/la/dense_vector_range.hh>
#include <honei/la/dense_vector_slice.hh>
#include <honei/la/dense_matrix-impl.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/matrix_error.hh>
#include <honei/util/shared_array-impl.hh>
#include <honei/util/stringify.hh>

#include <cmath>
#include <iterator>
#include <limits>
#include <ostream>

namespace honei
{
    /**
     * \brief DenseMatrixTile is a matrix with O(column * row) non-zero elements.
     * \brief Its purpose is to provide access to a part of a DenseMatrix for more
     * \brief efficient handling.
     *
     * \ingroup grpmatrix
     */
    template <typename DataType_> class DenseMatrixTile
    {
        private:
            /// Pointer to our elements.
            SharedArray<DataType_> _elements;

            /// Our rows.
            unsigned long _rows;

            /// Our columns.
            unsigned long _columns;

            /// Our row-vectors.
            SharedArray<std::tr1::shared_ptr<DenseVectorRange<DataType_> > > _row_vectors;

            /// Our column-vectors.
            SharedArray<std::tr1::shared_ptr<DenseVectorSlice<DataType_> > > _column_vectors;

            /// Our row-offset.
            unsigned long _row_offset;

            /// Our column-offset.
            unsigned long _column_offset;

            /// The number of columns of our source matrix.
            unsigned long _source_columns;

            /// Calculates our position from a given index.
            inline unsigned long _calc_pos(const unsigned long index) const
            {
                return (_row_offset + index / _columns) * _source_columns + _column_offset + index % _columns;
            }


        public:
            friend class honei::ConstElementIterator<storage::Dense, container::MatrixTile, DataType_>;
            friend class honei::ElementIterator<storage::Dense, container::MatrixTile, DataType_>;
            /// \todo remove impl friends after _calc_pos is moved
            friend struct honei::Implementation<honei::ConstElementIterator<storage::Dense, container::MatrixTile, DataType_> >;
            friend struct honei::Implementation<honei::ElementIterator<storage::Dense, container::MatrixTile, DataType_> >;

            /// Type of the const iterator over our elements.
            typedef honei::ConstElementIterator<storage::Dense, container::MatrixTile, DataType_> ConstElementIterator;

            /// Type of the iterator over our elements.
            typedef honei::ElementIterator<storage::Dense, container::MatrixTile, DataType_> ElementIterator;

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
                _elements(source._imp->elements),
                _rows(rows),
                _columns(columns),
                _row_vectors(rows),
                _column_vectors(columns),
                _row_offset(row_offset),
                _column_offset(column_offset),
                _source_columns(source.columns())
            {
                CONTEXT("When creating DenseMatrixTile from DenseMatrix:");
                ASSERT(rows > 0, "number of rows is zero!");
                ASSERT(columns > 0, "number of columns is zero!");
            }

            DenseMatrixTile(const DenseMatrixTile<DataType_> & source, const unsigned long rows, const unsigned long columns,
                                const unsigned long row_offset, const unsigned long column_offset) :
                _elements(source._elements),
                _rows(rows),
                _columns(columns),
                _row_vectors(rows),
                _column_vectors(columns),
                _row_offset(source._row_offset + row_offset),
                _column_offset(source._column_offset + column_offset),
                _source_columns(source._source_columns)
            {
                CONTEXT("When creating DenseMatrixTile from DenseMatrixTile:");
                ASSERT(rows > 0, "number of rows is zero!");
                ASSERT(columns > 0, "number of columns is zero!");
            }

            /// Destructor
            virtual ~DenseMatrixTile()
            {
            }

            /// Returns iterator pointing to the first element of the matrix.
            virtual ConstElementIterator begin_elements() const
            {
                return ConstElementIterator(*this, 0);
            }

            /// Returns iterator pointing behind the last element of the matrix.
            virtual ConstElementIterator end_elements() const
            {
                return ConstElementIterator(*this, _rows * _columns);
            }

            /// Returns iterator pointing to the first element of the matrix.
            virtual ElementIterator begin_elements()
            {
                return ElementIterator(*this, 0);
            }

            /// Returns iterator pointing behind the last element of the matrix.
            virtual ElementIterator end_elements()
            {
                return ElementIterator(*this, _rows * _columns);
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

    template <typename DataType_> struct Implementation<ConstElementIterator<storage::Dense, container::MatrixTile, DataType_> >
    {
            /// Our parent matrix.
            DenseMatrixTile<DataType_> _tile;

            /// Our index.
            unsigned long _index;

            /// Our position inside the matrix's elements.
            unsigned long _pos;


        Implementation(const DenseMatrixTile<DataType_> & tile, unsigned long index) :
            _tile(tile),
            _index(index),
            _pos(tile._calc_pos(index))
        {
        }

        Implementation(const Implementation<ElementIterator<storage::Dense, container::MatrixTile, DataType_> > & other) :
            _tile(other._tile),
            _index(other._index),
            _pos(other._pos)
        {
        }
    };

    template <typename DataType_>
    ConstElementIterator<storage::Dense, container::MatrixTile, DataType_>::ConstElementIterator(const DenseMatrixTile<DataType_> & matrix, unsigned long index) :
        PrivateImplementationPattern<ConstElementIterator<storage::Dense, container::MatrixTile, DataType_>, Single>(
                new Implementation<ConstElementIterator<storage::Dense, container::MatrixTile, DataType_> >(matrix, index))
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::Dense, container::MatrixTile, DataType_>::ConstElementIterator(const ConstElementIterator & other) :
        PrivateImplementationPattern<ConstElementIterator<storage::Dense, container::MatrixTile, DataType_>, Single>(
                new Implementation<ConstElementIterator<storage::Dense, container::MatrixTile, DataType_> >(*other._imp))
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::Dense, container::MatrixTile, DataType_>::ConstElementIterator(
            const ElementIterator<storage::Dense, container::MatrixTile, DataType_> & other) :
        PrivateImplementationPattern<ConstElementIterator<storage::Dense, container::MatrixTile, DataType_>, Single>(
                new Implementation<ConstElementIterator<storage::Dense, container::MatrixTile, DataType_> >(*other._imp))
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::Dense, container::MatrixTile, DataType_>::~ConstElementIterator()
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::Dense, container::MatrixTile, DataType_> &
    ConstElementIterator<storage::Dense, container::MatrixTile, DataType_>::operator= (
            const ConstElementIterator<storage::Dense, container::MatrixTile, DataType_> & other)
    {
        if (&other == this)
            return *this;

        this->_imp->_tile = other._imp->_tile;
        this->_imp->_index = other._imp->_index;
        this->_imp->_pos = other._imp->_pos;

        return *this;
    }

    template <typename DataType_>
    ConstElementIterator<storage::Dense, container::MatrixTile, DataType_> &
    ConstElementIterator<storage::Dense, container::MatrixTile, DataType_>::operator++ ()
    {
        CONTEXT("When incrementing ConstElementIterator<Dense, MatrixTile> by one:");

        ++this->_imp->_index;
        if (this->_imp->_index % this->_imp->_tile._columns)
        {
            ++this->_imp->_pos;
        }
        else
        {
            this->_imp->_pos += this->_imp->_tile._source_columns - this->_imp->_tile._columns + 1;
        }

        return *this;
    }

    template <typename DataType_>
    ConstElementIterator<storage::Dense, container::MatrixTile, DataType_> &
    ConstElementIterator<storage::Dense, container::MatrixTile, DataType_>::operator+= (const unsigned long step)
    {
        CONTEXT("When incrementing ConstElementIterator<Dense, MatrixTile> by '" + stringify(step) + "':");

        this->_imp->_index += step;
        this->_imp->_pos = this->_imp->_tile._calc_pos(this->_imp->_index);

        return *this;
    }

    template <typename DataType_>
    const DataType_ &
    ConstElementIterator<storage::Dense, container::MatrixTile, DataType_>::operator* () const
    {
        CONTEXT("When accessing unassignable element at index '" + stringify(this->_imp->_index) + "':");

        return this->_imp->_tile._elements[this->_imp->_pos];
    }

    template <typename DataType_>
    bool
    ConstElementIterator<storage::Dense, container::MatrixTile, DataType_>::operator< (
            const ConstElementIterator<storage::Dense, container::MatrixTile, DataType_> & other) const
    {
        return this->_imp->_index < other._imp->_index;
    }

    template <typename DataType_>
    bool
    ConstElementIterator<storage::Dense, container::MatrixTile, DataType_>::operator== (
            const ConstElementIterator<storage::Dense, container::MatrixTile, DataType_> & other) const
    {
        return ((this->_imp->_tile._elements.get() == other._imp->_tile._elements.get()) && (this->_imp->_index == other._imp->_index));
    }

    template <typename DataType_>
    bool
    ConstElementIterator<storage::Dense, container::MatrixTile, DataType_>::operator!= (
            const ConstElementIterator<storage::Dense, container::MatrixTile, DataType_> & other) const
    {
        return ((this->_imp->_tile._elements.get() != other._imp->_tile._elements.get()) || (this->_imp->_index != other._imp->_index));
    }

    template <typename DataType_>
    unsigned long
    ConstElementIterator<storage::Dense, container::MatrixTile, DataType_>::column() const
    {
        return this->_imp->_index % this->_imp->_tile._columns;
    }

    template <typename DataType_>
    unsigned long
    ConstElementIterator<storage::Dense, container::MatrixTile, DataType_>::row() const
    {
        return this->_imp->_index / this->_imp->_tile._columns;
    }

    template <typename DataType_>
    unsigned long
    ConstElementIterator<storage::Dense, container::MatrixTile, DataType_>::index() const
    {
        return this->_imp->_index;
    }

    template <typename DataType_> struct Implementation<ElementIterator<storage::Dense, container::MatrixTile, DataType_> >
    {
            /// Our parent matrix.
            DenseMatrixTile<DataType_> _tile;

            /// Our index.
            unsigned long _index;

            /// Our position inside the matrix's elements.
            unsigned long _pos;

        Implementation(DenseMatrixTile<DataType_> & tile, unsigned long index) :
            _tile(tile),
            _index(index),
            _pos(tile._calc_pos(index))
        {
        }
    };

    template <typename DataType_>
    ElementIterator<storage::Dense, container::MatrixTile, DataType_>::ElementIterator(DenseMatrixTile<DataType_> & matrix, unsigned long index) :
        PrivateImplementationPattern<ElementIterator<storage::Dense, container::MatrixTile, DataType_>, Single>(
                new Implementation<ElementIterator<storage::Dense, container::MatrixTile, DataType_> >(matrix, index))
    {
    }

    template <typename DataType_>
    ElementIterator<storage::Dense, container::MatrixTile, DataType_>::ElementIterator(const ElementIterator & other) :
        PrivateImplementationPattern<ElementIterator<storage::Dense, container::MatrixTile, DataType_>, Single>(
                new Implementation<ElementIterator<storage::Dense, container::MatrixTile, DataType_> >(*other._imp))
    {
    }

    template <typename DataType_>
    ElementIterator<storage::Dense, container::MatrixTile, DataType_>::~ElementIterator()
    {
    }

    template <typename DataType_>
    ElementIterator<storage::Dense, container::MatrixTile, DataType_> &
    ElementIterator<storage::Dense, container::MatrixTile, DataType_>::operator= (
            const ElementIterator<storage::Dense, container::MatrixTile, DataType_> & other)
    {
        if (&other == this)
            return *this;

        this->_imp->_tile = other._imp->_tile;
        this->_imp->_index = other._imp->_index;
        this->_imp->_pos = other._imp->_pos;

        return *this;
    }

    template <typename DataType_>
    ElementIterator<storage::Dense, container::MatrixTile, DataType_> &
    ElementIterator<storage::Dense, container::MatrixTile, DataType_>::operator++ ()
    {
        CONTEXT("When incrementing ElementIterator<Dense, MatrixTile> by one:");

        ++this->_imp->_index;
        if (this->_imp->_index % this->_imp->_tile._columns)
        {
            ++this->_imp->_pos;
        }
        else
        {
            this->_imp->_pos += this->_imp->_tile._source_columns - this->_imp->_tile._columns + 1;
        }

        return *this;
    }

    template <typename DataType_>
    ElementIterator<storage::Dense, container::MatrixTile, DataType_> &
    ElementIterator<storage::Dense, container::MatrixTile, DataType_>::operator+= (const unsigned long step)
    {
        CONTEXT("When incrementing ElementIterator<Dense, MatrixTile> by '" + stringify(step) + "':");

        this->_imp->_index += step;
        this->_imp->_pos = this->_imp->_tile._calc_pos(this->_imp->_index);

        return *this;
    }

    template <typename DataType_>
    DataType_ &
    ElementIterator<storage::Dense, container::MatrixTile, DataType_>::operator* () const
    {
        CONTEXT("When accessing unassignable element at index '" + stringify(this->_imp->_index) + "':");

        return this->_imp->_tile._elements[this->_imp->_pos];
    }

    template <typename DataType_>
    bool
    ElementIterator<storage::Dense, container::MatrixTile, DataType_>::operator< (
            const ElementIterator<storage::Dense, container::MatrixTile, DataType_> & other) const
    {
        return this->_imp->_index < other._imp->_index;
    }

    template <typename DataType_>
    bool
    ElementIterator<storage::Dense, container::MatrixTile, DataType_>::operator== (
            const ElementIterator<storage::Dense, container::MatrixTile, DataType_> & other) const
    {
        return ((this->_imp->_tile._elements.get() == other._imp->_tile._elements.get()) && (this->_imp->_index == other._imp->_index));
    }

    template <typename DataType_>
    bool
    ElementIterator<storage::Dense, container::MatrixTile, DataType_>::operator!= (
            const ElementIterator<storage::Dense, container::MatrixTile, DataType_> & other) const
    {
        return ((this->_imp->_tile._elements.get() != other._imp->_tile._elements.get()) || (this->_imp->_index != other._imp->_index));
    }

    template <typename DataType_>
    unsigned long
    ElementIterator<storage::Dense, container::MatrixTile, DataType_>::column() const
    {
        return this->_imp->_index % this->_imp->_tile._columns;
    }

    template <typename DataType_>
    unsigned long
    ElementIterator<storage::Dense, container::MatrixTile, DataType_>::row() const
    {
        return this->_imp->_index / this->_imp->_tile._columns;
    }

    template <typename DataType_>
    unsigned long
    ElementIterator<storage::Dense, container::MatrixTile, DataType_>::index() const
    {
        return this->_imp->_index;
    }

    template <typename DataType_>
    std::ostream &
    operator<< (std::ostream & lhs, const DenseMatrixTile<DataType_> & b)
    {
        lhs << "[" << std::endl;
        for (typename DenseMatrixTile<DataType_>::ConstElementIterator i(b.begin_elements()), i_end(b.end_elements()) ;
                i != i_end ; ++i)
        {
            lhs << "  " << *i;
        }
        lhs << "]" << std::endl;

        return lhs;
    }


    template <typename DataType_>
    bool
    operator== (const DenseMatrixTile<DataType_> & a, const DenseMatrixTile<DataType_> & b)
    {
        CONTEXT("When comparing two tiles of dense matrices:");

        if (a.columns() != b.columns())
        {
            throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
        }

        if (a.rows() != b.rows())
        {
            throw MatrixRowsDoNotMatch(b.rows(), a.rows());
        }

        for (typename DenseMatrixTile<DataType_>::ConstElementIterator i(a.begin_elements()), i_end(a.end_elements()), j(b.begin_elements()) ;
                i != i_end ; ++i, ++j)
        {
            if (std::fabs(*i - *j) > std::numeric_limits<DataType_>::epsilon())
                return false;
        }

        return true;
    }
}

#endif
