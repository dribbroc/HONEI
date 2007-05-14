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

#ifndef LIBLA_GUARD_MATRIX_HH
#define LIBLA_GUARD_MATRIX_HH 1

#include <iostream>

#include <libla/element_iterator.hh>
#include <libla/vector.hh>
#include <libutil/shared_array.hh>

#include <string.h>
#include <iterator>

namespace pg512 ///< \todo Namespace name?
{
    /**
     * A Matrix is the abstract baseclass for all matrix-like types used.
     **/
    template <typename DataType_> class Matrix
    {
        public:
            /// Type of the iterator over our elements.
            typedef ElementIteratorWrapper<Matrix<DataType_>, DataType_> ElementIterator;

            /// Type of the iterator over our row/column/band/diagonal vectors.
//            typedef VectorIteratorBase<Matrix<DataType_>, Vector<DataType_> > VectorIterator;

            /// Returns iterator pointing to the first element of the matrix.
            virtual ElementIterator begin_elements() const = 0;

            /// Returns iterator pointing behind the last element of the matrix.
            virtual ElementIterator end_elements() const = 0;
#if 0
            /// Returns iterator pointing to the first column of the matrix.
            virtual VectorIterator begin_columns() const = 0;

            /// Returns iterator pointing behind the last column of the matrix.
            virtual VectorIterator end_columns() const = 0;

            /// Returns iterator pointing to the first row of the matrix.
            virtual VectorIterator begin_rows() const = 0;

            /// Returns iterator pointing behind the last row of the matrix.
            virtual VectorIterator end_rows() const = 0;
#endif
            /// Returns our number of columns.
            virtual unsigned long columns() const = 0;

            /// Returns our number of rows.
            virtual unsigned long rows() const = 0;

            /// Retrieves element by index, zero-based, unassignable.
            virtual const Vector<DataType_> & operator[] (unsigned long row) const = 0;

            /// Retrieves element by index, zero-based, assignable
            virtual Vector<DataType_> & operator[] (unsigned long row) = 0;
    };

    /// Output our Matrix to an ostream.
    template <typename DataType_> std::ostream & operator<< (std::ostream & lhs, const Matrix<DataType_> & m)
    {
        lhs << "[ " << std::endl;
        for (unsigned long r(0) ; r < m.rows() ; ++r)
        {
            const Vector<DataType_> & v(m[r]);

            for (unsigned long c(0) ; c < m.columns() ; ++c)
            {
                lhs << v[c] << " ";
            }

            lhs << std::endl;
        }
        lhs << "]";

        return lhs;
    }

    /**
     * A DenseMatrix is a matrix with O(column * row) non-zero elements which keeps its data
     * sequential.
     **/
    template <typename DataType_> class DenseMatrix :
        public Matrix<DataType_>
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

        public:
            /// Our implementation of ElementIterator.
            template <typename ElementType_> class ElementIteratorImpl;
            friend class ElementIteratorImpl<DataType_>;

            /// Type of the iterator over our elements.
            typedef ElementIteratorWrapper<Matrix<DataType_>, DataType_> ElementIterator;

            /// Type of the iterator over our row/column/band/diagonal vectors.
//            typedef VectorIteratorBase<Matrix<DataType_>, Vector<DataType_> > VectorIterator;

            /// Constructor.
            DenseMatrix(unsigned long columns, unsigned long rows) :
                _elements(new DataType_[rows * columns]),
                _columns(columns),
                _rows(rows),
                _row_vectors(new std::tr1::shared_ptr<DenseVector<DataType_> >[rows])
            {
            }

            /// Constructor.
            DenseMatrix(unsigned long columns, unsigned long rows, DataType_ value) :
                _elements(new DataType_[rows * columns]),
                _columns(columns),
                _rows(rows),
                _row_vectors(new std::tr1::shared_ptr<DenseVector<DataType_> >[rows])
            {
                for (unsigned long i(0) ; i < rows * columns ; ++i)
                    _elements[i] = value;
            }

            /// Returns iterator pointing to the first element of the matrix.
            virtual ElementIterator begin_elements() const
            {
                return ElementIterator(new ElementIteratorImpl<DataType_>(*this, 0));
            }

            /// Returns iterator pointing behind the last element of the vector.
            virtual ElementIterator end_elements() const
            {
                return ElementIterator(new ElementIteratorImpl<DataType_>(*this, _rows * _columns));
            }

            /// Returns our columns.
            virtual unsigned long columns() const
            {
                return _columns;
            }

            /// Returns our rows.
            virtual unsigned long rows() const
            {
                return _rows;
            }

            /// Retrieves element by index, zero-based, assignable.
            virtual const Vector<DataType_> & operator[] (unsigned long row) const
            {
                if (! _row_vectors[row])
                    _row_vectors[row].reset(new DenseVector<DataType_>(_columns, _elements, row * _columns, 1));

                return *_row_vectors[row];
            }

            /// Retrieves element by index, zero-based, assignable.
            virtual Vector<DataType_> & operator[] (unsigned long row)
            {
                if (! _row_vectors[row])
                    _row_vectors[row].reset(new DenseVector<DataType_>(_columns, _elements, row * _columns, 1));

                return *_row_vectors[row];
            }
    };

    /**
     * A DenseMatrix::ElementIteratorImpl is a simple iterator implementation for dense matrices.
     **/
    template <> template <typename DataType_> class DenseMatrix<DataType_>::ElementIteratorImpl<DataType_> :
        public ElementIteratorImplBase<Matrix<DataType_>, DataType_>
    {
        private:
            /// Our matrix.
            const DenseMatrix<DataType_> & _matrix;

            /// Our index.
            unsigned long _index;

        public:
            /// Constructor.
            ElementIteratorImpl(const DenseMatrix<DataType_> & matrix, unsigned long index) :
                _matrix(matrix),
                _index(index)
            {
            }

            /// Constructor.
            ElementIteratorImpl(ElementIteratorImpl<DataType_> const & other) :
                _matrix(other._matrix),
                _index(other._index)
            {
            }

            /// Preincrement operator.
            virtual ElementIteratorImpl<DataType_> & operator++ ()
            {
                ++_index;
                return *this;
            }

            /// Postincrement operator.
            virtual ElementIteratorImpl<DataType_> operator++ (int)
            {
                ElementIteratorImpl<DataType_> result(*this);
                ++_index;
                return result;
            }

            /// Equality operator.
            virtual bool operator== (const ElementIteratorImplBase<Matrix<DataType_>, DataType_> & other) const
            {
                return (&_matrix == other.parent()) && (_index == other.index());
            }

            /// Inequality operator.
            virtual bool operator!= (const ElementIteratorImplBase<Matrix<DataType_>, DataType_> & other) const
            {
                return ((&_matrix != other.parent()) || (_index != other.index()));
            }

            /// Dereference operator 
            virtual DataType_ & operator* () const
            {
                return _matrix._elements[_index];
            }

            /// Our index.
            virtual const unsigned long index() const
            {
                return _index;
            }

            /// Our column.
            virtual const unsigned long column() const
            {
                return _index % _matrix._rows;
            }

            /// Our row.
            virtual const unsigned long row() const
            {
                return _index / _matrix._rows;
            }

            /// Our parent.
            virtual const Matrix<DataType_> * parent() const
            {
                return &_matrix;
            }
    };
}

#endif
