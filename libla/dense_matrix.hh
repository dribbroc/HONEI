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
#include <libla/dense_vector.hh>
#include <libutil/shared_array.hh>

#include <string.h>
#include <iterator>

/**
 * \file
 *
 * Implementation of DenseMatrix and related classes.
 *
 * \ingroup grpmatrix
 **/
namespace pg512 ///< \todo Namespace name?
{
    /**
     * A DenseMatrix is a matrix with O(column * row) non-zero elements which keeps its data
     * sequential.
     *
     * \ingroup grpmatrix
     **/
    template <typename DataType_> class DenseMatrix :
        public Matrix<DataType_>,
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

        public:
            /// Our implementation of ElementIterator.
            template <typename ElementType_> class ElementIteratorImpl;
            friend class ElementIteratorImpl<DataType_>;
            friend class ElementIteratorImpl<const DataType_>;

            /// Type of the const iterator over our elements.
            typedef ElementIteratorWrapper<Matrix<DataType_>, DataType_, const DataType_> ConstElementIterator;

            /// Type of the const iterator over our vectors.
            typedef VectorIteratorWrapper<DataType_, const DataType_> ConstVectorIterator;

            /// Type of the iterator over our elements.
            typedef ElementIteratorWrapper<Matrix<DataType_>, DataType_> ElementIterator;

            /// Type of the iterator over our vectors.
            typedef VectorIteratorWrapper<DataType_, DataType_> VectorIterator;

            /**
             * Constructor.
             *
             * \param columns Number of columns of the new matrix.
             * \param rows Number of rows of the new matrix.
             **/
            DenseMatrix(unsigned long columns, unsigned long rows) :
                _elements(new DataType_[rows * columns]),
                _columns(columns),
                _rows(rows),
                _row_vectors(new std::tr1::shared_ptr<DenseVector<DataType_> >[rows])
            {
            }

            /**
             * Constructor.
             *
             * \param columns Number of columns of the new matrix.
             * \param rows Number of rows of the new matrix.
             * \param value Default value of each of the new matrice's elements.
             **/
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
            virtual ConstElementIterator begin_elements() const
            {
                return ConstElementIterator(new ElementIteratorImpl<const DataType_>(*this, 0));
            }

            /// Returns iterator pointing behind the last element of the matrix.
            virtual ConstElementIterator end_elements() const
            {
                return ConstElementIterator(new ElementIteratorImpl<const DataType_>(*this, _rows * _columns));
            }

            /// Returns iterator pointing to the first element of the matrix.
            virtual ElementIterator begin_elements()
            {
                return ElementIterator(new ElementIteratorImpl<DataType_>(*this, 0));
            }

            /// Returns iterator pointing behind the last element of the matrix.
            virtual ElementIterator end_elements()
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

            /// Retrieves element by index, zero-based, unassignable.
            virtual const Vector<DataType_> & operator[] (unsigned long row) const
            {
                if (! _row_vectors[row])
                    _row_vectors[row].reset(new DenseVector<DataType_>(_columns, _elements, (row-1) * _columns, 1));

                return *_row_vectors[row];
            }

            /// Retrieves element by index, zero-based, assignable.
            virtual Vector<DataType_> & operator[] (unsigned long row)
            {
                if (! _row_vectors[row])
                    _row_vectors[row].reset(new DenseVector<DataType_>(_columns, _elements, (row-1) * _columns, 1));

                return *_row_vectors[row];
            }

            /// Retrieves element by index, zero-based, unassignable.
            virtual const Vector<DataType_> & column(unsigned long column) const
            {
                if (! _column_vectors[column])
                    _column_vectors[column].reset(new DenseVector<DataType_>(_rows, _elements, (column-1), _columns));

                return *_column_vectors[column];
            }

            /// Retrieves element by index, zero-based, assignable.
            virtual Vector<DataType_> & column(unsigned long column)
            {
                if (! _column_vectors[column])
                    _column_vectors[column].reset(new DenseVector<DataType_>(_rows, _elements, (column-1), _columns));

                return *_column_vectors[column];
            }
    };

    /**
     * A DenseMatrix::ElementIteratorImpl is a simple iterator implementation for dense matrices.
     *
     * \ingroup grpmatrix
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
            /**
             * Constructor.
             *
             * \param matrix The parent matrix that is referenced by the iterator.
             * \param index The index into the matrix.
             **/
            ElementIteratorImpl(const DenseMatrix<DataType_> & matrix, unsigned long index) :
                _matrix(matrix),
                _index(index)
            {
            }

            /// Copy-constructor.
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

            /// Dereference operator.
            virtual DataType_ & operator* () const
            {
                return _matrix._elements[_index];
            }

            /// Returns our index.
            virtual const unsigned long index() const
            {
                return _index;
            }

            /// Returns our column.
            virtual const unsigned long column() const
            {
                return _index % _matrix._columns;
            }

            /// Returns our row.
            virtual const unsigned long row() const
            {
                return _index / _matrix._columns;
            }

            /// Returns our parent matrix.
            virtual const Matrix<DataType_> * parent() const
            {
                return &_matrix;
            }
    };

    template <> template <typename DataType_> class DenseMatrix<DataType_>::ElementIteratorImpl<const DataType_> :
        public ElementIteratorImplBase<Matrix<DataType_>, DataType_, const DataType_>
    {
        private:
            /// Our matrix.
            const DenseMatrix<DataType_> & _matrix;

            /// Our index.
            unsigned long _index;

        public:
            /**
             * Constructor.
             *
             * \param matrix The parent matrix that is referenced by the iterator.
             * \param index The index into the matrix.
             **/
            ElementIteratorImpl(const DenseMatrix<DataType_> & matrix, unsigned long index) :
                _matrix(matrix),
                _index(index)
            {
            }

            /// Copy-constructor.
            ElementIteratorImpl(ElementIteratorImpl<const DataType_> const & other) :
                _matrix(other._matrix),
                _index(other._index)
            {
            }

            /// Preincrement operator.
            virtual ElementIteratorImpl<const DataType_> & operator++ ()
            {
                ++_index;
                return *this;
            }

            /// Postincrement operator.
            virtual ElementIteratorImpl<const DataType_> operator++ (int)
            {
                ElementIteratorImpl<const DataType_> result(*this);
                ++_index;
                return result;
            }

            /// Equality operator.
            virtual bool operator== (const ElementIteratorImplBase<Matrix<DataType_>, DataType_, const DataType_> & other) const
            {
                return (&_matrix == other.parent()) && (_index == other.index());
            }

            /// Inequality operator.
            virtual bool operator!= (const ElementIteratorImplBase<Matrix<DataType_>, DataType_, const DataType_> & other) const
            {
                return ((&_matrix != other.parent()) || (_index != other.index()));
            }

            /// Dereference operator.
            virtual const DataType_ & operator* () const
            {
                return _matrix._elements[_index];
            }

            /// Returns our index.
            virtual const unsigned long index() const
            {
                return _index;
            }

            /// Returns our column.
            virtual const unsigned long column() const
            {
                return _index % _matrix._columns;
            }

            /// Returns our row.
            virtual const unsigned long row() const
            {
                return _index / _matrix._columns;
            }

            /// Returns our parent matrix.
            virtual const Matrix<DataType_> * parent() const
            {
                return &_matrix;
            }
    };

}

#endif
