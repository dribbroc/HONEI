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

#ifndef LIBLA_GUARD_DIAGONAL_MATRIX_HH
#define LIBLA_GUARD_DIAGONAL_MATRIX_HH 1

#include <libla/element_iterator.hh>
#include <libla/matrix.hh>
#include <libla/vector.hh>
#include <libutil/shared_array.hh>

#include <iostream>

#include <string>
#include <vector>

namespace pg512 ///< \todo Namespace name?
{
    /**
     * A DiagonalMatrix is a square matrix whose diagonal elements are the only non-zero elements.
     **/
    template <typename DataType_> class DiagonalMatrix :
        public Matrix<DataType_>
    {
        private:
            /// Pointer to our diagonal vector.
            DenseVector<DataType_> _diagonal;

            /// Our columns.
            unsigned long _size;

            /// Our zero-element.
            static DataType_ _zero_element;

            /// Our row-vectors.
            std::vector<std::tr1::shared_ptr<SparseVector<DataType_> > > _row_vectors;

        public:
            /// Our implementation of ElementIterator.
            template <typename ElementType_> class ElementIteratorImpl;
            friend class ElementIteratorImpl<DataType_>;

            /// Type of the iterator over our elements.
            typedef ElementIteratorWrapper<Matrix<DataType_>, DataType_> ElementIterator;

            /// Constructor.
            DiagonalMatrix(const DenseVector<DataType_> & vector) :
                _diagonal(vector),
                _size(vector.size())
            {
            }

            /// Returns iterator pointing to the first element of the matrix.
            virtual ElementIterator begin_elements() const
            {
                return ElementIterator(new ElementIteratorImpl<DataType_>(*this, 0));
            }

            /// Returns iterator pointing behind the last element of the vector.
            virtual ElementIterator end_elements() const
            {
                return ElementIterator(new ElementIteratorImpl<DataType_>(*this, _size));
            }

            /// Returns our columns.
            virtual unsigned long columns() const
            {
                return _size;
            }

            /// Returns our rows.
            virtual unsigned long rows() const
            {
                return _size;
            }

            /// Retrieves row-vector by row-index, zero-based, unassignable.
            virtual const Vector<DataType_> & operator[] (unsigned long row) const
            {
                if (! _row_vectors[row])
                {
                    _row_vectors[row].reset(new SparseVector<DataType_>(_size, 1));
                    (*_row_vectors[row])[row] = _diagonal[row];
                }

                return *_row_vectors[row];
            }

            /// Returns our diagonal as a Vector.
            DenseVector<DataType_> & diagonal()
            {
                return _diagonal;
            }
    };

    /**
     * A DiagonalMatrix::ElementIteratorImpl is an iterator implementation for diagonal matrices.
     **/
    template <> template <typename DataType_> class DiagonalMatrix<DataType_>::ElementIteratorImpl<DataType_> :
        public ElementIteratorImplBase<Matrix<DataType_>, DataType_>
    {
        private:
            /// Our matrix.
            const DiagonalMatrix<DataType_> & _matrix;

            /// Our index.
            unsigned long _index;

        public:
            /// Constructor.
            ElementIteratorImpl(const DiagonalMatrix<DataType_> & matrix, unsigned long index) :
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
            virtual const DataType_ & operator* ()
            {
                std::cout << "c = " << column() << ", r = " << row() << std::endl;
                if (column() == row())
                    return _matrix._diagonal[row()];
                else
                    return _zero_element;
            }

            /// Our index.
            virtual const unsigned long index() const
            {
                return _index;
            }

            /// Our column.
            virtual const unsigned long column() const
            {
                return _index % _matrix._size;
            }

            /// Our row.
            virtual const unsigned long row() const
            {
                return _index / _matrix._size;
            }

            /// Our parent.
            virtual const Matrix<DataType_> * parent() const
            {
                return &_matrix;
            }
    };

}

#endif
