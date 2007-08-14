/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_SPARSE_MATRIX_HH
#define LIBLA_GUARD_SPARSE_MATRIX_HH 1

#include <libla/element_iterator.hh>
#include <libla/matrix.hh>
#include <libla/sparse_vector.hh>
#include <libutil/shared_array.hh>

#include <iterator>
#include <vector>

/**
 * \file
 *
 * Implementation of SparseMatrix and related classes.
 *
 * \ingroup grpmatrix
 **/
namespace pg512 ///< \todo Namespace name?
{
    /**
     * \brief SparseMatrix is a matrix with O(rows) sparse row-vectors.
     *
     * \ingroup grpmatrix
     **/
    template <typename DataType_> class SparseMatrix :
        public RowAccessMatrix<DataType_>,
        public MutableMatrix<DataType_>
    {
        private:
            /// Our row-vectors' initial capacity.
            unsigned long _capacity;

            /// Our columns.
            unsigned long _columns;

            /// Our rows.
            unsigned long _rows;

            /// Our row-vectors.
            SharedArray<std::tr1::shared_ptr<SparseVector<DataType_> > > _row_vectors;

            /// Our zero-vector.
            const SparseVector<DataType_> _zero_vector;

            /// Our implementation of ElementIteratorBase.
            template <typename ElementType_> class SparseElementIterator;
            /// Our smart implementation of ElementIteratorBase.
            template <typename ElementType_> class NonZeroElementIterator;            

            typedef typename Matrix<DataType_>::MatrixElementIterator MatrixElementIterator;

        public:
            friend class SparseElementIterator<DataType_>;
            friend class NonZeroElementIterator<DataType_>;

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
             **/
            SparseMatrix(unsigned long columns, unsigned long rows, unsigned long capacity = 1) :
                _capacity(capacity),
                _columns(columns),
                _rows(rows),
                _row_vectors(new std::tr1::shared_ptr<SparseVector<DataType_> >[rows + 1]),
                _zero_vector(columns, 1)
            {
            }

            ~SparseMatrix()
            {
            }

            /// \}

            /// Returns iterator pointing to the first element of the matrix.
            virtual ConstElementIterator begin_elements() const
            {
                return ConstElementIterator(new SparseElementIterator<DataType_>(*this, 0));
            }

            /// Returns iterator pointing behind the last element of the matrix.
            virtual ConstElementIterator end_elements() const
            {
                return ConstElementIterator(new SparseElementIterator<DataType_>(*this, _rows * _columns));
            }

            /// Returns iterator pointing to the first element of the matrix.
            virtual ElementIterator begin_elements()
            {
                return ElementIterator(new SparseElementIterator<DataType_>(*this, 0));
            }

            /// Returns iterator pointing behind the last element of the matrix.
            virtual ElementIterator end_elements()
            {
                return ElementIterator(new SparseElementIterator<DataType_>(*this, _rows * _columns));
            }

            /// Returns const iterator pointing to the first non-zero element of the matrix.
            virtual ConstElementIterator begin_non_zero_elements() const
            {
                return ConstElementIterator(new NonZeroElementIterator<DataType_>(*this, 0));
            }

            /// Returns const iterator pointing behind the last element of the matrix.
            virtual ConstElementIterator end_non_zero_elements() const
            {
                return ConstElementIterator(new NonZeroElementIterator<DataType_>(*this, _rows * _columns));
            }

            /// Returns iterator pointing to the first non-zero element of the matrix.
            virtual ElementIterator begin_non_zero_elements()
            {
                return ElementIterator(new NonZeroElementIterator<DataType_>(*this, 0));
            }

            /// Returns iterator pointing behind the last element of the matrix.
            virtual ElementIterator end_non_zero_elements()
            {
                return ElementIterator(new NonZeroElementIterator<DataType_>(*this, _rows * _columns));
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
            virtual const SparseVector<DataType_> & operator[] (unsigned long row) const
            {
                if (! _row_vectors[row])
                    return _zero_vector;

                return *_row_vectors[row];
            }

            /// Retrieves row vector by index, zero-based, assignable.
            virtual SparseVector<DataType_> & operator[] (unsigned long row)
            {
                if (! _row_vectors[row])
                    _row_vectors[row].reset(new SparseVector<DataType_>(_columns, _capacity));

                return *_row_vectors[row];
            }

            /// Returns a copy of the matrix.
            virtual SparseMatrix * copy() const
            {
                SparseMatrix * result(new SparseMatrix(_columns, _rows, _capacity));

                for (unsigned long i(0) ; i < _rows ; ++i)
                {
                    if (_row_vectors[i])
                        result->_row_vectors[i].reset(_row_vectors[i]->copy());
                }

                return result;
            }
    };

    /**
     * \brief SparseMatrix::SparseElementIterator is a simple iterator implementation for sparse matrices.
     *
     * \ingroup grpmatrix
     **/
    template <> template <typename DataType_> class SparseMatrix<DataType_>::SparseElementIterator<DataType_> :
        public MatrixElementIterator
    {
        private:
            /// Our parent matrix.
            const SparseMatrix<DataType_> & _matrix;

            /// Our index.
            unsigned long _index;

            /// Our row-vector's iterator.
            typename Vector<DataType_>::ElementIterator _iter;

            /// Our row index.
            unsigned long _row;

            static typename Vector<DataType_>::ElementIterator _get_iterator(const SparseMatrix<DataType_> & matrix,
                    unsigned long index)
            {
                if (! matrix._row_vectors[index / matrix._columns])
                    matrix._row_vectors[index / matrix._columns].reset(new SparseVector<DataType_>(matrix._columns,
                                matrix._capacity));

                return matrix._row_vectors[index / matrix._columns]->begin_elements();
            }

        public:
            /// \name Constructors and destructor
            /// \{

            /**
             * Constructor.
             *
             * \param matrix The parent matrix that is referenced by the iterator.
             * \param index The index into the matrix.
             **/
            SparseElementIterator(const SparseMatrix<DataType_> & matrix, unsigned long index) :
                _matrix(matrix),
                _index(index),
                _iter(_get_iterator(matrix, index)),
                _row(index / matrix._columns)
            {
            }

            /// Copy-constructor.
            SparseElementIterator(SparseElementIterator<DataType_> const & other) :
                _matrix(other._matrix),
                _index(other._index),
                _iter(other._iter),
                _row(other._row)
            {
            }

            /// Destructor.
            virtual ~SparseElementIterator()
            {
            }

            /// \}

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            virtual SparseElementIterator<DataType_> & operator++ ()
            {
                ++_index;
                unsigned long row(_index / _matrix._columns);
                if (row != _row)
                {
                    if (! _matrix._row_vectors[row])
                        _matrix._row_vectors[row].reset(new SparseVector<DataType_>(_matrix._columns,
                                    _matrix._capacity));

                    _iter = _matrix._row_vectors[row]->begin_elements();
                    _row = row;
                }
                else
                {
                    ++_iter;
                }

                return *this;
            }

            /// Dereference operator that returns assignable reference.
            virtual DataType_ & operator* ()
            {
                return *_iter;
            }

            /// Dereference operator that returns umassignable reference.
            virtual const DataType_ & operator* () const
            {
                typename Vector<DataType_>::ConstElementIterator const_iter(_iter);

                return *const_iter;
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
                return _iter.index();
            }

            /// Returns our row index.
            virtual unsigned long row() const
            {
                return _row;
            }

            /// Returns a pointer to our parent container.
            virtual const Matrix<DataType_> * parent() const
            {
                return &_matrix;
            }

            /// \}
    };

    /**
     * \brief SparseMatrix::NonZeroElementIterator is a smart iterator implementation that iterates over non-zero
     * \brief elements of sparse matrices.
     *
     * \ingroup grpmatrix
     **/
    template <> template <typename DataType_> class SparseMatrix<DataType_>::NonZeroElementIterator<DataType_> :
        public MatrixElementIterator
    {
        private:
            /// Our parent matrix.
            const SparseMatrix<DataType_> & _matrix;

            /// Our index.
            unsigned long _index;

            /// Our row-vector's iterator.
            typename Vector<DataType_>::ElementIterator _iter;

            /// Our row index.
            unsigned long _row;

            static typename Vector<DataType_>::ElementIterator _get_iterator(const SparseMatrix<DataType_> & matrix,
                    unsigned long index, unsigned long row)
            {
                if (! matrix._row_vectors[row])
                    matrix._row_vectors[row].reset(new SparseVector<DataType_>(matrix._columns,
                                matrix._capacity));                       
                return matrix._row_vectors[row]->begin_non_zero_elements();
            }

            static typename Vector<DataType_>::ElementIterator _get_first_iterator(const SparseMatrix<DataType_> & matrix, 
                    unsigned long row)
            {
                while (! matrix._row_vectors[row])
                {
                    ++row;
                    if (matrix._rows == row)
                    {
                        matrix._row_vectors[row].reset(new SparseVector<DataType_>(matrix._columns,
                                matrix._capacity));
                        return matrix._row_vectors[row]->begin_non_zero_elements();
                    }
                }
                return matrix._row_vectors[row]->begin_non_zero_elements();
            }
        public:
            /// \name Constructors and destructor
            /// \{

            /**
             * Constructor.
             *
             * \param matrix The parent matrix that is referenced by the iterator.
             * \param index The index into the matrix.
             **/
            NonZeroElementIterator(const SparseMatrix<DataType_> & matrix, unsigned long index) :
                _matrix(matrix),
                _index(index),
                _row(index / matrix._columns),
                _iter(_get_first_iterator(matrix, _row))
            {
            
            }

            /// Copy-constructor.
            NonZeroElementIterator(NonZeroElementIterator<DataType_> const & other) :
                _matrix(other._matrix),
                _index(other._index),
                _iter(other._iter),
                _row(other._row)
            {
            }

            /// Destructor.
            virtual ~NonZeroElementIterator()
            {
            }

            /// \}

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            virtual NonZeroElementIterator<DataType_> & operator++ ()
            {
                ++_index;
                ++_iter;
                if (_iter != _matrix._row_vectors[_row]->end_non_zero_elements())
                {
                    return *this; 
                }
                ++ _row;
                while (!_matrix._row_vectors[_row])
                {
                    if (_matrix._rows == _row)
                    {
                        _matrix._row_vectors[_row].reset(new SparseVector<DataType_>(_matrix._columns,
                                _matrix._capacity));
                        _iter = _matrix._row_vectors[_row]->begin_non_zero_elements();
                        return *this;
                    }
                    ++ _row;
                }
                _iter = _matrix._row_vectors[_row]->begin_non_zero_elements();
                return *this;
            }

            /// Dereference operator that returns assignable reference.
            virtual DataType_ & operator* ()
            {
                return *_iter;
            }

            /// Dereference operator that returns umassignable reference.
            virtual const DataType_ & operator* () const
            {
                typename Vector<DataType_>::ConstElementIterator const_iter(_iter);

                return *const_iter;
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
                return _iter.index();
            }

            /// Returns our row index.
            virtual unsigned long row() const
            {
                return _row;
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
