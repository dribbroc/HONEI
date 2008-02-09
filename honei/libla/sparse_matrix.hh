/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Sven Mallach <sven.mallach@honei.org>
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

#include <honei/libla/element_iterator.hh>
#include <honei/libla/vector_iterator.hh>
#include <honei/libla/matrix.hh>
#include <honei/libla/sparse_vector-impl.hh>
#include <honei/libutil/shared_array-impl.hh>

#include <iterator>
#include <vector>

namespace honei
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
            class SparseElementIterator;
            /// Our smart implementation of ElementIteratorBase.
            class NonZeroElementIterator;
            /// Our smart implementation of VectorIteratorBase.
            class NonZeroRowIterator;

            typedef typename Matrix<DataType_>::MatrixElementIterator MatrixElementIterator;

        public:
            friend class SparseElementIterator;
            friend class NonZeroElementIterator;
            friend class NonZeroRowIterator;

            /// Type of the const iterator over our elements.
            typedef typename Matrix<DataType_>::ConstElementIterator ConstElementIterator;

            /// Type of the iterator over our elements.
            typedef typename MutableMatrix<DataType_>::ElementIterator ElementIterator;

            /// Type of the const iterator over our vectors.
            typedef VectorIteratorWrapper<DataType_, SparseVector<DataType_> > ConstRowIterator;

            /// Type of the iterator over our vectors.
            typedef VectorIteratorWrapper<DataType_, SparseVector<DataType_> > RowIterator;

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
                _zero_vector(columns, 1)
            {
                CONTEXT("When creating SparseMatrix:");
                ASSERT(rows > 0, "number of rows is zero!");
                ASSERT(columns > 0, "number of columns is zero!");

                _row_vectors[rows].reset(new SparseVector<DataType_>(columns, 1));
            }

            ~SparseMatrix()
            {
            }

            /// \}

            /// Returns iterator pointing to the first element of the matrix.
            virtual ConstElementIterator begin_elements() const
            {
                return ConstElementIterator(new SparseElementIterator(*this, 0));
            }

            /// Returns iterator pointing behind the last element of the matrix.
            virtual ConstElementIterator end_elements() const
            {
                return ConstElementIterator(new SparseElementIterator(*this, _rows * _columns));
            }

            /// Returns iterator pointing to the first element of the matrix.
            virtual ElementIterator begin_elements()
            {
                return ElementIterator(new SparseElementIterator(*this, 0));
            }

            /// Returns iterator pointing behind the last element of the matrix.
            virtual ElementIterator end_elements()
            {
                return ElementIterator(new SparseElementIterator(*this, _rows * _columns));
            }

            /// Returns const iterator pointing to the first non-zero element of the matrix.
            virtual ConstElementIterator begin_non_zero_elements() const
            {
                return ConstElementIterator(new NonZeroElementIterator(*this));
            }

            /// Returns const iterator pointing behind the last element of the matrix.
            virtual ConstElementIterator end_non_zero_elements() const
            {
                return ConstElementIterator(new NonZeroElementIterator(*this, 0 /* Dummy */));
            }

            /// Returns iterator pointing to the first non-zero element of the matrix.
            virtual ElementIterator begin_non_zero_elements()
            {
                return ElementIterator(new NonZeroElementIterator(*this));
            }

            /// Returns iterator pointing behind the last element of the matrix.
            virtual ElementIterator end_non_zero_elements()
            {
                return ElementIterator(new NonZeroElementIterator(*this, 0 /* Dummy */));
            }

            /// Returns const iterator pointing to the first row of the matrix.
            ConstRowIterator begin_non_zero_rows() const
            {
                return ConstRowIterator(new NonZeroRowIterator(*this, 0));
            }

            /// Returns const iterator pointing behind the last row of the matrix.
            ConstRowIterator end_non_zero_rows() const
            {
                return ConstRowIterator(new NonZeroRowIterator(*this, _rows));
            }

            /// Returns iterator pointing to the first row of the matrix.
            RowIterator begin_non_zero_rows()
            {
                return RowIterator(new NonZeroRowIterator(*this, 0));
            }

            /// Returns iterator pointing behind the last row of the matrix.
            RowIterator end_non_zero_rows()
            {
                return RowIterator(new NonZeroRowIterator(*this, _rows));
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

            /// Retrieves element at (row, column), unassignable.
            inline virtual const DataType_ & operator() (unsigned long row, unsigned long column) const
            {
                if (! _row_vectors[row])
                    /// \todo access element directly
                    return _zero_vector[0];

                /// \todo access element directly
                return (*_row_vectors[row])[column];
            }

            /// Retrieves element at (row, column), assignable.
            inline virtual DataType_ & operator() (unsigned long row, unsigned long column)
            {
                if (! _row_vectors[row])
                    _row_vectors[row].reset(new SparseVector<DataType_>(_columns, _capacity));

                /// \todo access element directly
                return (*_row_vectors[row])[column];
            }

            /// Returns a copy of the matrix.
            virtual SparseMatrix copy() const
            {
                SparseMatrix result(_columns, _rows, _capacity);

                for (unsigned long i(0) ; i < _rows ; ++i)
                {
                    if (_row_vectors[i])
                        result._row_vectors[i].reset(new SparseVector<DataType_>(_row_vectors[i]->copy()));
                }

                return result;
            }
    };

    /**
     * \brief SparseMatrix::SparseElementIterator is a simple iterator implementation for sparse matrices.
     *
     * \ingroup grpmatrix
     **/
    template <> template <typename DataType_> class SparseMatrix<DataType_>::SparseElementIterator :
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
            SparseElementIterator(SparseElementIterator const & other) :
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
            virtual SparseElementIterator & operator++ ()
            {
                CONTEXT("When incrementing iterator by one:");

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

            /// In-place-add operator.
            virtual SparseElementIterator & operator+= (const unsigned long step)
            {
                CONTEXT("When incrementing iterator by '" + stringify(step) + "':");

                _index += step;
                if (_matrix._rows * _matrix._columns >= _index)
                    _index = _matrix._rows + _matrix._columns;

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
                CONTEXT("When accessing assignable element at index '" + stringify(_index) + "':");

                return *_iter;
            }

            /// Dereference operator that returns umassignable reference.
            virtual const DataType_ & operator* () const
            {
                CONTEXT("When accessing unassignable element at index '" + stringify(_index) + "':");

                typename Vector<DataType_>::ConstElementIterator const_iter(_iter);

                return *const_iter;
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
    template <> template <typename DataType_> class SparseMatrix<DataType_>::NonZeroElementIterator :
        public MatrixElementIterator
    {
        private:
            /// Our parent matrix.
            const SparseMatrix<DataType_> & _matrix;

            /// Our index.
            unsigned long _index;

            /// Our row-vector's iterator for the current position.
            typename Vector<DataType_>::ElementIterator _iter;

            /// Our row-vector's iterator for the end.
            typename Vector<DataType_>::ElementIterator _end;

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

        public:
            /// \name Constructors and destructor
            /// \{

            /**
             * Constructor, creates a begin iterator.
             *
             * \param matrix The parent matrix that is referenced by the iterator.
             **/
            NonZeroElementIterator(const SparseMatrix<DataType_> & matrix) :
                _matrix(matrix),
                _index(0),
                _column(0),
                _row(0),
                _iter(matrix._row_vectors[matrix._rows]->begin_non_zero_elements()), // Dummy
                _end(matrix._row_vectors[matrix._rows]->end_non_zero_elements()) // Dummy
            {
                _find_next_row();
            }

            /**
             * Constructor, creates an end iterator.
             *
             * \param matrix The parent matrix that is referenced by the iterator.
             */
            NonZeroElementIterator(const SparseMatrix<DataType_> & matrix, const DataType_ &) :
                _matrix(matrix),
                _index(matrix._rows * matrix._columns  ),
                _column(0),
                _row(matrix._rows),
                _iter(matrix._row_vectors[matrix._rows]->begin_non_zero_elements()),
                _end(matrix._row_vectors[matrix._rows]->end_non_zero_elements())
            {
                _find_next_row();
            }

            /// Copy-constructor.
            NonZeroElementIterator(NonZeroElementIterator const & other) :
                _matrix(other._matrix),
                _index(other._index),
                _iter(other._iter),
                _end(other._end),
                _column(other._columns),
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
            virtual NonZeroElementIterator & operator++ ()
            {
                    ++_iter;
                    _column= _iter.index();

                if (_iter == _end)
                {
                    ++_row;
                    _column = 0;
                    _find_next_row();
                }

                _index = _row * _matrix._columns + _column;

                return *this;
            }

            /// In-place-add operator.
            virtual NonZeroElementIterator & operator+= (const unsigned long step)
            {
                /// \todo

                return *this;
            }

            /// Dereference operator that returns assignable reference.
            virtual DataType_ & operator* ()
            {
                CONTEXT("When accessing assignable element at index '" + stringify(_index) + "':");

                return *_iter;
            }

            /// Dereference operator that returns umassignable reference.
            virtual const DataType_ & operator* () const
            {
                CONTEXT("When accessing unassignable element at index '" + stringify(_index) + "':");

                typename Vector<DataType_>::ConstElementIterator const_iter(_iter);

                return *const_iter;
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
                return _column;
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
     * \brief SparseMatrix::NonZeroRowIterator is a smart iterator implementation that iterates over non-zero
     * \brief rows of sparse matrices.
     *
     * \ingroup grpmatrix
     **/
    template <> template <typename DataType_> class SparseMatrix<DataType_>::NonZeroRowIterator :
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

            /**
             * Constructor, creates a begin iterator.
             *
             * \param matrix The parent matrix that is referenced by the iterator.
             **/
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
    };
}

#endif
