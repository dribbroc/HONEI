/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Michael Abshoff <michael.abshoff@fsmath.mathematik.uni-dortmund.de>
 * Copyright (c) 2007 Sven Mallach <sven.mallach@uni-dortmund.de>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef HONEI_GUARD_LA_DENSE_MATRIX_IMPL_HH
#define HONEI_GUARD_LA_DENSE_MATRIX_IMPL_HH 1

#include <honei/la/element_iterator.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/dense_vector-impl.hh>
#include <honei/la/dense_vector_range-impl.hh>
#include <honei/la/dense_vector_slice-impl.hh>
#include <honei/la/matrix_error.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/util/private_implementation_pattern-impl.hh>
#include <honei/util/shared_array-impl.hh>
#include <honei/util/stringify.hh>
#include <honei/util/type_traits.hh>

namespace honei
{
    template <typename DataType_> struct Implementation<DenseMatrix<DataType_> >
    {
        /// Pointer to our elements.
        SharedArray<DataType_> elements;

        /// Our columns.
        unsigned long columns;

        /// Our rows.
        unsigned long rows;

        /// Our size.
        unsigned long size;

        /// Constructor
        Implementation(unsigned long rows, unsigned long columns) :
            elements(rows * columns),
            columns(columns),
            rows(rows),
            size(rows * columns)
        {
        }

        DataType_ & operator() (unsigned long row, unsigned long column)
        {
            CONTEXT("When retrieving DenseMatrix element (internal):");
            ASSERT(row < rows, "row index is out of bounds!");
            ASSERT(column < columns, "column number is out of bounds!");

            return elements[columns * row + column];
        }
    };

    template <typename DataType_>
    DenseMatrix<DataType_>::DenseMatrix(unsigned long rows, unsigned long columns) :
        PrivateImplementationPattern<DenseMatrix<DataType_>, Shared>(new Implementation<DenseMatrix<DataType_> >(rows, columns))
    {
        CONTEXT("When creating DenseMatrix:");
        ASSERT(rows > 0, "number of rows is zero!");
        ASSERT(columns > 0, "number of columns is zero!");
    }

    template <typename DataType_>
    DenseMatrix<DataType_>::DenseMatrix(unsigned long rows, unsigned long columns, DataType_ value) :
        PrivateImplementationPattern<DenseMatrix<DataType_>, Shared>(new Implementation<DenseMatrix<DataType_> >(rows, columns))
    {
        CONTEXT("When creating DenseMatrix:");
        ASSERT(rows > 0, "number of rows is zero!");
        ASSERT(columns > 0, "number of columns is zero!");

        TypeTraits<DataType_>::fill(this->_imp->elements.get(), this->_imp->size, value);
    }

    template <typename DataType_>
    DenseMatrix<DataType_>::DenseMatrix(const SparseMatrix<DataType_> & other) :
        PrivateImplementationPattern<DenseMatrix<DataType_>, Shared>(new Implementation<DenseMatrix<DataType_> >(other.rows(), other.columns()))
    {
        CONTEXT("When creating DenseMatrix form SparseMatrix:");

        TypeTraits<DataType_>::fill(this->_imp->elements.get(), this->_imp->size, DataType_(0));

        for (typename SparseMatrix<DataType_>::NonZeroConstElementIterator i(other.begin_non_zero_elements()),
                i_end(other.end_non_zero_elements()) ; i != i_end ; ++i)
        {
            (*this->_imp)(i.row(), i.column()) = *i;
        }
    }

    template <typename DataType_>
    DenseMatrix<DataType_>::DenseMatrix(const DenseMatrix<DataType_> & source, unsigned long column_offset, unsigned long columns,
            unsigned long row_offset, unsigned long rows) :
        PrivateImplementationPattern<DenseMatrix<DataType_>, Shared>(new Implementation<DenseMatrix<DataType_> >(rows, columns))
    {
        ASSERT(rows > 0, "number of rows is zero!");
        ASSERT(columns > 0, "number of columns is zero!");

        if (column_offset + columns > source.columns())
        {
            throw MatrixColumnsDoNotMatch(column_offset + columns, source.columns());
        }

        if (row_offset + rows > source.rows())
        {
            throw MatrixRowsDoNotMatch(row_offset + rows, source.rows());
        }

        for (unsigned long i(0) ; i < rows ; ++i)
        {
            for (unsigned long j(0); j < columns ; ++j)
            {
                this->_imp->elements[j + columns * i] = source._imp->elements[j + column_offset  +
                    ((i + row_offset) * source.columns())];
            }
        }
    }

    template <typename DataType_>
    DenseMatrix<DataType_>::~DenseMatrix()
    {
    }

    template <typename DataType_>
    typename DenseMatrix<DataType_>::ConstElementIterator
    DenseMatrix<DataType_>::begin_elements() const
    {
        return typename DenseMatrix<DataType_>::ConstElementIterator(*this, 0);
    }

    template <typename DataType_>
    typename DenseMatrix<DataType_>::ConstElementIterator
    DenseMatrix<DataType_>::element_at(unsigned long index) const
    {
         return typename DenseMatrix<DataType_>::ConstElementIterator(*this, index);
    }

    template <typename DataType_>
    typename DenseMatrix<DataType_>::ConstElementIterator
    DenseMatrix<DataType_>::end_elements() const
    {
        return typename DenseMatrix<DataType_>::ConstElementIterator(*this, this->_imp->size);
    }

    template <typename DataType_>
    typename DenseMatrix<DataType_>::ElementIterator
    DenseMatrix<DataType_>::begin_elements()
    {
        return typename DenseMatrix<DataType_>::ElementIterator(*this, 0);
    }

    template <typename DataType_>
    typename DenseMatrix<DataType_>::ElementIterator
    DenseMatrix<DataType_>::element_at(unsigned long index)
    {
         return typename DenseMatrix<DataType_>::ElementIterator(*this, index);
    }

    template <typename DataType_>
    typename DenseMatrix<DataType_>::ElementIterator
    DenseMatrix<DataType_>::end_elements()
    {
        return typename DenseMatrix<DataType_>::ElementIterator(*this, this->_imp->size);
    }

    template <typename DataType_>
    unsigned long
    DenseMatrix<DataType_>::columns() const
    {
        return this->_imp->columns;
    }

    template <typename DataType_>
    unsigned long
    DenseMatrix<DataType_>::rows() const
    {
        return this->_imp->rows;
    }

    template <typename DataType_>
    unsigned long
    DenseMatrix<DataType_>::size() const
    {
        return this->_imp->size;
    }

    template <typename DataType_>
    typename DenseMatrix<DataType_>::ConstRow
    DenseMatrix<DataType_>::operator[] (unsigned long row) const
    {
        return DenseVectorRange<DataType_>(this->_imp->elements, this->_imp->columns, row * this->_imp->columns);
    }

    template <typename DataType_>
    typename DenseMatrix<DataType_>::Row
    DenseMatrix<DataType_>::operator[] (unsigned long row)
    {
        return DenseVectorRange<DataType_>(this->_imp->elements, this->_imp->columns, row * this->_imp->columns);
    }

    template <typename DataType_>
    const DataType_ &
    DenseMatrix<DataType_>::operator() (unsigned long row, unsigned long column) const
    {
        return (*this->_imp)(row, column);
    }

    template <typename DataType_>
    DataType_ &
    DenseMatrix<DataType_>::operator() (unsigned long row, unsigned long column)
    {
        return (*this->_imp)(row, column);
    }

    template <typename DataType_>
    const DataType_ &
    DenseMatrix<DataType_>::operator() (unsigned long index) const
    {
        CONTEXT("When retrieving DenseMatrix element, unassignable:");
        ASSERT(index < this->_imp->size, "index is out of bounds!");

        return this->_imp->elements[index];
    }

    template <typename DataType_>
    DataType_ &
    DenseMatrix<DataType_>::operator() (unsigned long index)
    {
        CONTEXT("When retrieving DenseMatrix element, assignable:");
        ASSERT(index < this->_imp->size, "index is out of bounds!");

        return this->_imp->elements[index];
    }

    template <typename DataType_>
    typename DenseMatrix<DataType_>::ConstColumn
    DenseMatrix<DataType_>::column(unsigned long column) const
    {
        return DenseVectorSlice<DataType_>(this->_imp->elements, this->_imp->rows, column, this->_imp->columns);
    }

    template <typename DataType_>
    typename DenseMatrix<DataType_>::Column
    DenseMatrix<DataType_>::column(unsigned long column)
    {
        return DenseVectorSlice<DataType_>(this->_imp->elements, this->_imp->rows, column, this->_imp->columns);
    }

    template <typename DataType_>
    DataType_ *
    DenseMatrix<DataType_>::elements() const
    {
        return this->_imp->elements.get();
    }

    template <typename DataType_>
    inline
    void *
    DenseMatrix<DataType_>::address() const
    {
        return this->_imp->elements.get();
    }

    template <typename DataType_>
    inline
    void *
    DenseMatrix<DataType_>::memid() const
    {
        return this->_imp->elements.get();
    }

    template <typename DataType_>
    void *
    DenseMatrix<DataType_>::lock(LockMode mode, tags::TagValue memory) const
    {
        return MemoryArbiter::instance()->lock(mode, memory, this->memid(), this->address(), this->size() * sizeof(DataType_));
    }

    template <typename DataType_>
    void
    DenseMatrix<DataType_>::unlock(LockMode mode) const
    {
        MemoryArbiter::instance()->unlock(mode, this->memid());
    }

    template <typename DataType_>
    DenseMatrix<DataType_>
    DenseMatrix<DataType_>::copy() const
    {
        DenseMatrix result(this->_imp->rows, this->_imp->columns);
        this->lock(lm_read_only);

        TypeTraits<DataType_>::copy(this->_imp->elements.get(), result._imp->elements.get(), this->_imp->size);

        this->unlock(lm_read_only);
        return result;
    }

    template <typename DataType_> struct Implementation<ConstElementIterator<storage::Dense, container::Matrix, DataType_> >
    {
        SharedArray<DataType_> elements;

        unsigned long columns;

        unsigned long index;

        Implementation(const SharedArray<DataType_> & elements, unsigned long columns, unsigned long index) :
            elements(elements),
            columns(columns),
            index(index)
        {
        }

        Implementation(const Implementation<ElementIterator<storage::Dense, container::Matrix, DataType_> > & other) :
            elements(other.elements),
            columns(other.columns),
            index(other.index)
        {
        }
    };

    template <typename DataType_>
    ConstElementIterator<storage::Dense, container::Matrix, DataType_>::ConstElementIterator(const DenseMatrix<DataType_> & matrix, unsigned long index) :
        PrivateImplementationPattern<ConstElementIterator<storage::Dense, container::Matrix, DataType_>, Single>(
                new Implementation<ConstElementIterator<storage::Dense, container::Matrix, DataType_> >(matrix._imp->elements, matrix._imp->columns, index))
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::Dense, container::Matrix, DataType_>::ConstElementIterator(const ConstElementIterator & other) :
        PrivateImplementationPattern<ConstElementIterator<storage::Dense, container::Matrix, DataType_>, Single>(
                new Implementation<ConstElementIterator<storage::Dense, container::Matrix, DataType_> >(*other._imp))
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::Dense, container::Matrix, DataType_>::ConstElementIterator(
            const ElementIterator<storage::Dense, container::Matrix, DataType_> & other) :
        PrivateImplementationPattern<ConstElementIterator<storage::Dense, container::Matrix, DataType_>, Single>(
                new Implementation<ConstElementIterator<storage::Dense, container::Matrix, DataType_> >(*other._imp))
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::Dense, container::Matrix, DataType_>::~ConstElementIterator()
    {
    }

    template <typename DataType_>
    ConstElementIterator<storage::Dense, container::Matrix, DataType_> &
    ConstElementIterator<storage::Dense, container::Matrix, DataType_>::operator= (
            const ConstElementIterator<storage::Dense, container::Matrix, DataType_> & other)
    {
        if (&other == this)
            return *this;

        this->_imp->elements = other._imp->elements;
        this->_imp->columns = other._imp->columns;
        this->_imp->index = other._imp->index;

        return *this;
    }

    template <typename DataType_>
    ConstElementIterator<storage::Dense, container::Matrix, DataType_> &
    ConstElementIterator<storage::Dense, container::Matrix, DataType_>::operator++ ()
    {
        CONTEXT("When incrementing ConstElementIterator<Dense, Matrix> by one:");

        ++this->_imp->index;

        return *this;
    }

    template <typename DataType_>
    ConstElementIterator<storage::Dense, container::Matrix, DataType_> &
    ConstElementIterator<storage::Dense, container::Matrix, DataType_>::operator+= (const unsigned long step)
    {
        CONTEXT("When incrementing ConstElementIterator<Dense, Matrix> by '" + stringify(step) + "':");

        this->_imp->index += step;

        return *this;
    }

    template <typename DataType_>
    const DataType_ &
    ConstElementIterator<storage::Dense, container::Matrix, DataType_>::operator* () const
    {
        CONTEXT("When accessing unassignable element at index '" + stringify(this->_imp->index) + "':");

        return this->_imp->elements[this->_imp->index];
    }

    template <typename DataType_>
    bool
    ConstElementIterator<storage::Dense, container::Matrix, DataType_>::operator< (
            const ConstElementIterator<storage::Dense, container::Matrix, DataType_> & other) const
    {
        return this->_imp->index < other._imp->index;
    }

    template <typename DataType_>
    bool
    ConstElementIterator<storage::Dense, container::Matrix, DataType_>::operator== (
            const ConstElementIterator<storage::Dense, container::Matrix, DataType_> & other) const
    {
        return ((this->_imp->elements.get() == other._imp->elements.get()) && (this->_imp->index == other._imp->index));
    }

    template <typename DataType_>
    bool
    ConstElementIterator<storage::Dense, container::Matrix, DataType_>::operator!= (
            const ConstElementIterator<storage::Dense, container::Matrix, DataType_> & other) const
    {
        return ((this->_imp->elements.get() != other._imp->elements.get()) || (this->_imp->index != other._imp->index));
    }

    template <typename DataType_>
    unsigned long
    ConstElementIterator<storage::Dense, container::Matrix, DataType_>::column() const
    {
        return this->_imp->index % this->_imp->columns;
    }

    template <typename DataType_>
    unsigned long
    ConstElementIterator<storage::Dense, container::Matrix, DataType_>::index() const
    {
        return this->_imp->index;
    }

    template <typename DataType_>
    unsigned long
    ConstElementIterator<storage::Dense, container::Matrix, DataType_>::row() const
    {
        return this->_imp->index / this->_imp->columns;
    }

    template <typename DataType_> struct Implementation<ElementIterator<storage::Dense, container::Matrix, DataType_> >
    {
        SharedArray<DataType_> elements;

        unsigned long columns;

        unsigned long index;

        Implementation(const SharedArray<DataType_> & elements, unsigned long columns, unsigned long index) :
            elements(elements),
            columns(columns),
            index(index)
        {
        }
    };

    template <typename DataType_>
    ElementIterator<storage::Dense, container::Matrix, DataType_>::ElementIterator(const DenseMatrix<DataType_> & matrix, unsigned long index) :
        PrivateImplementationPattern<ElementIterator<storage::Dense, container::Matrix, DataType_>, Single>(
                new Implementation<ElementIterator<storage::Dense, container::Matrix, DataType_> >(matrix._imp->elements, matrix._imp->columns, index))
    {
    }

    template <typename DataType_>
    ElementIterator<storage::Dense, container::Matrix, DataType_>::ElementIterator(const ElementIterator & other) :
        PrivateImplementationPattern<ElementIterator<storage::Dense, container::Matrix, DataType_>, Single>(
                new Implementation<ElementIterator<storage::Dense, container::Matrix, DataType_> >(*other._imp))
    {
    }

    template <typename DataType_>
    ElementIterator<storage::Dense, container::Matrix, DataType_>::~ElementIterator()
    {
    }

    template <typename DataType_>
    ElementIterator<storage::Dense, container::Matrix, DataType_> &
    ElementIterator<storage::Dense, container::Matrix, DataType_>::operator= (
            const ElementIterator<storage::Dense, container::Matrix, DataType_> & other)
    {
        if (&other == this)
            return *this;

        this->_imp->elements = other._imp->elements;
        this->_imp->columns = other._imp->columns;
        this->_imp->index = other._imp->index;

        return *this;
    }

    template <typename DataType_>
    ElementIterator<storage::Dense, container::Matrix, DataType_> &
    ElementIterator<storage::Dense, container::Matrix, DataType_>::operator++ ()
    {
        CONTEXT("When incrementing ElementIterator<Dense, Matrix> by one:");

        ++this->_imp->index;

        return *this;
    }

    template <typename DataType_>
    ElementIterator<storage::Dense, container::Matrix, DataType_> &
    ElementIterator<storage::Dense, container::Matrix, DataType_>::operator+= (const unsigned long step)
    {
        CONTEXT("When incrementing ElementIterator<Dense, Matrix> by '" + stringify(step) + "':");

        this->_imp->index += step;

        return *this;
    }

    template <typename DataType_>
    DataType_ &
    ElementIterator<storage::Dense, container::Matrix, DataType_>::operator* () const
    {
        CONTEXT("When accessing assignable element at index '" + stringify(this->_imp->index) + "':");

        return this->_imp->elements[this->_imp->index];
    }

    template <typename DataType_>
    bool
    ElementIterator<storage::Dense, container::Matrix, DataType_>::operator< (
            const ElementIterator<storage::Dense, container::Matrix, DataType_> & other) const
    {
        return this->_imp->index < other._imp->index;
    }

    template <typename DataType_>
    bool
    ElementIterator<storage::Dense, container::Matrix, DataType_>::operator== (
            const ElementIterator<storage::Dense, container::Matrix, DataType_> & other) const
    {
        return ((this->_imp->elements.get() == other._imp->elements.get()) && (this->_imp->index == other._imp->index));
    }

    template <typename DataType_>
    bool
    ElementIterator<storage::Dense, container::Matrix, DataType_>::operator!= (
            const ElementIterator<storage::Dense, container::Matrix, DataType_> & other) const
    {
        return ((this->_imp->elements.get() != other._imp->elements.get()) || (this->_imp->index != other._imp->index));
    }

    template <typename DataType_>
    unsigned long
    ElementIterator<storage::Dense, container::Matrix, DataType_>::column() const
    {
        return this->_imp->index % this->_imp->columns;
    }

    template <typename DataType_>
    unsigned long
    ElementIterator<storage::Dense, container::Matrix, DataType_>::index() const
    {
        return this->_imp->index;
    }

    template <typename DataType_>
    unsigned long
    ElementIterator<storage::Dense, container::Matrix, DataType_>::row() const
    {
        return this->_imp->index / this->_imp->columns;
    }

    template <typename DataType_>
    bool
    operator== (const DenseMatrix<DataType_> & left, const DenseMatrix<DataType_> & right)
    {
        CONTEXT("When comparing two dense matrices:");

        bool result(true);

        if (left.columns() != right.columns())
        {
            throw MatrixColumnsDoNotMatch(right.columns(), left.columns());
        }

        if (left.rows() != right.rows())
        {
            throw MatrixRowsDoNotMatch(right.rows(), left.rows());
        }

        for (typename DenseMatrix<DataType_>::ConstElementIterator i(left.begin_elements()), i_end(left.end_elements()),
                j(right.begin_elements()) ; i != i_end ; ++i)
        {
            CONTEXT("When comparing elements at index '" + stringify(i.index()) + "':");

            if (fabs((*i - *j)) <= std::numeric_limits<DataType_>::epsilon())
            {
                ++j;
                continue;
            }

            result = false;
            break;
        }

        return result;
    }

    template <typename DataType_>
    std::ostream &
    operator<< (std::ostream & lhs, const DenseMatrix<DataType_> & m)
    {
        lhs << "[" << std::endl;
        for (unsigned long r(0) ; r < m.rows() ; ++r) ///< \todo Add row-iteration to DenseMatrix.
        {
            lhs << " " << m[r];
        }
        lhs << "]" << std::endl;

        return lhs;
    }
}

#endif
