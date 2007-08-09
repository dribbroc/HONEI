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

#include <libla/element_iterator.hh>
#include <libla/vector.hh>
#include <libla/matrix_error.hh>
#include <libla/vector_iterator.hh>
#include <libutil/shared_array.hh>

#include <ostream>
#include <cmath>

/**
 * \file
 *
 * Interface declarations for Matrix-based types.
 *
 * \ingroup grpmatrix
 **/
namespace pg512 ///< \todo Namespace name?
{
    /**
     * \brief Matrix is the abstract base class for all matrix-like types used.
     *
     * \ingroup grpmatrix
     **/
    template <typename DataType_> class Matrix
    {
        public:
            template <typename ElementType_> class ElementIteratorBase;
            typedef ElementIteratorBase<DataType_> MatrixElementIterator;

            /// Type of the const iterator over our elements.
            template <typename ElementType_> class ConstElementIteratorWrapper;
            typedef ConstElementIteratorWrapper<DataType_> ConstElementIterator;

            /// Returns iterator pointing to the first element of the matrix.
            virtual ConstElementIterator begin_elements() const = 0;

            /// Returns iterator pointing behind the last element of the matrix.
            virtual ConstElementIterator end_elements() const = 0;

            /// Returns the number of our columns.
            virtual unsigned long columns() const = 0;

            /// Returns the number of our rows.
            virtual unsigned long rows() const = 0;

            /// Returns a copy of the matrix.
            virtual Matrix * copy() const = 0;

            /// Returns true if the matrix is square.
            virtual inline bool square() const
            {
                return rows() == columns();
            }
    };

    /**
     * \brief MutableMatrix is the abstract interface for matrices with non-const iterators.
     *
     * \ingroup grpmatrix
     **/
    template <typename DataType_> class MutableMatrix
    {
        public:
            typedef typename Matrix<DataType_>::MatrixElementIterator MatrixElementIterator;

            /// Type of the iterator over our elements.
            template <typename ElementType_> class ElementIteratorWrapper;
            typedef ElementIteratorWrapper<DataType_> ElementIterator;

            /// Returns iterator pointing to the first element of the matrix.
            virtual ElementIterator begin_elements() = 0;

            /// Returns iterator pointing behind the last element of the matrix.
            virtual ElementIterator end_elements() = 0;
    };


    /**
     * \brief RowAccessMatrix is the abstract base class for all matrix-like types
     * \brief that offer random access to their rows.
     *
     * \ingroup grpmatrix
     **/
    template <typename DataType_> class RowAccessMatrix :
        public Matrix<DataType_>
    {
        public:
            /// Retrieves row by index, zero-based, unassignable.
            virtual const Vector<DataType_> & operator[] (unsigned long row) const = 0;

            /// Retrieves row by index, zero-based, assignable
            virtual Vector<DataType_> & operator[] (unsigned long row) = 0;

            /// \todo Iteration over rows.
    };

    /**
     * \brief Matrix::ElementIteratorBase declares the minimal interface for any ElementIterator implementation
     * \brief for matrix-like types.
     *
     * \ingroup grpmatrix
     */
    template <> template <typename DataType_> class Matrix<DataType_>::ElementIteratorBase :
        public IteratorBase<DataType_, Matrix<DataType_> >
    {
    };

    /**
     * \brief MutableMatrix::ElementIteratorWrapper provides a covariant mutable forward iterator that wraps the actual
     * \brief ElementIteratorBase implementations of any of Matrix's descendants.
     *
     * \ingroup grpmatrix
     */
    template <> template <typename DataType_> class MutableMatrix<DataType_>::ElementIteratorWrapper<DataType_> :
        public std::iterator<std::forward_iterator_tag, DataType_>
    {
        private:
            std::tr1::shared_ptr<MatrixElementIterator> _iterator;

        public:
            friend class Matrix<DataType_>::ConstElementIterator;

            /// \name Constructors
            /// \{

            /**
             * Constructor.
             *
             * \param iterator An instance of one of ElementIteratorBase's descendants that shall be wrapped.
             */
            ElementIteratorWrapper(MatrixElementIterator * iterator) :
                _iterator(iterator)
            {
                if (! iterator)
                    throw std::string("Eek. Iterator is 0, that should not happen!");
            }

            /// Copy-constructor.
            ElementIteratorWrapper(const ElementIteratorWrapper<DataType_> & other) :
                _iterator(other._iterator)
            {
            }

            /// \}

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            virtual ElementIteratorWrapper<DataType_> & operator++ ()
            {
                ++(*_iterator);

                return *this;
            }

            /// Dereference operator that returns an assignable reference.
            virtual DataType_ & operator* ()
            {
                return *(*_iterator);
            }

            /// Comparison operator for equality.
            virtual bool operator== (const ElementIteratorWrapper<DataType_> & other) const
            {
                return (*_iterator == *other._iterator);
            }

            /// Comparison operator for inequality.
            virtual bool operator!= (const ElementIteratorWrapper<DataType_> & other) const
            {
                return (*_iterator != *other._iterator);
            }

            /// \}

            /// \name IteratorTraits interface
            /// \{

            /// Returns our index.
            unsigned long index() const
            {
                return _iterator->index();
            }

            /// Returns our column index.
            unsigned long column() const
            {
                return _iterator->column();
            }

            /// Returns our row index.
            unsigned long row() const
            {
                return _iterator->row();
            }

            /// Returns a pointer to our parent container.
            const Matrix<DataType_> * parent() const
            {
                return _iterator->parent();
            }

            /// \}
    };

    /**
     * \brief ConstElementIteratorWrapper provides a covariant const forward iterator that wraps the actual
     * \brief ElementIteratorBase implementations of any of Vector's descendants.
     *
     * \ingroup grpmatrix
     */
    template <> template <typename DataType_> class Matrix<DataType_>::ConstElementIteratorWrapper<DataType_> :
        public std::iterator<std::forward_iterator_tag, const DataType_>
    {
        private:
            std::tr1::shared_ptr<ElementIteratorBase<DataType_> > _iterator;

        public:
            /// \name Constructors
            /// \{

            /**
             * Constructor.
             *
             * \param iterator An instance of one of ElementIteratorBase's descendants that shall be wrapped.
             */
            ConstElementIteratorWrapper(ElementIteratorBase<DataType_> * iterator) :
                _iterator(iterator)
            {
                if (! iterator)
                    throw std::string("Eek. Iterator is 0, that should not happen!");
            }

            /// Copy-constructor.
            ConstElementIteratorWrapper(const ConstElementIteratorWrapper<DataType_> & other) :
                _iterator(other._iterator)
            {
            }

            /// Conversion-constructor from mutable ElementIteratorWrapper.
            ConstElementIteratorWrapper(const typename MutableMatrix<DataType_>::ElementIterator & other) :
                _iterator(other._iterator)
            {
            }

            /// \}

            /// \name Forward iterator interface
            /// \{

            /// Preincrement operator.
            virtual ConstElementIteratorWrapper<DataType_> & operator++ ()
            {
                ++(*_iterator);

                return *this;
            }

            /// Dereference operator that returns an unassignable reference.
            virtual const DataType_ & operator* () const
            {
                const ElementIteratorBase<DataType_> & iterator(*_iterator);

                return *iterator;
            }

            /// Comparison operator for equality.
            virtual bool operator== (const ConstElementIteratorWrapper<DataType_> & other) const
            {
                return (*_iterator == *other._iterator);
            }

            /// Comparison operator for inequality.
            virtual bool operator!= (const ConstElementIteratorWrapper<DataType_> & other) const
            {
                return (*_iterator != *other._iterator);
            }

            /// \}

            /// \name IteratorTraits interface
            /// \{

            /// Returns our index.
            virtual unsigned long index() const
            {
                return _iterator->index();
            }

            /// Returns our column index.
            unsigned long column() const
            {
                return _iterator->column();
            }

            /// Returns our row index.
            unsigned long row() const
            {
                return _iterator->row();
            }

            /// Returns a pointer to our parent container.
            virtual const Matrix<DataType_> * parent() const
            {
                return _iterator->parent();
            }

            /// \}
    };

    /**
     * \brief Output our Matrix to an ostream.
     *
     * \ingroup grpmatrixoperations
     **/
    template <typename DataType_> std::ostream & operator<< (std::ostream & lhs, const Matrix<DataType_> & m)
    {
        unsigned long row(0);

        lhs << "[ " << std::endl;
        for (typename Matrix<DataType_>::ConstElementIterator i(m.begin_elements()), i_end(m.end_elements()) ;
                i != i_end ; ++i)
        {
            if (row != i.row())
            {
                lhs << std::endl;
                row = i.row();
            }
            lhs << " " << *i;
        }

        lhs << "]" << std::endl << "]";

        return lhs;
    }

    /**
     * \brief Compare two Matrices for equality.
     *
     * \ingroup grpmatrixoperations
     **/
    template <typename DataType_> bool operator== (const Matrix<DataType_> & left, const Matrix<DataType_> & right)
    {
        bool result(true);
        if (left.columns() != right.columns())
        {
            throw MatrixColumnsDoNotMatch(right.columns(), left.columns());
        }

        if (left.rows() != right.rows())
        {
            throw MatrixRowsDoNotMatch(right.rows(), left.rows());
        }

        for (typename Matrix<DataType_>::ConstElementIterator i(left.begin_elements()), i_end(left.end_elements()),
                j(right.begin_elements()) ; i != i_end ; ++i)
        {
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

    /**
     * \brief Output our RowAccessMatrix to an ostream.
     *
     * \ingroup grpmatrixoperations
     **/
    template <typename DataType_> std::ostream & operator<< (std::ostream & lhs, const RowAccessMatrix<DataType_> & m)
    {
        lhs << "[ " << std::endl;
        for (unsigned long r(0) ; r < m.rows() ; ++r) ///< \todo Add row-iteration to RowAccessMatrix.
        {
            const Vector<DataType_> & v(m[r]);

            lhs << " " << v << std::endl;
        }
        lhs << "]";

        return lhs;
    }
}

#endif
