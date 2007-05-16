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
#include <libutil/shared_array.hh>

#include <ostream>
#include <string.h>

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

            /// Returns iterator pointing to the first element of the matrix.
            virtual ElementIterator begin_elements() const = 0;

            /// Returns iterator pointing behind the last element of the matrix.
            virtual ElementIterator end_elements() const = 0;

            /// Returns our number of columns.
            virtual unsigned long columns() const = 0;

            /// Returns our number of rows.
            virtual unsigned long rows() const = 0;

    };

    /// Output our Matrix to an ostream.
    template <typename DataType_> std::ostream & operator<< (std::ostream & lhs, const Matrix<DataType_> & m)
    {
        unsigned long row(0);

        lhs << "[ " << std::endl;
        for (typename Matrix<DataType_>::ElementIterator i(m.begin_elements()), i_end(m.end_elements()) ;
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
     * A RowAccessMatrix is the abstract baseclass for all matrix-like types
     * that offere randome access to their rows.
     **/
    template <typename DataType_> class RowAccessMatrix :
        public Matrix<DataType_>
    {
        public:
            /// Retrieves element by index, zero-based, unassignable.
            virtual const Vector<DataType_> & operator[] (unsigned long row) const = 0;

            /// Retrieves element by index, zero-based, assignable
            virtual Vector<DataType_> & operator[] (unsigned long row) = 0;

            /// \todo Iteration over rows.
    };

    /// Output our RowAccessMatrix to an ostream.
    template <typename DataType_> std::ostream & operator<< (std::ostream & lhs, const RowAccessMatrix<DataType_> & m)
    {
        lhs << "[ " << std::endl;
        for (unsigned long r(0) ; r < m.rows() ; ++r) ///< \todo Add row-iteration to RowAccessMatrix.
        {
            const Vector<DataType_> & v(m[r]);

            lhs << " "<< v << std::endl;
        }
        lhs << "]";

        return lhs;
    }
}

#endif
