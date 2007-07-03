/* vim: set sw=4 sts=4 et nofoldenable : */

/*
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

#ifndef LIBLA_GUARD_MATRIX_ELEMENT_SUM_HH
#define LIBLA_GUARD_MATRIX_ELEMENT_SUM_HH 1

#include <libutil/tags.hh>
#include <libla/dense_matrix.hh>
#include <libla/banded_matrix.hh>
#include <libla/matrix_error.hh>

/**
 * \file
 *
 * Templatized definitions of matrix element sums.<br/>
 *
 * \ingroup grpmatrixoperations
 **/
namespace pg512
{

    /**
     * MatrixElementSum is the class template for the sum of all elements of a matrix.
	 * The used matrix will be invariant under this operation.
     * \ingroup grpmatrixoperations
     **/
    template <typename Tag_ = tags::CPU> struct MatrixElementSum
    {
        /**
         * Returns the the resulting scalar of the sum of all elements of a given dense matrix instance.
         *
         * \param left Reference to dense matrix which elements will be accumulated.
         **/
        template <typename DataType_> static DataType_ value(const RowAccessMatrix<DataType_> & matrix)
        {
            DataType_ result(0);
			for (typename Matrix<DataType_>::ConstElementIterator l(matrix.begin_elements()),
                    l_end(matrix.end_elements()) ; l != l_end ; ++l)
            {
                result += *l;
            }

            return result;
        }

        /**
         * Returns the the resulting scalar of the sum of all elements of a given banded matrix instance.
         *
         * \param left Reference to banded matrix which elements will be accumulated.
         **/
        template <typename DataType_> static DataType_ value(const BandedMatrix<DataType_> & matrix)
        {
            DataType_ result(0);
			for (typename Matrix<DataType_>::ConstElementIterator l(matrix.begin_elements()),
                    l_end(matrix.end_elements()) ; l != l_end ; ++l)
            {
                result += *l;
            }

            return result;
        }

    };
}
#endif
