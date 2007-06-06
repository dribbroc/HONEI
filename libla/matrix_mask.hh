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

#ifndef LIBLA_GUARD_MATRIX_MASK_HH
#define LIBLA_GUARD_MATRIX_MASK_HH 1

#include <libla/tags.hh>
#include <libla/dense_matrix.hh>
#include <libla/banded_matrix.hh>
#include <libla/matrix_error.hh>

/**
 * \file
 *
 * Implementation of MatrixMask.
 *
 * \ingroup grpmatrixoperations
 **/
namespace pg512
{
    /**
     * \brief MatrixMask masks a given matrix by setting all elements to zero, for which the mask matrix contains the value "false".
     *
     * \ingroup grpmatrixoperations
     **/
    template <typename DataType_, typename Tag_ = tags::CPU> struct MatrixMask
    {
        static DenseMatrix<DataType_> value(DenseMatrix<DataType_> & matrix, const Matrix<bool> & mask)
        {
            if (matrix.columns() != mask.columns())
            {
                throw MatrixColumnsDoNotMatch(matrix.columns(), mask.columns());
            }

            if (matrix.rows() != mask.rows())
            {
                throw MatrixRowsDoNotMatch(matrix.rows(), mask.rows());
            }

            for (typename MutableMatrix<DataType_>::ElementIterator l(matrix.begin_elements()), l_end(matrix.end_elements()),
                    Matrix<bool>::ConstElementIterator r(mask.begin_elements()), r_end(mask.end_elements()) ; l != l_end ; ++l)
            {
                *l = (*r ? *l : static_cast<DataType_>(0));
                ++r;
            }

            return matrix;
        }

		static BandedMatrix<DataType_> value(BandedMatrix<DataType_> & matrix, const Matrix<bool> & mask)
        {
            if (matrix.columns() != mask.columns())
            {
                throw MatrixColumnsDoNotMatch(matrix.columns(), mask.columns());
            }

            if (matrix.rows() != mask.rows())
            {
                throw MatrixRowsDoNotMatch(matrix.rows(), mask.rows());
            }

            for (typename MutableMatrix<DataType_>::ElementIterator l(matrix.begin_elements()), l_end(matrix.end_elements()),
                    Matrix<bool>::ConstElementIterator r(mask.begin_elements()), r_end(mask.end_elements()) ; l != l_end ; ++l)
            {
                *l = (*r ? *l : static_cast<DataType_>(0));
                ++r;
            }

            return matrix;
        }

        /// \todo overload MaskMatrix::value for sparse matrices.
    };
}
#endif
