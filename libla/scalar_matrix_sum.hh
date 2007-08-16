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

#ifndef LIBLA_GUARD_SCALAR_MATRIX_SUM_HH
#define LIBLA_GUARD_SCALAR_MATRIX_SUM_HH 1

#include <libutil/tags.hh>
#include <libla/dense_matrix.hh>
#include <libla/matrix.hh>
#include <libla/sparse_matrix.hh>

/**
 * \file
 *
 * Templatized definitions of scalar-matrix sums.
 *
 * \ingroup grpmatrixoperations
 **/
namespace pg512
{
    /**
     * ScalarMatrixSum is the class template for adding a scalar to a matrix.
     * \brief With a negative scalar, ScalarMatrixSum can also be used to compute a ScalarMatrixDifference.
     * \ingroup grpmatrixoperations
     **/
    template <typename Tag_ = tags::CPU> struct ScalarMatrixSum
    {
        /**
         * Returns the resulting matrix after add a scalar to a given DenseMatrix instance.
         *
         * \param matrix The DenseMatrix to be used.
         * \param scalar The scalar to be added.
         *
         * \retval matrix The referenced matrix is changed by adding the given scalar to each of its elements.
         **/
        template <typename DT1_, typename DT2_>
        static DenseMatrix<DT2_> & value(const DT1_ scalar, DenseMatrix<DT2_> & matrix)
        {
            for (typename MutableMatrix<DT2_>::ElementIterator l(matrix.begin_elements()),
                    l_end(matrix.end_elements()) ; l != l_end ; ++l)
            {
                *l += scalar;
            }

            return matrix;
        }
    };
}

#endif
