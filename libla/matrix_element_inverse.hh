/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#ifndef LIBLA_GUARD_MATRIX_ELEMENT_INVERSE_HH
#define LIBLA_GUARD_MATRIX_ELEMENT_INVERSE_HH 1

#include <libutil/tags.hh>
#include <libla/dense_matrix.hh>
#include <libla/banded_matrix.hh>

/**
 * \file
 *
 * Templatized definitions of inverting a matrix's elements.<br/>
 *
 * \ingroup grpmatrixoperations
 **/
namespace pg512
{
    /**
     * VectorAbsolute is the class template that inverts a matrix's elements.
     *
     * \ingroup grpmatrixoperations
     **/
    template <typename Tag_ = tags::CPU> struct MatrixElementInverse
    {
        /**
         * Return a matrix's inverse elements. All elements that equal zero will
         * be invariant under this operation.
         *
         * \param matrix DenseMatrix whose non-zero elements shall be inverted.
         **/
        template <typename DataType_> static DenseMatrix<DataType_> & value(DenseMatrix<DataType_> & matrix)
        {
            for (typename MutableMatrix<DataType_>::ElementIterator i(matrix.begin_elements()), i_end(matrix.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (*i == static_cast<DataType_>(0))
                    continue;

                *i = DataType_(1) / *i;
            }

            return matrix;
        }

        template <typename DataType_> static BandedMatrix<DataType_> & value(BandedMatrix<DataType_> & matrix)
        {
            for (typename MutableMatrix<DataType_>::ElementIterator i(matrix.begin_elements()), i_end(matrix.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (*i == static_cast<DataType_>(0))
                    continue;

                *i = DataType_(1) / *i;
            }

            return matrix;
        }

    };
}

#endif
