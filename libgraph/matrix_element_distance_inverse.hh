/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Sven Mallach <sven.mallach@uni-dortmund.de>
 * Copyright (c) 2007 Nina Harmuth <nina.harmuth@uni-dortmund.de>
 *
 * This file is part of the Graph C++ library. LibGraph is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibGraph is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef LIBGRAPH_GUARD_MATRIX_ELEMENT_DISTANCE_INVERSE_HH
#define LIBGRAPH_GUARD_MATRIX_ELEMENT_DISTANCE_INVERSE_HH 1

#include <libla/tags.hh>
#include <libla/dense_matrix.hh>
#include <libla/banded_matrix.hh>
#include <libla/dense_vector.hh>
#include <libla/matrix_error.hh>

/**
 * \file
 *
 * Implementation of MatrixElementDistanceInverse. </br>
 *
 * \ingroup grplibgraph
 **/
namespace pg512
{
    /**
     * \brief MatrixElementDistanceInverse is used in the algorithm of Fruchterman-Reingold.
     * \brief MatrixElementDistanceInverse computes the inverse distance between nodes,
     * \brief which positions are given in a position matrix.
     *
     * \ingroup grplibgraph
     **/
    template <typename DataType_, typename Tag_ = tags::CPU> struct MatrixElementDistanceInverse
    {
        /**
         * Returns the resulting inverse distance matrix.
         * \param pos_matrix The matrix with the positions of nodes. It may only have two rows for x- and y-coordinates.
         **/
        static DenseMatrix<DataType_> value(const RowAccessMatrix<DataType_> & pos_matrix)
        {

            if (pos_matrix.rows() != 2)
                throw MatrixRowsDoNotMatch(2, pos_matrix.rows());

	    DenseMatrix<DataType_> result(pos_matrix.columns(), pos_matrix.columns());
	    typename MutableMatrix<DataType_>::ElementIterator e(result.begin_elements());

            for (typename Vector<DataType_>::ConstElementIterator i(pos_matrix[0].begin_elements()),
                    i_end(pos_matrix[0].end_elements()), k(pos_matrix[1].begin_elements()),
                    k_end(pos_matrix[1].end_elements()) ; i != i_end ; ++i, ++k)
            {
                for (typename Vector<DataType_>::ConstElementIterator j(pos_matrix[0].begin_elements()),
                    j_end(pos_matrix[0].end_elements()), l(pos_matrix[1].begin_elements()),
                    l_end(pos_matrix[1].end_elements()) ; j != j_end ; ++j, ++l)
                {  
                    if (*i != *j)
                    {
                        *e = (1 / ((*i -*j) * (*i -*j) + (*k -*l) * (*k-*l)));
                    }
                    else 
                        *e = 0;
                    ++e;
                }
            }

            return result;
        }

    };
}
#endif
