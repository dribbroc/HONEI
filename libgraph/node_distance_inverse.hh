/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Sven Mallach <sven.mallach@uni-dortmund.de>
 * Copyright (c) 2007 Nina Harmuth <nina.harmuth@uni-dortmund.de>
 * Copyright (c) 2007 Thorsten Deinert <thorsten.deinert@uni-dortmund.de>
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

#ifndef LIBGRAPH_GUARD_NODE_DISTANCE_INVERSE_HH
#define LIBGRAPH_GUARD_NODE_DISTANCE_INVERSE_HH 1

#include <libutil/tags.hh>
#include <libla/dense_matrix.hh>
#include <libla/banded_matrix.hh>
#include <libla/dense_vector.hh>
#include <libla/matrix_error.hh>
#include <libla/vector_norm.hh>
#include <libla/vector_difference.hh>

/**
 * \file
 *
 * Implementation of NodeDistanceInverse. </br>
 *
 * \ingroup grplibgraph
 **/
namespace pg512
{
    /**
     * \brief NodeDistanceInverse is used in the algorithm of Fruchterman-Reingold.
     * \brief Hardly surprising, NodeDistanceInverse computes the inverse distance between nodes,
     * \brief which positions are given in a position matrix.
     *
     * \ingroup grplibgraph
     **/
    template <typename Tag_ = tags::CPU>
    struct NodeDistanceInverse
    {
        /**
         * Returns the resulting inverse distance matrix.
         * \param pos_matrix The matrix with the positions of nodes. It may have as many rows as needed for all dimensions
         */
        template <typename DataType_>
        static DenseMatrix<DataType_> value(const DenseMatrix<DataType_> & pos_matrix)
        {
            // Create the result matrix
            DenseMatrix<DataType_> result(pos_matrix.columns(), pos_matrix.columns());

            // Retrieving an ElementIterator e initialized to the first element in the result matrix
            typename MutableMatrix<DataType_>::ElementIterator e(result.begin_elements());

            // Iterate over all nodes (=columns) in the given position matrix
            for (int i =  0; i < pos_matrix.columns(); ++i)
            {
                // The column-vector represents all n coordinates of node i, so we grab them from the matrix
                DenseVector<DataType_> v = pos_matrix.column(i);

                // Now, iterate over all nodes to calculate the distance between node i and node j. For the resulting
                // distance matrix is symmetric, result(i,j) equals result(j,i) and so were waiting for symmetric matrices
                // to gain performance. So long, the trivial loops.
                for (int j = 0; j < pos_matrix.columns(); ++j, ++e)
                {
                    DenseVector<DataType_> w = pos_matrix.column(j);
                    // Now calculate difference tmp = v - w for each pair of nodes and
                    // then, calculate d = tmp1^2 + tmp2^2 + ... + tmpN^2 which is the
                    // l2-norm of tmp without a root. (see template parameter "root")
                    // If d != 0, we put 1/d into the resulting distance matrix, otherwise we write 0.
                    DataType_ d = VectorNorm<DataType_>::value(VectorDifference<>::value(*v.copy(), w));
                    *e = d == 0 ? 0 : 1 / d;
                }
            }

            return result;
        }
    };
}

#endif
