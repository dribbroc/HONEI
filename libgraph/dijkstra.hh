/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Mathias Kadolsky <mathkad@gmx.de>
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

#ifndef LIBGRAPH_GUARD_DIJKSTRA_HH
#define LIBGRAPH_GUARD_DIJKSTRA_HH 1

#include <libutil/tags.hh>
#include <libla/dense_matrix.hh>
#include <libla/banded_matrix.hh>
#include <libla/dense_vector.hh>
#include <libla/matrix_error.hh>
#include <libla/vector_masked_min_index.hh>

/**
 * \file
 *
 * Implementation of Dijkstra's algorithm.
 *
 * \ingroup grplibgraph
 **/
namespace pg512
{
    /**
     * \brief Dijkstra computes the graph distances between nodes.
     *
     * \ingroup grplibgraph
     **/
    template <typename Tag_ = tags::CPU> struct Dijkstra
    {
        /**
         * Returns the resulting graph distance matrix.
         * \param cost_matrix The matrix with the costs of edges.
         */
        template <typename DataType_>
        static DenseMatrix<DataType_> value(const DenseMatrix<DataType_> & cost_matrix)
        {
           if (cost_matrix.rows() != cost_matrix.columns())
               throw MatrixRowsDoNotMatch(cost_matrix.columns(), cost_matrix.rows());

            // Create the result matrix
            DenseMatrix<DataType_> result(cost_matrix.columns(), cost_matrix.columns(), 0);

            // Create a boolean vector for the active set of nodes in Dijkstra
            /// \todo Please use a real set here, i.e. std::set!
            DenseVector<bool> set_of_nodes(cost_matrix.column(0).size(), false);

            // Index of the current node
            int Index;

            // Auxiliary variables
            DataType_ help1;

            // Iterate over all columns in the given cost matrix
            for (int i =  0; i < cost_matrix.columns(); ++i)
            {
                // The column-vector represents the graph distances between node i and the other nodes, so we grab them from the matrix
                DenseVector<DataType_> v(cost_matrix.column(i));

                // Set the boolean vector representing the active set of nodes
                for (int j = 0; j < v.size(); ++j)
                {
                    i == j ? set_of_nodes[j] = false : set_of_nodes[j] = true; /// \todo Check if 1 means true :-)
                }

                v[i] = 0;

                for (int j = 0; j < v.size() - 1; ++j)
                {
                    // Calculate the node with the shortest distance to the active set of nodes
                    Index = VectorMaskedMinIndex<>::value(v, set_of_nodes);

                    // Delete the current node from the active set of nodes
                    set_of_nodes[Index] = 0;

                    // Update the graph distances from node i to each other node by using the costs of the current node
                    for (int n = 0; n < v.size(); ++n)
                    {
                        help1 = v[Index] + cost_matrix[n][Index];
                        if (set_of_nodes[n] == 1 && v[n] > help1)
                        {
                            v[n] = help1;
                        }
                    }
                }

                 // Write the graph distances from node i to each other node represented by vector v into the corresponding column of the result matrix
                 for (int j = 0; j < v.size() ; ++j)
                 {
                     result[j][i] = v[j];
                 }
            }

            return result;
        }
    };
}
#endif
