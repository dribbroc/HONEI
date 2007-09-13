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
#include <libla/dense_vector-impl.hh>
#include <libla/matrix_error.hh>
#include <libgraph/vector_masked_min_index.hh>

/**
 * \file
 *
 * Implementation of Dijkstra's algorithm.
 *
 * \ingroup grplibgraph
**/
namespace honei
{
    /**
     * \brief Dijkstra computes the graph distances between nodes (and the previous nodes of the shortest path).
     *
     * \ingroup grplibgraph
     **/

    /// \todo improving performance of Dijkstra, e.g. using priority queue

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
            DenseVector<bool> set_of_nodes(cost_matrix.column(0).size(), false);

            // Index of the current node
            int Index;

            // Auxiliary variables
            DataType_ help1;

            // Iterate over all columns in the given cost matrix
            for (int i(0); i < cost_matrix.columns(); ++i)
            {
                // The column-vector represents the graph distances between node i and the other nodes, so we grab them from the matrix
                DenseVector<DataType_> v(cost_matrix.column(i));

                // Set the boolean vector representing the active set of nodes
                for (int j(0); j < v.size(); ++j)
                {
                    i == j ? set_of_nodes[j] = false : set_of_nodes[j] = true;
                }

                v[i] = 0;

                for (int j(0); j < v.size() - 1; ++j)
                {
                    // Calculate the node with the shortest distance to the active set of nodes
                    Index = VectorMaskedMinIndex<>::value(v, set_of_nodes);

                    // Delete the current node from the active set of nodes
                    set_of_nodes[Index] = false;

                    // Update the graph distances from node i to each other node by using the costs of the current node
                    for (int n(0); n < v.size(); ++n)
                    {
                        help1 = v[Index] + cost_matrix[n][Index];
                        if (set_of_nodes[n] == 1 && v[n] > help1)
                        {
                            v[n] = help1;
                        }
                    }
                }

                // Write the graph distances from node i to each other node represented by vector v into the corresponding column of the result matrix
                for (int j(0); j < v.size() ; ++j)
                {
                    result[j][i] = v[j];
                }
            }
            return result;
        }

        /**
         * Computes the resulting graph distance matrix and the previuos nodes matrix.
         * \param cost_matrix The matrix with the costs of edges.
         * \param previous_nodes An empty matrix with the size of cost_matrix.
         */
        template <typename DataType_>
        static void value(DenseMatrix<DataType_> & cost_matrix, DenseMatrix<int> & previous_nodes)
        {
            if (cost_matrix.rows() != cost_matrix.columns())
                throw MatrixRowsDoNotMatch(cost_matrix.columns(), cost_matrix.rows());

            if (cost_matrix.columns() != previous_nodes.columns())
                throw MatrixColumnsDoNotMatch(previous_nodes.columns(), cost_matrix.columns());

            if (cost_matrix.rows() != previous_nodes.rows())
                throw MatrixRowsDoNotMatch(previous_nodes.rows(), cost_matrix.rows());

            // Auxiliary variable
            DataType_ help1;

            // Create the graph_distance_matrix and initialize previous_nodes
            DenseMatrix<DataType_>  graph_distance_matrix(*cost_matrix.copy());            

            for (typename MutableMatrix<int>::ElementIterator e(previous_nodes.begin_elements()),
                    e_end(previous_nodes.end_elements()); e != e_end ; ++e)
                    {
                        *e = e.column();
                    }

            // Create a boolean vector for the active set of nodes in Dijkstra
            DenseVector<bool> set_of_nodes(cost_matrix.column(0).size(), false);

            // Index of the current node
            int Index;

            // Iterate over all columns in the given cost matrix
            for (int i(0); i < cost_matrix.columns(); ++i)
            {
                // The column-vector represents the graph distances between node i and the other nodes, so we grab them from the matrix
                DenseVector<DataType_> v(graph_distance_matrix.column(i));

                // The column-vector represents the previous nodes. If i, p1, p2, j is the shortest path from node i to node j then p(j) = p2.
                DenseVector<int> p(previous_nodes.column(i));

                // Set the boolean vector representing the active set of nodes
                for (int j(0); j < v.size(); ++j)
                {
                    i == j ? set_of_nodes[j] = false : set_of_nodes[j] = true;
                }

                v[i] = 0;

                for (int j(0); j < v.size() - 1; ++j)
                {
                    // Calculate the node with the shortest distance to the active set of nodes
                    Index = VectorMaskedMinIndex<>::value(v, set_of_nodes);

                    // Delete the current node from the active set of nodes
                    set_of_nodes[Index] = false;

                    // Update the graph distances from node i to each other node by using the costs of the current node
                    for (int n(0); n < v.size(); ++n)
                    {
                        help1 = v[Index] + cost_matrix[n][Index];
                        if (set_of_nodes[n] == true && v[n] > help1)
                        {
                            v[n] = help1;
                            p[n] = Index;
                        }
                    }
                }
            }

            // Write the graph distances from node i to each other node represented by graph_distance_matrix into the cost_matrix
            for (typename MutableMatrix<DataType_>::ElementIterator e(graph_distance_matrix.begin_elements()),
                    e_end(graph_distance_matrix.end_elements()), f(cost_matrix.begin_elements()); e != e_end ; ++e, ++f)
                    {
                        *f = *e;
                    }
        }
    };
}
#endif
