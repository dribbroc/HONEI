/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Mathias Kadolsky <mathkad@gmx.de>
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

#ifndef LIBGRAPH_GUARD_BREADTHFRISTSEARCH_HH
#define LIBGRAPH_GUARD_BREADTHFRISTSEARCH_HH 1

#include <honei/util/tags.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/dense_vector-impl.hh>
#include <honei/la/matrix_error.hh>
#include <honei/la/vector_error.hh>
#include <honei/graph/graph_error.hh>
#include <honei/graph/abstract_graph.hh>

#include <queue>
#include <iostream>

namespace honei
{
    /**
     * \brief BFS computes the graph distances between nodes.
     *
     * \ingroup grplibgraph
     **/


    template <typename Tag_ = tags::CPU> struct BreadthFirstSearch;

    template <> struct BreadthFirstSearch<tags::CPU>
    {
        /**
         * Computes the resulting weighted graph distance matrix.
         * \param distance_matrix Empty matrix, which will contains the result.
         * \param node_weights Vector with the weights of nodes.
         * \param edge_weights Matrix with the weights of edges.
         */
        template <typename DataType_>
        static bool value(DenseMatrix<DataType_> & distance_matrix, const DenseVector<DataType_> & node_weights,
        const SparseMatrix<DataType_> & edge_weights)
        {
            if (edge_weights.rows() != edge_weights.columns())
                throw MatrixIsNotSquare(edge_weights.rows(), edge_weights.columns());

            if (distance_matrix.rows() != distance_matrix.columns())
                throw MatrixIsNotSquare(distance_matrix.rows(), distance_matrix.columns());

            if (edge_weights.rows() != distance_matrix.rows())
                throw MatrixRowsDoNotMatch(distance_matrix.rows(), edge_weights.rows());

            if (node_weights.size() != distance_matrix.rows())
                throw VectorSizeDoesNotMatch(node_weights.size(), distance_matrix.rows());

            bool invalid(true);
            for (typename DenseVector<DataType_>::ConstElementIterator e(node_weights.begin_elements()),
                e_end(node_weights.end_elements()); e != e_end ; ++e)
            {
                if (*e <= DataType_(0)) invalid = false;
            }

            if (! invalid)
                throw GraphError("Node weights vector contained invalid value");

            // Create the result
            bool result(true);

            // Number of nodes watched
            unsigned long number_of_nodes(0);

            // Visited nodes
            DenseMatrix<float> visited_nodes(edge_weights.rows(), edge_weights.columns(), false);

            // Node queue
            std::queue<unsigned long> node_queue;

            // Start BFS for each node of the graph to calculate all distances between this starting-node and all other nodes
            for (unsigned long i(0); i < visited_nodes.rows(); i++)
            {
                // Insert current starting-node into the queue; queue defines the order in which all nodes have to be visited
                node_queue.push(i);
                visited_nodes(i, i) = true;
                distance_matrix(i, i) = DataType_(0);

                while (!node_queue.empty())
                {
                    // Return the first element of the queue
                    unsigned long current_node(node_queue.front());

                    // Delete the first element of the queue
                    node_queue.pop();

                    // Auxiliary variable to check if the graph is coherent
                    number_of_nodes++;

                    // Select all nodes that are adjacent to the starting-node
                    for (typename Vector<DataType_>::ConstElementIterator e((edge_weights[current_node]).begin_non_zero_elements()),
                    e_end((edge_weights[current_node]).end_non_zero_elements()); e != e_end ; ++e)
                    {
                        // Check that node has not been visited yet
                        if (!visited_nodes(i, e.index()))
                        {
                            // Insert current node (at the end of the queue)
                            node_queue.push(e.index());

                            // Mark node as visited
                            visited_nodes(i, e.index()) = true;

                            // Calculate the distance between the starting-node and the current node by dint of its predecessor, involving the weights of the edges
                            distance_matrix(i, e.index()) = distance_matrix(i, current_node) +
                            (sqrt(node_weights[current_node] * node_weights[e.index()]) / *e);
                        }
                    }
                }

                if (i == 0 && number_of_nodes != edge_weights.rows())
                {
                    result = false;
                }
            }

            return result;
        }

        template <typename DataType_>
        static bool value(DenseMatrix<DataType_> & distance_matrix, const DenseVector<DataType_> & node_weights,
        const SparseMatrix<DataType_> & edge_weights, AbstractGraph<DataType_> & graph)
        {
            if (edge_weights.rows() != edge_weights.columns())
                throw MatrixIsNotSquare(edge_weights.rows(), edge_weights.columns());

            if (distance_matrix.rows() != distance_matrix.columns())
                throw MatrixIsNotSquare(distance_matrix.rows(), distance_matrix.columns());

            if (edge_weights.rows() != distance_matrix.rows())
                throw MatrixRowsDoNotMatch(distance_matrix.rows(), edge_weights.rows());

            if (node_weights.size() != distance_matrix.rows())
                throw VectorSizeDoesNotMatch(node_weights.size(), distance_matrix.rows());

            bool invalid(true);
            for (typename DenseVector<DataType_>::ConstElementIterator e(node_weights.begin_elements()),
                e_end(node_weights.end_elements()); e != e_end ; ++e)
            {
                if (*e <= DataType_(0)) invalid = false;
            }

            if (! invalid)
                throw GraphError("Node weights vector contained invalid value");

            // Create the result
            bool result(true);

            // Number of nodes watched
            unsigned long number_of_nodes(0);

            // Visited nodes
            DenseMatrix<float> visited_nodes(edge_weights.rows(), edge_weights.columns(), false);

            // Node queue
            std::queue<unsigned long> node_queue;

            // Start BFS for each node of the graph to calculate all distances between this starting-node and all other nodes
            for (unsigned long i(0); i < visited_nodes.rows(); i++)
            {
                // Insert current starting-node into the queue; queue defines the order in which all nodes have to be visited
                node_queue.push(i);
                visited_nodes(i, i) = true;
                distance_matrix(i, i) = DataType_(0);

                while (!node_queue.empty()) 
                {
                    // Return the first element of the queue
                    unsigned long current_node(node_queue.front());

                    // Delete the first element of the queue
                    node_queue.pop();

                    // Auxiliary variable to check if the graph is coherent
                    number_of_nodes++;

                    // Select all nodes that are adjacent to the starting-node
                    for (typename Vector<DataType_>::ConstElementIterator e((edge_weights[current_node]).begin_non_zero_elements()),
                    e_end((edge_weights[current_node]).end_non_zero_elements()); e != e_end ; ++e)
                    {
                        // Check that node has not been visited yet
                        if (!visited_nodes(i, e.index()))
                        {
                            // Insert current node (at the end of the queue)
                            node_queue.push(e.index());

                            // Mark node as visited
                            visited_nodes(i, e.index()) = true;

                            // Calculate the distance between the starting-node and the current node by dint of its predecessor, involving the weights of the edges
                            distance_matrix(i, e.index()) = distance_matrix(i, current_node) +
                            (graph.same_timeslice(current_node, e.index()) ?
                            (sqrt(node_weights[current_node] * node_weights[e.index()]) / *e):
                            0);
                        }
                    }
                }

                if (i == 0 && number_of_nodes != edge_weights.rows())
                {
                    result = false;
                }
            }

            return result;
        }
    };

    template <> struct BreadthFirstSearch<tags::CPU::SSE> : public BreadthFirstSearch<tags::CPU>{};
    template <> struct BreadthFirstSearch<tags::Cell> : public BreadthFirstSearch<tags::CPU>{};

    template <> struct BreadthFirstSearch<tags::CPU::MultiCore> : public BreadthFirstSearch<tags::CPU>{};
    template <> struct BreadthFirstSearch<tags::CPU::MultiCore::SSE> : public BreadthFirstSearch<tags::CPU>{};

}
#endif
