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

#include <libutil/tags.hh>
#include <libla/dense_matrix.hh>
#include <libla/sparse_matrix.hh>
#include <libla/dense_vector.hh>
#include <libla/dense_vector-impl.hh>
#include <libla/matrix_error.hh>
#include <libla/vector_error.hh>
#include <libgraph/graph_error.hh>
#include <libgraph/abstract_graph.hh>
#include <libgraph/breadth_first_search-mc.hh>

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
         * Returns the resulting graph distance matrix.
         * \param distance_matrix Empty matrix, which will contains the result.
         * \param adjacency_matrix Which nodes are neighbours?.
         */

        template <typename DataType_>
        static bool value(DenseMatrix<DataType_> & distance_matrix, const SparseMatrix<bool> & adjacency_matrix)
        {
            if (! adjacency_matrix.square())
                throw MatrixIsNotSquare(adjacency_matrix.rows(), adjacency_matrix.columns());

            bool symmetric(true);
            bool trace(true);
            for (typename Matrix<bool>::ConstElementIterator e(adjacency_matrix.begin_non_zero_elements()),
                        e_end(adjacency_matrix.end_non_zero_elements()); e != e_end ; ++e)
            {
                if (*e != adjacency_matrix(e.column(), e.row())) symmetric = false;
                if (e.column() == e.row()) trace = false;
            }

            if (! symmetric)
                throw GraphError("Adjacency matrix must be symmetric");

            if (! trace)
                throw GraphError("Adjacency matrix has non vanishing trace");

            if (! distance_matrix.square())
                throw MatrixIsNotSquare(distance_matrix.rows(), distance_matrix.columns());

            if (adjacency_matrix.rows() != distance_matrix.rows())
                throw MatrixRowsDoNotMatch(distance_matrix.rows(), adjacency_matrix.rows());

            // Create the result
            bool result(true);

            // Number of nodes watched
            unsigned long number_of_nodes(0);

            // Visited nodes
            DenseMatrix<bool> visited_nodes(adjacency_matrix.rows(), adjacency_matrix.columns(), false);

            // Node queue
            std::queue<unsigned long> node_queue;

            // Start BFS for each node of the graph to calculate all distances between this starting-node and all other nodes
            for (unsigned long i(0); i < visited_nodes.rows(); i++)
            {
                //Insert current starting-node into the queue; queue defines the order in which all nodes have to be visited
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
                    for (typename Vector<bool>::ConstElementIterator e((adjacency_matrix[current_node]).begin_non_zero_elements()),
                    e_end((adjacency_matrix[current_node]).end_non_zero_elements()); e != e_end ; ++e)
                    {
                        // Check that node has not been visited yet
                        if (!visited_nodes(i, e.index()))
                        {
                            // Insert current node (at the end of the queue)
                            node_queue.push(e.index());

                            // Mark node as visited
                            visited_nodes(i, e.index()) = true;

                            // Calculate the distance between the starting-node and the current node by dint of its predecessor
                            distance_matrix(i, e.index()) = distance_matrix(i, current_node) + 1;
                        }
                    }
                }

                if (i == 0 && number_of_nodes != adjacency_matrix.rows())
                {
                    result = false;
                }
            }

            return result;
        }

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
            if (! edge_weights.square())
                throw MatrixIsNotSquare(edge_weights.rows(), edge_weights.columns());

            bool symmetric(true);
            bool trace(true);
            for (typename Matrix<DataType_>::ConstElementIterator e(edge_weights.begin_non_zero_elements()),
                e_end(edge_weights.end_non_zero_elements()); e != e_end ; ++e)
            {
                if ((*e - edge_weights(e.column(), e.row()) > std::numeric_limits<DataType_>::epsilon())
                || (*e - edge_weights(e.column(), e.row()) < - std::numeric_limits<DataType_>::epsilon())) symmetric = false;
                if (e.column() == e.row()) trace = false;
            }

            if (! symmetric)
                throw GraphError("Edge weights matrix must be symmetric");

            if (! trace)
                throw GraphError("Edge weights matrix has non vanishing trace");

            if (! distance_matrix.square())
                throw MatrixIsNotSquare(distance_matrix.rows(), distance_matrix.columns());

            if (edge_weights.rows() != distance_matrix.rows())
                throw MatrixRowsDoNotMatch(distance_matrix.rows(), edge_weights.rows());

            if (node_weights.size() != distance_matrix.rows())
                throw VectorSizeDoesNotMatch(node_weights.size(), distance_matrix.rows());

            bool invalid(true);
            for (typename Vector<DataType_>::ConstElementIterator e(node_weights.begin_elements()),
                e_end(node_weights.end_elements()); e != e_end ; ++e)
            {
                if (*e < std::numeric_limits<DataType_>::epsilon()) invalid = false;
            }

            if (! invalid)
                throw GraphError("Node weights vector contained invalid value");

            // Create the result
            bool result(true);

            // Number of nodes watched
            unsigned long number_of_nodes(0);

            // Visited nodes
            DenseMatrix<bool> visited_nodes(edge_weights.rows(), edge_weights.columns(), false);

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
        const SparseMatrix<DataType_> & edge_weights, const AbstractGraph<DataType_> & graph)
        {
            if (! edge_weights.square())
                throw MatrixIsNotSquare(edge_weights.rows(), edge_weights.columns());

            bool symmetric(true);
            bool trace(true);
            for (typename Matrix<DataType_>::ConstElementIterator e(edge_weights.begin_non_zero_elements()),
                e_end(edge_weights.end_non_zero_elements()); e != e_end ; ++e)
            {
                if ((*e - edge_weights(e.column(), e.row()) > std::numeric_limits<DataType_>::epsilon())
                || (*e - edge_weights(e.column(), e.row()) < - std::numeric_limits<DataType_>::epsilon())) symmetric = false;
                if (e.column() == e.row()) trace = false;
            }

            if (! symmetric)
                throw GraphError("Edge weights matrix must be symmetric");

            if (! trace)
                throw GraphError("Edge weights matrix has non vanishing trace");

            if (! distance_matrix.square())
                throw MatrixIsNotSquare(distance_matrix.rows(), distance_matrix.columns());

            if (edge_weights.rows() != distance_matrix.rows())
                throw MatrixRowsDoNotMatch(distance_matrix.rows(), edge_weights.rows());

            if (node_weights.size() != distance_matrix.rows())
                throw VectorSizeDoesNotMatch(node_weights.size(), distance_matrix.rows());

            bool invalid(true);
            for (typename Vector<DataType_>::ConstElementIterator e(node_weights.begin_elements()),
                e_end(node_weights.end_elements()); e != e_end ; ++e)
            {
                if (*e < std::numeric_limits<DataType_>::epsilon()) invalid = false;
            }

            if (! invalid)
                throw GraphError("Node weights vector contained invalid value");

            // Create the result
            bool result(true);

            // Number of nodes watched
            unsigned long number_of_nodes(0);

            // Visited nodes
            DenseMatrix<bool> visited_nodes(edge_weights.rows(), edge_weights.columns(), false);

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
                            if (graph.sameTimeslice(i, e.index()))
                            {
                                distance_matrix(i, e.index()) = distance_matrix(i, current_node) +
                                (sqrt(node_weights[current_node] * node_weights[e.index()]) / *e);
                            }
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

    template <> struct BreadthFirstSearch<tags::CPU::MultiCore> : public MCBreadthFirstSearch<tags::CPU::MultiCore>{};
}
#endif
