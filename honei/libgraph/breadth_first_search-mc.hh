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

#ifndef LIBGRAPH_GUARD_BREADTHFRISTSEARCH_MC_HH
#define LIBGRAPH_GUARD_BREADTHFRISTSEARCH_MC_HH 1

#include <honei/libla/dense_matrix_tile.hh>
#include <honei/util/wrapper.hh>
#include <honei/util/thread_pool.hh>
#include <honei/libgraph/breadth_first_search.hh>

#include <queue>

// For optimization purposes
#define MIN_CHUNK_SIZE 256
#define PARTS 8

namespace honei
{
    /**
     * \brief BFS computes the graph distances between nodes.
     *
     * \ingroup grplibgraph
     **/

    template <typename Tag_ > struct BreadthFirstSearch;

    template <typename Tag_> struct MCBreadthFirstSearch
    {
        template <typename DataType_>
        static inline DataType_ calculate_parts(const DataType_ ref)
        {
            if (ref / 16 < PARTS) return ref / 16;
            return PARTS;
        }

        template <typename DataType_>
        static bool value(DenseMatrix<DataType_> & distance_matrix, const DenseVector<DataType_> & node_weights,
        const SparseMatrix<DataType_> & edge_weights)
        {
            if (! edge_weights.square())
                throw MatrixIsNotSquare(edge_weights.rows(), edge_weights.columns());

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
                if (*e <= DataType_(0)) invalid = false;
            }

            if (! invalid)
                throw GraphError("Node weights vector contained invalid value");

            // Create the result
            bool result(true);

            ThreadPool * tp(ThreadPool::instance());

            std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

            unsigned long num_parts(MCBreadthFirstSearch<Tag_>::calculate_parts(distance_matrix.rows()));

            if (num_parts)
            {
                unsigned long rows(distance_matrix.rows() / num_parts);
                unsigned long rest_of_rows(distance_matrix.rows() % num_parts);
                unsigned long matrix_row(0);
                unsigned long i(0);

                for ( ; i < num_parts - 1; ++i)
                {
                    DenseMatrixTile<DataType_> tile(distance_matrix, rows, distance_matrix.columns(), matrix_row, 0);
                    ResultFourArgWrapper< MCBreadthFirstSearch<Tag_>, bool, DenseMatrixTile<DataType_>, const DenseVector<DataType_>, 
                    const SparseMatrix<DataType_>, unsigned long> 
                    wrapper(result, tile, node_weights, edge_weights, matrix_row);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(wrapper));
                    dispatched_tasks.push_back(ptr);
                    matrix_row += rows;
                }
                DenseMatrixTile<DataType_> tile(distance_matrix, rows + rest_of_rows, distance_matrix.columns(), matrix_row, 0);
                ResultFourArgWrapper< MCBreadthFirstSearch<Tag_>, bool, DenseMatrixTile<DataType_>, const DenseVector<DataType_>, 
                const SparseMatrix<DataType_>, unsigned long> 
                wrapper(result, tile, node_weights, edge_weights, matrix_row);
                std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(wrapper));
                dispatched_tasks.push_back(ptr);
            }
            else
            {
                result = BreadthFirstSearch<typename Tag_::DelegateTo>::value(distance_matrix, node_weights, edge_weights);
            }

            while(! dispatched_tasks.empty())
            {
                dispatched_tasks.front()->wait_on();
                dispatched_tasks.pop_front();
            }

            return result;
        }

        template <typename DataType_>
        static bool value(DenseMatrix<DataType_> & distance_matrix, const DenseVector<DataType_> & node_weights,
        const SparseMatrix<DataType_> & edge_weights, AbstractGraph<DataType_> & graph)
        {
            if (! edge_weights.square())
                throw MatrixIsNotSquare(edge_weights.rows(), edge_weights.columns());

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
                if (*e <= DataType_(0)) invalid = false;
            }

            if (! invalid)
                throw GraphError("Node weights vector contained invalid value");

            // Create the result
            bool result(true);

            ThreadPool * tp(ThreadPool::instance());

            std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

            unsigned long num_parts(MCBreadthFirstSearch<Tag_>::calculate_parts(distance_matrix.rows()));

            if (num_parts)
            {
                unsigned long rows(distance_matrix.rows() / num_parts);
                unsigned long rest_of_rows(distance_matrix.rows() % num_parts);
                unsigned long matrix_row(0);
                unsigned long i(0);


                for ( ; i < num_parts - 1; ++i)
                {
                    DenseMatrixTile<DataType_> tile(distance_matrix, rows, distance_matrix.columns(), matrix_row, 0);
                    ResultFiveArgWrapper< MCBreadthFirstSearch<Tag_>, bool, DenseMatrixTile<DataType_>, const DenseVector<DataType_>, 
                    const SparseMatrix<DataType_>, const AbstractGraph<DataType_>, unsigned long> 
                    wrapper(result, tile, node_weights, edge_weights, graph, matrix_row);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(wrapper));
                    dispatched_tasks.push_back(ptr);
                    matrix_row += rows;
                }
                DenseMatrixTile<DataType_> tile(distance_matrix, rows + rest_of_rows, distance_matrix.columns(), matrix_row, 0);
                ResultFiveArgWrapper< MCBreadthFirstSearch<Tag_>, bool, DenseMatrixTile<DataType_>, const DenseVector<DataType_>, 
                const SparseMatrix<DataType_>, const AbstractGraph<DataType_>, unsigned long> 
                wrapper(result, tile, node_weights, edge_weights, graph, matrix_row);
                std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(wrapper));
                dispatched_tasks.push_back(ptr);
            }
            else
            {
                result = BreadthFirstSearch<typename Tag_::DelegateTo>::value(distance_matrix, node_weights, edge_weights, graph);
            }

            while(! dispatched_tasks.empty())
            {
                dispatched_tasks.front()->wait_on();
                dispatched_tasks.pop_front();
            }

            return result;
        }

        template <typename DataType_>
        static bool value(DenseMatrixTile<DataType_> & distance_matrix, const DenseVector<DataType_> & node_weights,
        const SparseMatrix<DataType_> & edge_weights, unsigned long start)
        {
            // Create the result
            bool result(true);

            // Number of nodes watched
            unsigned long number_of_nodes(0);

            // Node queue
            std::queue<unsigned long> node_queue;

            // Current index of distance matrix
            unsigned long current_index(0);

            // Visited nodes
            DenseMatrix<bool> visited_nodes(distance_matrix.rows(), edge_weights.columns(), false);

            // Start BFS for each node of the graph to calculate all distances between this starting-node and all other nodes
            for (unsigned long i(0); i < visited_nodes.rows(); i++)
            {
                //Insert current starting-node into the queue; queue defines the order in which all nodes have to be visited
                current_index = i + start;
                node_queue.push(current_index);
                visited_nodes(i, current_index) = true;
                distance_matrix(i, current_index) = DataType_(0);

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

                            // Calculate the distance between the starting-node and the current node by dint of its predecessor
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
        static bool value(DenseMatrixTile<DataType_> & distance_matrix, const DenseVector<DataType_> & node_weights,
        const SparseMatrix<DataType_> & edge_weights, AbstractGraph<DataType_> & graph, unsigned long start)
        {
            // Create the result
            bool result(true);

            // Number of nodes watched
            unsigned long number_of_nodes(0);

            // Node queue
            std::queue<unsigned long> node_queue;

            // Current index of distance matrix
            unsigned long current_index(0);

            // Visited nodes
            DenseMatrix<bool> visited_nodes(distance_matrix.rows(), edge_weights.columns(), false);

            // Start BFS for each node of the graph to calculate all distances between this starting-node and all other nodes
            for (unsigned long i(0); i < visited_nodes.rows(); i++)
            {
                //Insert current starting-node into the queue; queue defines the order in which all nodes have to be visited
                current_index = i + start;
                node_queue.push(current_index);
                visited_nodes(i, current_index) = true;
                distance_matrix(i, current_index) = DataType_(0);

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

                            // Calculate the distance between the starting-node and the current node by dint of its predecessor
                            distance_matrix(i, e.index()) = distance_matrix(i, current_node) +
                            (graph.sameTimeslice(current_node, e.index()) ?
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
}
#endif
