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

#ifndef LIBGRAPH_GUARD_NODE_DISTANCE_MC_HH
#define LIBGRAPH_GUARD_NODE_DISTANCE_MC_HH 1

#include <honei/libla/dense_matrix_tile.hh>
#include <honei/libutil/wrapper.hh>
#include <honei/libutil/thread_pool.hh>
#include <honei/libgraph/graph_error.hh>
#include <honei/libgraph/node_distance.hh>
#include <honei/libutil/tags.hh>
#include <iostream>


// For optimization purposes
#define MIN_CHUNK_SIZE 256
#define PARTS 8

namespace honei
{
    /**
     * \brief NodeDistance is used in the algorithm of Kamada-Kawai.
     * \brief Hardly surprising, NodeDistance computes the distance between nodes,
     * \brief which positions are given in a position matrix.
     *
     * \ingroup grplibgraph
     **/

    template <typename Tag_> struct NodeDistance;

    namespace intern
    {
        template <typename DataType_>
        static inline DataType_ calculate_parts(const DataType_ ref)
        {
            if (ref / 16 < 8) return ref / 16;
            return 8;
        }
    }

    template <typename Tag_> struct MCNodeDistance
    {
        template <typename DataType_>
        static DenseMatrix<DataType_> value(const DenseMatrix<DataType_> & pos_matrix)
        {
            // Create the result
            DenseMatrix<DataType_> result(pos_matrix.rows(), pos_matrix.rows(), DataType_(0));

            ThreadPool * tp(ThreadPool::instance());

            std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

            unsigned long num_parts(intern::calculate_parts(result.rows()));

            if (num_parts)
            {
                unsigned long rows(result.rows() / num_parts);
                unsigned long rest_of_rows(result.rows() % num_parts);
                unsigned long matrix_row(0);
                unsigned long i(0);

                for ( ; i < num_parts - 1; ++i)
                {
                    DenseMatrixTile<DataType_> tile(result, rows, result.columns(), matrix_row, 0);
                    ThreeArgWrapper< MCNodeDistance<Tag_>, DenseMatrixTile<DataType_>, const DenseMatrix<DataType_>, unsigned long>
                    wrapper(tile, pos_matrix, matrix_row);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(wrapper));
                    dispatched_tasks.push_back(ptr);
                    matrix_row += rows;
                }
                DenseMatrixTile<DataType_> tile(result, rows + rest_of_rows, result.columns(), matrix_row, 0);
                ThreeArgWrapper< MCNodeDistance<Tag_>, DenseMatrixTile<DataType_>, const DenseMatrix<DataType_>, unsigned long>
                wrapper(tile, pos_matrix, matrix_row);
                std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(wrapper));
                dispatched_tasks.push_back(ptr);
            }
            else
            {
                result = NodeDistance<typename Tag_::DelegateTo>::value(pos_matrix);
            }

            while(! dispatched_tasks.empty())
            {
                dispatched_tasks.front()->wait_on();
                dispatched_tasks.pop_front();
            }

            return result;
        }

        template <typename DataType_>
        static void value(const DenseMatrix<DataType_> & pos_matrix, const DenseMatrix<DataType_> & edge_weights,
        DenseMatrix<DataType_> & square_dist, DenseMatrix<DataType_> & inv_square_dist,
        const DataType_ repulsive_force_range)
        {
            if (! edge_weights.square())
                throw MatrixIsNotSquare(edge_weights.rows(), edge_weights.columns());

            if (! square_dist.square())
                throw MatrixIsNotSquare(square_dist.rows(), square_dist.columns());

            if (! inv_square_dist.square())
                throw MatrixIsNotSquare(inv_square_dist.rows(), inv_square_dist.columns());

            if (edge_weights.rows() != pos_matrix.rows())
                throw MatrixRowsDoNotMatch(pos_matrix.rows(), edge_weights.rows());

            if (square_dist.rows() != pos_matrix.rows())
                throw MatrixRowsDoNotMatch(pos_matrix.rows(), square_dist.rows());

            if (inv_square_dist.rows() != pos_matrix.rows())
                throw MatrixRowsDoNotMatch(pos_matrix.rows(), inv_square_dist.rows());

            if (repulsive_force_range < 0)
                throw GraphError("Repulsive-Force-Range must be positiv");

            ThreadPool * tp(ThreadPool::instance());

            std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

            unsigned long num_parts(intern::calculate_parts(edge_weights.rows()));

            if (num_parts)
            {
                unsigned long rows(edge_weights.rows() / num_parts);
                unsigned long rest_of_rows(edge_weights.rows() % num_parts);
                unsigned long matrix_row(0);
                unsigned long i(0);

                for ( ; i < num_parts - 1; ++i)
                {
                    DenseMatrixTile<DataType_> tile_1(edge_weights, rows, edge_weights.columns(), matrix_row, 0);
                    DenseMatrixTile<DataType_> tile_2(square_dist, rows, square_dist.columns(), matrix_row, 0);
                    DenseMatrixTile<DataType_> tile_3(inv_square_dist, rows, inv_square_dist.columns(), matrix_row, 0);
                    SixArgWrapper< MCNodeDistance<Tag_>, const DenseMatrix<DataType_>, const DenseMatrixTile<DataType_>,
                    DenseMatrixTile<DataType_>, DenseMatrixTile<DataType_>, const DataType_, unsigned long>
                    wrapper(pos_matrix, tile_1, tile_2, tile_3, repulsive_force_range, matrix_row);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(wrapper));
                    dispatched_tasks.push_back(ptr);
                    matrix_row += rows;
                }
                DenseMatrixTile<DataType_> tile_1(edge_weights, rows + rest_of_rows, edge_weights.columns(), matrix_row, 0);
                DenseMatrixTile<DataType_> tile_2(square_dist, rows + rest_of_rows, square_dist.columns(), matrix_row, 0);
                DenseMatrixTile<DataType_> tile_3(inv_square_dist, rows + rest_of_rows, inv_square_dist.columns(), matrix_row, 0);
                SixArgWrapper< MCNodeDistance<Tag_>, const DenseMatrix<DataType_>, const DenseMatrixTile<DataType_>,
                DenseMatrixTile<DataType_>, DenseMatrixTile<DataType_>, const DataType_, unsigned long>
                wrapper(pos_matrix, tile_1, tile_2, tile_3, repulsive_force_range, matrix_row);
                std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(wrapper));
                dispatched_tasks.push_back(ptr);
            }
            else
            {
                NodeDistance<typename Tag_::DelegateTo>::value(pos_matrix, edge_weights, square_dist, inv_square_dist, repulsive_force_range);
            }

            while(! dispatched_tasks.empty())
            {
                dispatched_tasks.front()->wait_on();
                dispatched_tasks.pop_front();
            }
        }

        template <typename DataType_>
        static void value(DenseMatrixTile<DataType_> & result, const DenseMatrix<DataType_> & pos_matrix, unsigned long start)
        {
            CONTEXT("When calculating the distance beetween nodes (with MC):");

            // Retrieving an ElementIterator e initialized to the first element in the result matrix
            typename DenseMatrixTile<DataType_>::ElementIterator e(result.begin_elements());

            // Iterate over all nodes (=columns) in the given position matrix
            for (unsigned long i = 0; i < result.rows(); ++i)
            {
                // The column-vector represents all n coordinates of node i, so we grab them from the matrix
                DenseVectorRange<DataType_> v(pos_matrix[i + start]);

                // Now, iterate over all nodes to calculate the distance between node i and node j. For the resulting
                // distance matrix is symmetric, result(i,j) equals result(j,i) and so were waiting for symmetric matrices
                // to gain performance. So long, the trivial loops.
                for (unsigned long j = 0; j < pos_matrix.rows(); ++j, ++e)
                {
                    DenseVectorRange<DataType_> w(pos_matrix[j]);
                    // Now calculate difference tmp = v - w for each pair of nodes and
                    // then, calculate d = tmp1^2 + tmp2^2 + ... + tmpN^2 which is the
                    // l2-norm of tmp without a root. (see template parameter "root")
                    DenseVector<DataType_> v_copy(v.copy());
                    DataType_ d(Norm<vnt_l_two, false, tags::CPU>::value(Difference<tags::CPU>::value(v_copy, w)));
                    *e = d;
                }
            }
        }

        template <typename DataType_>
        static void value(const DenseMatrix<DataType_> & pos_matrix, const DenseMatrixTile<DataType_> & edge_weights,
        DenseMatrixTile<DataType_> & square_dist, DenseMatrixTile<DataType_> & inv_square_dist,
        const DataType_ repulsive_force_range, unsigned long start)
        {
            CONTEXT("When calculating the distance beetween nodes (with MC):");

            // Initialize ElementIterators
            typename DenseMatrixTile<DataType_>::ElementIterator e(inv_square_dist.begin_elements());
            typename DenseMatrixTile<DataType_>::ElementIterator f(square_dist.begin_elements());
            typename DenseMatrixTile<DataType_>::ConstElementIterator g(edge_weights.begin_elements());
            DataType_ square_force_range(repulsive_force_range * repulsive_force_range);

            // Iterate over all nodes (=columns) in the given position matrix
            for (unsigned long i = 0; i < edge_weights.rows(); ++i)
            {
                // The column-vector represents all n coordinates of node i, so we grab them from the matrix
                DenseVectorRange<DataType_> v(pos_matrix[i + start]);

                // Now, iterate over all nodes to calculate the distance between node i and node j. For the resulting
                // distance matrices are symmetric, and so we are waiting for symmetric matrices
                // to gain performance. So long, the trivial loops.
                for (unsigned long j = 0; j < pos_matrix.rows(); ++j, ++e, ++f, ++g)
                {
                    DenseVectorRange<DataType_> w(pos_matrix[j]);
                    // Now calculate difference tmp = v - w for each pair of nodes and
                    // then, calculate d = tmp1^2 + tmp2^2 + ... + tmpN^2 which is the
                    // l2-norm of tmp without a root. (see template parameter "root")
                    DenseVector<DataType_> v_copy(v.copy());
                    DataType_ d(Norm<vnt_l_two, false, tags::CPU>::value(Difference<tags::CPU>::value(v_copy, w)));
                    (d < square_force_range) && (d > std::numeric_limits<DataType_>::epsilon()) ? *e = 1 / d : *e = 0;
                    if (*g > std::numeric_limits<DataType_>::epsilon()) *f = d;
                }
            }
        }
    };

    template <typename Tag_> struct MCSSENodeDistance
    {
        template <typename DataType_>
        static DenseMatrix<DataType_> value(const DenseMatrix<DataType_> & pos_matrix)
        {
            // Create the result
            DenseMatrix<DataType_> result(pos_matrix.rows(), pos_matrix.rows(), DataType_(0));

            ThreadPool * tp(ThreadPool::instance());

            std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

            unsigned long num_parts(intern::calculate_parts(result.rows()));

            if (num_parts)
            {
                unsigned long rows(result.rows() / num_parts);
                unsigned long rest_of_rows(result.rows() % num_parts);
                unsigned long matrix_row(0);
                unsigned long i(0);

                for ( ; i < num_parts - 1; ++i)
                {
                    DenseMatrixTile<DataType_> tile(result, rows, result.columns(), matrix_row, 0);
                    ThreeArgWrapper< MCSSENodeDistance<Tag_>, DenseMatrixTile<DataType_>, const DenseMatrix<DataType_>, unsigned long>
                    wrapper(tile, pos_matrix, matrix_row);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(wrapper));
                    dispatched_tasks.push_back(ptr);
                    matrix_row += rows;
                }
                DenseMatrixTile<DataType_> tile(result, rows + rest_of_rows, result.columns(), matrix_row, 0);
                ThreeArgWrapper< MCSSENodeDistance<Tag_>, DenseMatrixTile<DataType_>, const DenseMatrix<DataType_>, unsigned long>
                wrapper(tile, pos_matrix, matrix_row);
                std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(wrapper));
                dispatched_tasks.push_back(ptr);
            }
            else
            {
                result = NodeDistance<typename Tag_::DelegateTo>::value(pos_matrix);
            }

            while(! dispatched_tasks.empty())
            {
                dispatched_tasks.front()->wait_on();
                dispatched_tasks.pop_front();
            }

            return result;
        }

        template <typename DataType_>
        static void value(const DenseMatrix<DataType_> & pos_matrix, const DenseMatrix<DataType_> & edge_weights,
        DenseMatrix<DataType_> & square_dist, DenseMatrix<DataType_> & inv_square_dist,
        const DataType_ repulsive_force_range)
        {
            if (! edge_weights.square())
                throw MatrixIsNotSquare(edge_weights.rows(), edge_weights.columns());

            if (! square_dist.square())
                throw MatrixIsNotSquare(square_dist.rows(), square_dist.columns());

            if (! inv_square_dist.square())
                throw MatrixIsNotSquare(inv_square_dist.rows(), inv_square_dist.columns());

            if (edge_weights.rows() != pos_matrix.rows())
                throw MatrixRowsDoNotMatch(pos_matrix.rows(), edge_weights.rows());

            if (square_dist.rows() != pos_matrix.rows())
                throw MatrixRowsDoNotMatch(pos_matrix.rows(), square_dist.rows());

            if (inv_square_dist.rows() != pos_matrix.rows())
                throw MatrixRowsDoNotMatch(pos_matrix.rows(), inv_square_dist.rows());

            if (repulsive_force_range < 0)
                throw GraphError("Repulsive-Force-Range must be positiv");

            ThreadPool * tp(ThreadPool::instance());

            std::list< std::tr1::shared_ptr<PoolTask> > dispatched_tasks;

            unsigned long num_parts(intern::calculate_parts(edge_weights.rows()));

            if (num_parts)
            {
                unsigned long rows(edge_weights.rows() / num_parts);
                unsigned long rest_of_rows(edge_weights.rows() % num_parts);
                unsigned long matrix_row(0);
                unsigned long i(0);

                for ( ; i < num_parts - 1; ++i)
                {
                    DenseMatrixTile<DataType_> tile_1(edge_weights, rows, edge_weights.columns(), matrix_row, 0);
                    DenseMatrixTile<DataType_> tile_2(square_dist, rows, square_dist.columns(), matrix_row, 0);
                    DenseMatrixTile<DataType_> tile_3(inv_square_dist, rows, inv_square_dist.columns(), matrix_row, 0);
                    SixArgWrapper< MCSSENodeDistance<Tag_>, const DenseMatrix<DataType_>, const DenseMatrixTile<DataType_>,
                    DenseMatrixTile<DataType_>, DenseMatrixTile<DataType_>, const DataType_, unsigned long>
                    wrapper(pos_matrix, tile_1, tile_2, tile_3, repulsive_force_range, matrix_row);
                    std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(wrapper));
                    dispatched_tasks.push_back(ptr);
                    matrix_row += rows;
                }
                DenseMatrixTile<DataType_> tile_1(edge_weights, rows + rest_of_rows, edge_weights.columns(), matrix_row, 0);
                DenseMatrixTile<DataType_> tile_2(square_dist, rows + rest_of_rows, square_dist.columns(), matrix_row, 0);
                DenseMatrixTile<DataType_> tile_3(inv_square_dist, rows + rest_of_rows, inv_square_dist.columns(), matrix_row, 0);
                SixArgWrapper< MCSSENodeDistance<Tag_>, const DenseMatrix<DataType_>, const DenseMatrixTile<DataType_>,
                DenseMatrixTile<DataType_>, DenseMatrixTile<DataType_>, const DataType_, unsigned long>
                wrapper(pos_matrix, tile_1, tile_2, tile_3, repulsive_force_range, matrix_row);
                std::tr1::shared_ptr<PoolTask> ptr(ThreadPool::instance()->dispatch(wrapper));
                dispatched_tasks.push_back(ptr);
            }
            else
            {
                NodeDistance<typename Tag_::DelegateTo>::value(pos_matrix, edge_weights, square_dist, inv_square_dist, repulsive_force_range);
            }

            while(! dispatched_tasks.empty())
            {
                dispatched_tasks.front()->wait_on();
                dispatched_tasks.pop_front();
            }
        }

        template <typename DataType_>
        static void value(DenseMatrixTile<DataType_> & result, const DenseMatrix<DataType_> & pos_matrix, const unsigned long start)
        {
            CONTEXT("When calculating the distance beetween nodes (with MC-SSE):");

            // Iterate over all nodes (=columns) in the given position matrix
            for (unsigned long i = 0, j = result.columns(); i < result.rows(); ++i)
            {
                NodeDistance<typename Tag_::DelegateTo>::node_distance(result[i].elements(), pos_matrix.elements(), pos_matrix(i  + start, 0), pos_matrix(i + start, 1), j);
            }
        }

        template <typename DataType_>
        static void value(const DenseMatrix<DataType_> & pos_matrix, const DenseMatrixTile<DataType_> & edge_weights,
        DenseMatrixTile<DataType_> & square_dist, DenseMatrixTile<DataType_> & inv_square_dist,
        const DataType_ repulsive_force_range, unsigned long start)
        {
           CONTEXT("When calculating the distance and the inverse distance between nodes (with MC-SSE):");

            // Iterate over all nodes (=columns) in the given position matrix
            for (unsigned long i = 0, j = inv_square_dist.columns(); i < inv_square_dist.rows(); ++i)
            {
                NodeDistance<typename Tag_::DelegateTo>::node_distance(inv_square_dist[i].elements(), pos_matrix.elements(), pos_matrix(i + start, 0), pos_matrix(i + start, 1), j);
            }

            NodeDistance<typename Tag_::DelegateTo>::set_distance(square_dist.elements(), inv_square_dist.elements(), edge_weights.elements(), square_dist.columns() * square_dist.rows());

            NodeDistance<typename Tag_::DelegateTo>::invert_distance(inv_square_dist.elements(), repulsive_force_range, inv_square_dist.columns() * inv_square_dist.rows());
        }
    };
}
#endif
