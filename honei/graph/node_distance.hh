/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Sven Mallach <sven.mallach@uni-dortmund.de>
 * Copyright (c) 2007 Nina Harmuth <nina.harmuth@uni-dortmund.de>
 * Copyright (c) 2007 Thorsten Deinert <thorsten.deinert@uni-dortmund.de>
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

#pragma once
#ifndef LIBGRAPH_GUARD_NODE_DISTANCE_HH
#define LIBGRAPH_GUARD_NODE_DISTANCE_HH 1

#include <honei/la/banded_matrix.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/difference.hh>
#include <honei/la/matrix_error.hh>
#include <honei/la/norm.hh>
#include <honei/util/tags.hh>
#include <honei/graph/graph_error.hh>


/**
 * \file
 *
 * Implementation of NodeDistance.
 *
 * \ingroup grplibgraph
 **/
namespace honei
{
    template <typename Tag_ = tags::CPU> struct NodeDistance;
    /**
     * \brief NodeDistance is used in the algorithm of Kamada-Kawai.
     * \brief Hardly surprising, NodeDistance computes the distance between nodes,
     * \brief which positions are given in a position matrix.
     *
     * \ingroup grplibgraph
     **/
    template <> struct NodeDistance<tags::CPU>
    {
        /**
         * Returns the resulting distance matrix.
         * \param pos_matrix The matrix with the positions of nodes. It may have as many rows as needed for all dimensions
         */
        template <typename DataType_>
        static DenseMatrix<DataType_> value(const DenseMatrix<DataType_> & pos_matrix)
        {
            CONTEXT("When calculating the distance beetween nodes:");
            // Create the result matrix
            DenseMatrix<DataType_> result(pos_matrix.rows(), pos_matrix.rows());

            // Retrieving an ElementIterator e initialized to the first element in the result matrix
            typename DenseMatrix<DataType_>::ElementIterator e(result.begin_elements());

            // Iterate over all nodes (=columns) in the given position matrix
            for (unsigned long i =  0; i < pos_matrix.rows(); ++i)
            {
                // The column-vector represents all n coordinates of node i, so we grab them from the matrix
                DenseVectorRange<DataType_> v(pos_matrix[i]);

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
            return result;
        }

        template <typename DataType_>
        static void value(const DenseMatrix<DataType_> & pos_matrix, const DenseMatrix<DataType_> & edge_weights,
        DenseMatrix<DataType_> & square_dist, DenseMatrix<DataType_> & inv_square_dist,
        const DataType_ repulsive_force_range)
        {
            if (edge_weights.rows() != edge_weights.columns())
                throw MatrixIsNotSquare(edge_weights.rows(), edge_weights.columns());

            if (square_dist.rows() != square_dist.columns())
                throw MatrixIsNotSquare(square_dist.rows(), square_dist.columns());

            if (inv_square_dist.rows() != inv_square_dist.columns())
                throw MatrixIsNotSquare(inv_square_dist.rows(), inv_square_dist.columns());

            if (edge_weights.rows() != pos_matrix.rows())
                throw MatrixRowsDoNotMatch(pos_matrix.rows(), edge_weights.rows());

            if (square_dist.rows() != pos_matrix.rows())
                throw MatrixRowsDoNotMatch(pos_matrix.rows(), square_dist.rows());

            if (inv_square_dist.rows() != pos_matrix.rows())
                throw MatrixRowsDoNotMatch(pos_matrix.rows(), inv_square_dist.rows());

            if (repulsive_force_range < 0)
                throw GraphError("Repulsive-Force-Range must be positiv");

            // Initialize ElementIterators
            typename DenseMatrix<DataType_>::ElementIterator e(inv_square_dist.begin_elements());
            typename DenseMatrix<DataType_>::ElementIterator f(square_dist.begin_elements());
            typename DenseMatrix<DataType_>::ConstElementIterator g(edge_weights.begin_elements());
            DataType_ square_force_range(repulsive_force_range * repulsive_force_range);

            // Iterate over all nodes (=columns) in the given position matrix
            for (unsigned long i =  0; i < pos_matrix.rows(); ++i)
            {
                // The column-vector represents all n coordinates of node i, so we grab them from the matrix
                DenseVectorRange<DataType_> v(pos_matrix[i]);

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



    template <> struct NodeDistance<tags::CPU::SSE>
    {
        static DenseMatrix<float> value(const DenseMatrix<float> & pos_matrix);
        static DenseMatrix<double> value(const DenseMatrix<double> & pos_matrix);

        static void value(const DenseMatrix<float> & pos_matrix, const DenseMatrix<float> & edge_weights,
                DenseMatrix<float> & square_dist, DenseMatrix<float> & inv_square_dist,
                const float repulsive_force_range);

        static void value(const DenseMatrix<double> & pos_matrix, const DenseMatrix<double> & edge_weights,
                DenseMatrix<double> & square_dist, DenseMatrix<double> & inv_square_dist,
                const double repulsive_force_range);

        static void node_distance(float * result, float * pos_matrix, float x, float y, unsigned long size);
        static void node_distance(double * result, double * pos_matrix, double x, double y, unsigned long size);

        static void set_distance(float * dist, float * dist_src, float * mask, unsigned long size);
        static void set_distance(double * dist, double * dist_src, double * mask, unsigned long size);

        static void invert_distance(float * dist, float square_force_range, unsigned long size);
        static void invert_distance(double * dist, double square_force_range, unsigned long size);
    };

    template <> struct NodeDistance<tags::Cell>
    {
        static DenseMatrix<float> value(const DenseMatrix<float> & pos_matrix);

        static void value(const DenseMatrix<float> & pos_matrix, const DenseMatrix<float> & edge_weights,
                DenseMatrix<float> & square_dist, DenseMatrix<float> & inv_square_dist,
                const float repulsive_force_range);
    };

    template <> struct NodeDistance<tags::CPU::MultiCore> : public NodeDistance<tags::CPU>{};
    template <> struct NodeDistance<tags::CPU::MultiCore::SSE> : public NodeDistance<tags::CPU>{};
}
#endif
