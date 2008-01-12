/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include "node_distance.hh"

namespace honei
{
    namespace intern
    {
        namespace sse
        {
            inline void node_distance(float * result, float * pos_matrix, float x, float y, unsigned long size)
            {
                for (unsigned long row(0) ; row < size ; ++row)
                {
                    result[row] = ((x - pos_matrix[row * 2]) * (x - pos_matrix[row * 2])
                            +  ((y - pos_matrix[row * 2 + 1]) * (y - pos_matrix[row * 2 + 1])));
                }
            }

            inline void node_distance(double * result, double * pos_matrix, double x, double y, unsigned long size)
            {
                for (unsigned long row(0) ; row < size ; ++row)
                {
                    result[row] = ((x - pos_matrix[row * 2]) * (x - pos_matrix[row * 2])
                            +  ((y - pos_matrix[row * 2 + 1]) * (y - pos_matrix[row * 2 + 1])));
                }
            }

            inline void invert_distance(float * dist, float square_force_range, unsigned long size)
            {
                for (unsigned long index(0) ; index < size ; ++index)
                {
                    if (dist[index] < square_force_range && dist[index] > std::numeric_limits<float>::epsilon())
                        dist[index] = float(1) / dist[index];
                    else
                        dist[index] = float(0);
                }
            }

            inline void invert_distance(double * dist, double square_force_range, unsigned long size)
            {
                for (unsigned long index(0) ; index < size ; ++index)
                {
                    if (dist[index] < square_force_range && dist[index] > std::numeric_limits<double>::epsilon())
                        dist[index] = float(1) / dist[index];
                    else
                        dist[index] = float(0);
                }
            }
        }
    }
}

using namespace honei;

DenseMatrix<float> NodeDistance<tags::CPU::SSE>::value(const DenseMatrix<float> & pos_matrix)
{
    CONTEXT("When calculating the distance beetween nodes (float) (with SSE):");
    if (pos_matrix.columns() != 2)
        throw MatrixColumnsDoNotMatch(2, pos_matrix.columns());

    DenseMatrix<float> result(pos_matrix.rows(), pos_matrix.rows());

    for (unsigned long i(0); i < pos_matrix.rows(); ++i)
    {
        intern::sse::node_distance(result[i].elements(), pos_matrix.elements(), pos_matrix(i, 0), pos_matrix(i, 1), result[i].size());
    }

    return result;
}

DenseMatrix<double> NodeDistance<tags::CPU::SSE>::value(const DenseMatrix<double> & pos_matrix)
{
    CONTEXT("When calculating the distance beetween nodes (double) (with SSE):");
    if (pos_matrix.columns() != 2)
        throw MatrixColumnsDoNotMatch(2, pos_matrix.columns());

    DenseMatrix<double> result(pos_matrix.rows(), pos_matrix.rows());

    for (unsigned long i(0); i < pos_matrix.rows(); ++i)
    {
        intern::sse::node_distance(result[i].elements(), pos_matrix.elements(), pos_matrix(i, 0), pos_matrix(i, 1), result[i].size());
    }

    return result;
}

void NodeDistance<tags::CPU::SSE>::value(const DenseMatrix<float> & pos_matrix, const SparseMatrix<bool> & neighbours,
        SparseMatrix<float> & square_dist, DenseMatrix<float> & inv_square_dist, const float repulsive_force_range)
{
    float square_force_range(repulsive_force_range * repulsive_force_range);

    for (unsigned long i(0); i < pos_matrix.rows(); ++i)
    {
        intern::sse::node_distance(inv_square_dist[i].elements(), pos_matrix.elements(), pos_matrix(i, 0), pos_matrix(i, 1), inv_square_dist[i].size());
    }

    for(SparseMatrix<bool>::ConstElementIterator g(neighbours.begin_non_zero_elements()), g_end(neighbours.end_non_zero_elements()) ;
            g != g_end ; ++g)
    {
        if (*g != false)
            square_dist(g.row(), g.column()) = inv_square_dist(g.row(), g.column());
    }

    intern::sse::invert_distance(inv_square_dist.elements(), square_force_range, inv_square_dist.columns() * inv_square_dist.rows());
}

void NodeDistance<tags::CPU::SSE>::value(const DenseMatrix<double> & pos_matrix, const SparseMatrix<bool> & neighbours,
        SparseMatrix<double> & square_dist, DenseMatrix<double> & inv_square_dist, const double repulsive_force_range)
{
    double square_force_range(repulsive_force_range * repulsive_force_range);

    for (unsigned long i(0); i < pos_matrix.rows(); ++i)
    {
        intern::sse::node_distance(inv_square_dist[i].elements(), pos_matrix.elements(), pos_matrix(i, 0), pos_matrix(i, 1), inv_square_dist[i].size());
    }

    for(SparseMatrix<bool>::ConstElementIterator g(neighbours.begin_non_zero_elements()), g_end(neighbours.end_non_zero_elements()) ;
            g != g_end ; ++g)
    {
        if (*g != false)
            square_dist(g.row(), g.column()) = inv_square_dist(g.row(), g.column());
    }

    intern::sse::invert_distance(inv_square_dist.elements(), square_force_range, inv_square_dist.columns() * inv_square_dist.rows());
}

void NodeDistance<tags::CPU::SSE>::value(const DenseMatrix<float> & pos_matrix, const SparseMatrix<float> & neighbours,
        SparseMatrix<float> & square_dist, DenseMatrix<float> & inv_square_dist, const float repulsive_force_range)
{
    float square_force_range(repulsive_force_range * repulsive_force_range);

    for (unsigned long i(0); i < pos_matrix.rows(); ++i)
    {
        intern::sse::node_distance(inv_square_dist[i].elements(), pos_matrix.elements(), pos_matrix(i, 0), pos_matrix(i, 1), inv_square_dist[i].size());
    }

    for(SparseMatrix<float>::ConstElementIterator g(neighbours.begin_non_zero_elements()), g_end(neighbours.end_non_zero_elements()) ;
            g != g_end ; ++g)
    {
        if (*g > std::numeric_limits<float>::epsilon())
            square_dist(g.row(), g.column()) = inv_square_dist(g.row(), g.column());
    }

    intern::sse::invert_distance(inv_square_dist.elements(), square_force_range, inv_square_dist.columns() * inv_square_dist.rows());
}

void NodeDistance<tags::CPU::SSE>::value(const DenseMatrix<double> & pos_matrix, const SparseMatrix<double> & neighbours,
        SparseMatrix<double> & square_dist, DenseMatrix<double> & inv_square_dist, const double repulsive_force_range)
{
    double square_force_range(repulsive_force_range * repulsive_force_range);

    for (unsigned long i(0); i < pos_matrix.rows(); ++i)
    {
        intern::sse::node_distance(inv_square_dist[i].elements(), pos_matrix.elements(), pos_matrix(i, 0), pos_matrix(i, 1), inv_square_dist[i].size());
    }

    for(SparseMatrix<double>::ConstElementIterator g(neighbours.begin_non_zero_elements()), g_end(neighbours.end_non_zero_elements()) ;
            g != g_end ; ++g)
    {
        if (*g > std::numeric_limits<double>::epsilon())
            square_dist(g.row(), g.column()) = inv_square_dist(g.row(), g.column());
    }

    intern::sse::invert_distance(inv_square_dist.elements(), square_force_range, inv_square_dist.columns() * inv_square_dist.rows());
}
