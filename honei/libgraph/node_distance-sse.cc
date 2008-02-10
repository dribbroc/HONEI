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

#include <honei/libutil/attributes.hh>
#include <honei/libgraph/node_distance.hh>

#if defined (__SSE3__)
#include <pmmintrin.h>
#else
#include <emmintrin.h>
#endif

namespace honei
{
    namespace intern
    {
        namespace sse
        {
#if defined (__SSE3__)
            inline void node_distance(float * result, float * pos_matrix, float x, float y, unsigned long size)
            {
                __m128 m1, m2, m3, m4, m5, m6, m7, m8;
                float HONEI_ALIGNED(16) coordinates[4];
                coordinates[0] = x;
                coordinates[1] = y;
                coordinates[2] = x;
                coordinates[3] = y;
                m8 = _mm_load_ps(coordinates);

                unsigned long result_address = (unsigned long)result;
                unsigned long result_offset = result_address % 16;

                unsigned long x_offset(result_offset / 4);
                x_offset = (4 - x_offset) % 4;

                unsigned long quad_start = x_offset;
                unsigned long quad_end(size - ((size - quad_start) % 12));

                if (size < 24)
                {
                    quad_end = 0;
                    quad_start = 0;
                }

                for (unsigned long index(0) ; index < quad_start ; ++index)
                {
                    result[index] = ((x - pos_matrix[index * 2]) * (x - pos_matrix[index * 2])
                            +  ((y - pos_matrix[index * 2 + 1]) * (y - pos_matrix[index * 2 + 1])));
                }
                for (unsigned long index(quad_start) ; index < quad_end ; index += 12)
                {
                    m1 = _mm_load_ps(pos_matrix + (index * 2));
                    m2 = _mm_load_ps(pos_matrix + (index * 2) + 4);
                    m3 = _mm_load_ps(pos_matrix + (index * 2) + 8);
                    m4 = _mm_load_ps(pos_matrix + (index * 2) + 12);
                    m5 = _mm_load_ps(pos_matrix + (index * 2) + 16);
                    m6 = _mm_load_ps(pos_matrix + (index * 2) + 20);
                    m1 = _mm_sub_ps(m8, m1);
                    m2 = _mm_sub_ps(m8, m2);
                    m3 = _mm_sub_ps(m8, m3);
                    m4 = _mm_sub_ps(m8, m4);
                    m5 = _mm_sub_ps(m8, m5);
                    m6 = _mm_sub_ps(m8, m6);
                    m1 = _mm_mul_ps(m1, m1);
                    m2 = _mm_mul_ps(m2, m2);
                    m3 = _mm_mul_ps(m3, m3);
                    m4 = _mm_mul_ps(m4, m4);
                    m5 = _mm_mul_ps(m5, m5);
                    m6 = _mm_mul_ps(m6, m6);
                    m1 = _mm_hadd_ps(m1, m2);
                    m3 = _mm_hadd_ps(m3, m4);
                    m5 = _mm_hadd_ps(m5, m6);
                    _mm_store_ps(result + index, m1);
                    _mm_store_ps(result + index + 4, m3);
                    _mm_store_ps(result + index + 8, m5);
                }
                for (unsigned long index(quad_end) ; index < size ; ++index)
                {
                    result[index] = ((x - pos_matrix[index * 2]) * (x - pos_matrix[index * 2])
                            +  ((y - pos_matrix[index * 2 + 1]) * (y - pos_matrix[index * 2 + 1])));
                }
            }

            inline void node_distance(double * result, double * pos_matrix, double x, double y, unsigned long size)
            {
                __m128d m1, m2, m3, m4, m5, m6, m8;
                double HONEI_ALIGNED(16) coordinates[4];
                coordinates[0] = x;
                coordinates[1] = y;
                coordinates[2] = x;
                coordinates[3] = y;
                m8 = _mm_load_pd(coordinates);

                unsigned long result_address = (unsigned long)result;
                unsigned long result_offset = result_address % 16;

                unsigned long x_offset(result_offset / 8);

                unsigned long quad_start = x_offset;
                unsigned long quad_end(size - ((size - quad_start) % 6));

                if (size < 24)
                {
                    quad_end = 0;
                    quad_start = 0;
                }

                for (unsigned long index(0) ; index < quad_start ; ++index)
                {
                    result[index] = ((x - pos_matrix[index * 2]) * (x - pos_matrix[index * 2])
                            +  ((y - pos_matrix[index * 2 + 1]) * (y - pos_matrix[index * 2 + 1])));
                }
                for (unsigned long index(quad_start) ; index < quad_end ; index += 6)
                {
                    m1 = _mm_load_pd(pos_matrix + (index * 2));
                    m2 = _mm_load_pd(pos_matrix + (index * 2) + 2);
                    m3 = _mm_load_pd(pos_matrix + (index * 2) + 4);
                    m4 = _mm_load_pd(pos_matrix + (index * 2) + 6);
                    m5 = _mm_load_pd(pos_matrix + (index * 2) + 8);
                    m6 = _mm_load_pd(pos_matrix + (index * 2) + 10);
                    m1 = _mm_sub_pd(m8, m1);
                    m2 = _mm_sub_pd(m8, m2);
                    m3 = _mm_sub_pd(m8, m3);
                    m4 = _mm_sub_pd(m8, m4);
                    m5 = _mm_sub_pd(m8, m5);
                    m6 = _mm_sub_pd(m8, m6);
                    m1 = _mm_mul_pd(m1, m1);
                    m2 = _mm_mul_pd(m2, m2);
                    m3 = _mm_mul_pd(m3, m3);
                    m4 = _mm_mul_pd(m4, m4);
                    m5 = _mm_mul_pd(m5, m5);
                    m6 = _mm_mul_pd(m6, m6);
                    m1 = _mm_hadd_pd(m1, m2);
                    m3 = _mm_hadd_pd(m3, m4);
                    m5 = _mm_hadd_pd(m5, m6);
                    _mm_store_pd(result + index, m1);
                    _mm_store_pd(result + index + 2, m3);
                    _mm_store_pd(result + index + 4, m5);
                }
                for (unsigned long index(quad_end) ; index < size ; ++index)
                {
                    result[index] = ((x - pos_matrix[index * 2]) * (x - pos_matrix[index * 2])
                            +  ((y - pos_matrix[index * 2 + 1]) * (y - pos_matrix[index * 2 + 1])));
                }
            }
#else
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
#endif
            inline void set_distance(float * dist, float * dist_src, float * mask, unsigned long size)
            {
                __m128 m1, m2, m3, m4, m5, m6, m7, m8;
                float HONEI_ALIGNED(16) eps(std::numeric_limits<float>::epsilon());
                m8 = _mm_load1_ps(&eps);

                unsigned long dist_address = (unsigned long)dist;
                unsigned long dist_offset = dist_address % 16;

                unsigned long x_offset(dist_offset / 4);
                x_offset = (4 - x_offset) % 4;

                unsigned long quad_start = x_offset;
                unsigned long quad_end(size - ((size - quad_start) % 12));

                if (size < 24)
                {
                    quad_end = 0;
                    quad_start = 0;
                }

                for (unsigned long index(0) ; index < quad_start ; ++index)
                {
                    if (mask[index] > std::numeric_limits<float>::epsilon())
                            dist[index] = dist_src[index];
                }
                for (unsigned long index(quad_start) ; index < quad_end; index += 12)
                {
                    m1 = _mm_load_ps(dist_src + index);
                    m2 = _mm_load_ps(mask + index);
                    m3 = _mm_load_ps(dist_src + index + 4);
                    m4 = _mm_load_ps(mask + index + 4);
                    m5 = _mm_load_ps(dist_src + index + 8);
                    m6 = _mm_load_ps(mask + index + 8);

                    m7 = _mm_cmpgt_ps(m2, m8);
                    m1 = _mm_and_ps(m1, m7);
                    m7 = _mm_cmpgt_ps(m4, m8);
                    m3 = _mm_and_ps(m3, m7);
                    m7 = _mm_cmpgt_ps(m6, m8);
                    m5 = _mm_and_ps(m5, m7);

                    _mm_store_ps(dist + index, m1);
                    _mm_store_ps(dist + index + 4, m3);
                    _mm_store_ps(dist + index + 8, m5);
                }
                for (unsigned long index(quad_end) ; index < size ; ++index)
                {
                    if (mask[index] > std::numeric_limits<float>::epsilon())
                            dist[index] = dist_src[index];
                }
            }

            inline void set_distance(double * dist, double * dist_src, double * mask, unsigned long size)
            {
                __m128d m1, m2, m3, m4, m5, m6, m7, m8;
                double HONEI_ALIGNED(16) eps(std::numeric_limits<double>::epsilon());
                m8 = _mm_load1_pd(&eps);

                unsigned long dist_address = (unsigned long)dist;
                unsigned long dist_offset = dist_address % 16;

                unsigned long x_offset(dist_offset / 8);

                unsigned long quad_start = x_offset;
                unsigned long quad_end(size - ((size - quad_start) % 6));

                if (size < 24)
                {
                    quad_end = 0;
                    quad_start = 0;
                }

                for (unsigned long index(0) ; index < quad_start ; ++index)
                {
                    if (mask[index] > std::numeric_limits<double>::epsilon())
                            dist[index] = dist_src[index];
                }
                for (unsigned long index(quad_start) ; index < quad_end; index += 6)
                {
                    m1 = _mm_load_pd(dist_src + index);
                    m2 = _mm_load_pd(mask + index);
                    m3 = _mm_load_pd(dist_src + index + 2);
                    m4 = _mm_load_pd(mask + index + 2);
                    m5 = _mm_load_pd(dist_src + index + 4);
                    m6 = _mm_load_pd(mask + index + 4);

                    m7 = _mm_cmpgt_pd(m2, m8);
                    m1 = _mm_and_pd(m1, m7);
                    m7 = _mm_cmpgt_pd(m4, m8);
                    m3 = _mm_and_pd(m3, m7);
                    m7 = _mm_cmpgt_pd(m6, m8);
                    m5 = _mm_and_pd(m5, m7);

                    _mm_store_pd(dist + index, m1);
                    _mm_store_pd(dist + index + 2, m3);
                    _mm_store_pd(dist + index + 4, m5);
                }
                for (unsigned long index(quad_end) ; index < size ; ++index)
                {
                    if (mask[index] > std::numeric_limits<double>::epsilon())
                            dist[index] = dist_src[index];
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

void NodeDistance<tags::CPU::SSE>::value(const DenseMatrix<float> & pos_matrix, const DenseMatrix<float> & neighbours,
        DenseMatrix<float> & square_dist, DenseMatrix<float> & inv_square_dist, const float repulsive_force_range)
{
    float square_force_range(repulsive_force_range * repulsive_force_range);

    for (unsigned long i(0); i < pos_matrix.rows(); ++i)
    {
        intern::sse::node_distance(inv_square_dist[i].elements(), pos_matrix.elements(), pos_matrix(i, 0), pos_matrix(i, 1), inv_square_dist[i].size());
    }

    intern::sse::set_distance(square_dist.elements(), inv_square_dist.elements(), neighbours.elements(), square_dist.columns() * square_dist.rows());

    intern::sse::invert_distance(inv_square_dist.elements(), square_force_range, inv_square_dist.columns() * inv_square_dist.rows());
}

void NodeDistance<tags::CPU::SSE>::value(const DenseMatrix<double> & pos_matrix, const DenseMatrix<double> & neighbours,
        DenseMatrix<double> & square_dist, DenseMatrix<double> & inv_square_dist, const double repulsive_force_range)
{
    double square_force_range(repulsive_force_range * repulsive_force_range);

    for (unsigned long i(0); i < pos_matrix.rows(); ++i)
    {
        intern::sse::node_distance(inv_square_dist[i].elements(), pos_matrix.elements(), pos_matrix(i, 0), pos_matrix(i, 1), inv_square_dist[i].size());
    }

    intern::sse::set_distance(square_dist.elements(), inv_square_dist.elements(), neighbours.elements(), square_dist.columns() * square_dist.rows());

    intern::sse::invert_distance(inv_square_dist.elements(), square_force_range, inv_square_dist.columns() * inv_square_dist.rows());
}
