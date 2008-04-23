/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the LA C++ library. LibLa is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibLa is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/la/element_product.hh>

#include <xmmintrin.h>
#include <emmintrin.h>

namespace honei
{
    namespace intern
    {
        namespace sse
        {
            inline void element_product(float * a, const float * b, unsigned long size)
            {
                __m128 m1, m2, m3, m4, m5, m6, m7, m8;

                unsigned long a_address = reinterpret_cast<unsigned long>(a);
                unsigned long a_offset = a_address % 16;
                unsigned long b_address = reinterpret_cast<unsigned long>(b);
                unsigned long b_offset = b_address % 16;

                unsigned long x_offset(a_offset / 4);
                x_offset = (4 - x_offset) % 4;

                unsigned long quad_start = x_offset;
                unsigned long quad_end(size - ((size - quad_start) % 16));

                if (size < 24)
                {
                    quad_end = 0;
                    quad_start = 0;
                }

                if (a_offset == b_offset)
                {
                    for (unsigned long index(quad_start) ; index < quad_end ; index += 16)
                    {
                        m1 = _mm_load_ps(a + index);
                        m3 = _mm_load_ps(a + index + 4);
                        m5 = _mm_load_ps(a + index + 8);
                        m7 = _mm_load_ps(a + index + 12);
                        m2 = _mm_load_ps(b + index);
                        m4 = _mm_load_ps(b + index + 4);
                        m6 = _mm_load_ps(b + index + 8);
                        m8 = _mm_load_ps(b + index + 12);

                        m1 = _mm_mul_ps(m1, m2);
                        m3 = _mm_mul_ps(m3, m4);
                        m5 = _mm_mul_ps(m5, m6);
                        m7 = _mm_mul_ps(m7, m8);

                        _mm_store_ps(a + index, m1);
                        _mm_store_ps(a + index + 4, m3);
                        _mm_store_ps(a + index + 8, m5);
                        _mm_store_ps(a + index + 12, m7);
                    }
                }
                else
                {
                    for (unsigned long index(quad_start) ; index < quad_end ; index += 16)
                    {
                        m1 = _mm_load_ps(a + index);
                        m3 = _mm_load_ps(a + index + 4);
                        m5 = _mm_load_ps(a + index + 8);
                        m7 = _mm_load_ps(a + index + 12);
                        m2 = _mm_loadu_ps(b + index);
                        m4 = _mm_loadu_ps(b + index + 4);
                        m6 = _mm_loadu_ps(b + index + 8);
                        m8 = _mm_loadu_ps(b + index + 12);

                        m1 = _mm_mul_ps(m1, m2);
                        m3 = _mm_mul_ps(m3, m4);
                        m5 = _mm_mul_ps(m5, m6);
                        m7 = _mm_mul_ps(m7, m8);

                        _mm_store_ps(a + index, m1);
                        _mm_store_ps(a + index + 4, m3);
                        _mm_store_ps(a + index + 8, m5);
                        _mm_store_ps(a + index + 12, m7);
                    }
                }

                for (unsigned long index(0) ; index < quad_start ; index++)
                {
                    a[index] *= b[index];
                }

                for (unsigned long index(quad_end) ; index < size ; index++)
                {
                    a[index] *= b[index];
                }
            }

            inline void element_product(double * a, const double * b, unsigned long size)
            {
                __m128d m1, m2, m3, m4, m5, m6, m7, m8;

                unsigned long a_address = reinterpret_cast<unsigned long>(a);
                unsigned long a_offset = a_address % 16;
                unsigned long b_address = reinterpret_cast<unsigned long>(b);
                unsigned long b_offset = b_address % 16;

                unsigned long x_offset(a_offset / 8);

                unsigned long quad_start = x_offset;
                unsigned long quad_end(size - ((size - quad_start) % 8));

                if (size < 16)
                {
                    quad_end = 0;
                    quad_start = 0;
                }

                if (a_offset == b_offset)
                {
                    for (unsigned long index(quad_start) ; index < quad_end ; index += 8)
                    {
                        m1 = _mm_load_pd(a + index);
                        m3 = _mm_load_pd(a + index + 2);
                        m5 = _mm_load_pd(a + index + 4);
                        m7 = _mm_load_pd(a + index + 6);
                        m2 = _mm_load_pd(b + index);
                        m4 = _mm_load_pd(b + index + 2);
                        m6 = _mm_load_pd(b + index + 4);
                        m8 = _mm_load_pd(b + index + 6);

                        m1 = _mm_mul_pd(m1, m2);
                        m3 = _mm_mul_pd(m3, m4);
                        m5 = _mm_mul_pd(m5, m6);
                        m7 = _mm_mul_pd(m7, m8);

                        _mm_store_pd(a + index, m1);
                        _mm_store_pd(a + index + 2, m3);
                        _mm_store_pd(a + index + 4, m5);
                        _mm_store_pd(a + index + 6, m7);
                    }
                }
                else
                {
                    for (unsigned long index(quad_start) ; index < quad_end ; index += 8)
                    {
                        m1 = _mm_load_pd(a + index);
                        m3 = _mm_load_pd(a + index + 2);
                        m5 = _mm_load_pd(a + index + 4);
                        m7 = _mm_load_pd(a + index + 6);
                        m2 = _mm_loadu_pd(b + index);
                        m4 = _mm_loadu_pd(b + index + 2);
                        m6 = _mm_loadu_pd(b + index + 4);
                        m8 = _mm_loadu_pd(b + index + 6);

                        m1 = _mm_mul_pd(m1, m2);
                        m3 = _mm_mul_pd(m3, m4);
                        m5 = _mm_mul_pd(m5, m6);
                        m7 = _mm_mul_pd(m7, m8);

                        _mm_store_pd(a + index, m1);
                        _mm_store_pd(a + index + 2, m3);
                        _mm_store_pd(a + index + 4, m5);
                        _mm_store_pd(a + index + 6, m7);
                    }
                }

                for (unsigned long index(0) ; index < quad_start ; index++)
                {
                    a[index] *= b[index];
                }

                for (unsigned long index(quad_end) ; index < size ; index++)
                {
                    a[index] *= b[index];
                }
            }
        }
    }
}

using namespace honei;

DenseVectorContinuousBase<float> & ElementProduct<tags::CPU::SSE>::value(DenseVectorContinuousBase<float> & a,
        const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When multiplying DenseVectorContinuousBase<float> and DenseVectorContinuousBase<float> elementwise "
            "with SSE:");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    intern::sse::element_product(a.elements(), b.elements(), a.size());

    return a;
}

DenseVectorContinuousBase<double> & ElementProduct<tags::CPU::SSE>::value(DenseVectorContinuousBase<double> & a,
        const DenseVectorContinuousBase<double> & b)
{
    CONTEXT("When multiplying DenseVectorContinuousBase<double> and DenseVectorContinuousBase<double> elementwise "
            "with SSE:");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    intern::sse::element_product(a.elements(), b.elements(), a.size());

    return a;
}

DenseMatrix<float> & ElementProduct<tags::CPU::SSE>::value(DenseMatrix<float> & a, const DenseMatrix<float> & b)
{
    CONTEXT("When multiplying DenseMatrix<float> and DenseMatrix<float> elementwise with SSE:");

    if (a.columns() != b.columns())
    {
        throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
    }

    if (a.rows() != b.rows())
    {
        throw MatrixRowsDoNotMatch(b.rows(), a.rows());
    }

    intern::sse::element_product(a.elements(), b.elements(), a.rows() * a.columns());

    return a;
}

DenseMatrix<double> & ElementProduct<tags::CPU::SSE>::value(DenseMatrix<double> & a, const DenseMatrix<double> & b)
{
    CONTEXT("When multiplying DenseMatrix<double> and DenseMatrix<double> elementwise with SSE:");

    if (a.columns() != b.columns())
    {
        throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
    }

    if (a.rows() != b.rows())
    {
        throw MatrixRowsDoNotMatch(b.rows(), a.rows());
    }

    intern::sse::element_product(a.elements(), b.elements(), a.rows() * a.columns());

    return a;
}

