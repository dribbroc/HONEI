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

#include <libla/sum.hh>

#include <xmmintrin.h>
#include <emmintrin.h>

namespace honei
{
    namespace intern
    {
        namespace sse
        {
            inline void sum(float * a, const float * b, unsigned long size)
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

                        m1 = _mm_add_ps(m1, m2);
                        m3 = _mm_add_ps(m3, m4);
                        m5 = _mm_add_ps(m5, m6);
                        m7 = _mm_add_ps(m7, m8);

                        _mm_stream_ps(a + index, m1);
                        _mm_stream_ps(a + index + 4, m3);
                        _mm_stream_ps(a + index + 8, m5);
                        _mm_stream_ps(a + index + 12, m7);
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

                        m1 = _mm_add_ps(m1, m2);
                        m3 = _mm_add_ps(m3, m4);
                        m5 = _mm_add_ps(m5, m6);
                        m7 = _mm_add_ps(m7, m8);

                        _mm_stream_ps(a + index, m1);
                        _mm_stream_ps(a + index + 4, m3);
                        _mm_stream_ps(a + index + 8, m5);
                        _mm_stream_ps(a + index + 12, m7);
                    }
                }

                for (unsigned long index(0) ; index < quad_start ; index++)
                {
                    a[index] += b[index];
                }

                for (unsigned long index(quad_end) ; index < size ; index++)
                {
                    a[index] += b[index];
                }

                _mm_sfence();
            }

            inline void sum(double * a, const double * b, unsigned long size)
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

                        m1 = _mm_add_pd(m1, m2);
                        m3 = _mm_add_pd(m3, m4);
                        m5 = _mm_add_pd(m5, m6);
                        m7 = _mm_add_pd(m7, m8);

                        _mm_stream_pd(a + index, m1);
                        _mm_stream_pd(a + index + 2, m3);
                        _mm_stream_pd(a + index + 4, m5);
                        _mm_stream_pd(a + index + 6, m7);
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

                        m1 = _mm_add_pd(m1, m2);
                        m3 = _mm_add_pd(m3, m4);
                        m5 = _mm_add_pd(m5, m6);
                        m7 = _mm_add_pd(m7, m8);

                        _mm_stream_pd(a + index, m1);
                        _mm_stream_pd(a + index + 2, m3);
                        _mm_stream_pd(a + index + 4, m5);
                        _mm_stream_pd(a + index + 6, m7);
                    }
                }

                for (unsigned long index(0) ; index < quad_start ; index++)
                {
                    a[index] += b[index];
                }

                for (unsigned long index(quad_end) ; index < size ; index++)
                {
                    a[index] += b[index];
                }

                _mm_sfence();
            }

            inline void sum(const float a, float * x, unsigned long size)
            {
                __m128 m1, m2, m3, m4, m5, m6, m7, m8;
                float __attribute__((aligned(16))) a_data;
                a_data= a;
                m8 = _mm_load_ps1(&a_data);

                unsigned long x_address = reinterpret_cast<unsigned long>(x);
                unsigned long x_offset = x_address % 16;

                unsigned long z_offset(x_offset / 4);
                z_offset = (4 - z_offset) % 4;

                unsigned long quad_start = z_offset;
                unsigned long quad_end(size - ((size - quad_start) % 28));

                if (size < 36)
                {
                    quad_end = 0;
                    quad_start = 0;
                }

                for (unsigned long index(quad_start) ; index < quad_end ; index += 28)
                {
                    m1 = _mm_load_ps(x + index);
                    m2 = _mm_load_ps(x + index + 4);
                    m3 = _mm_load_ps(x + index + 8);
                    m4 = _mm_load_ps(x + index + 12);
                    m5 = _mm_load_ps(x + index + 16);
                    m6 = _mm_load_ps(x + index + 20);
                    m7 = _mm_load_ps(x + index + 24);

                    m1 = _mm_add_ps(m1, m8);
                    m2 = _mm_add_ps(m2, m8);
                    m3 = _mm_add_ps(m3, m8);
                    m4 = _mm_add_ps(m4, m8);
                    m5 = _mm_add_ps(m5, m8);
                    m6 = _mm_add_ps(m6, m8);
                    m7 = _mm_add_ps(m7, m8);

                    _mm_stream_ps(x + index, m1);
                    _mm_stream_ps(x + index + 4, m2);
                    _mm_stream_ps(x + index + 8, m3);
                    _mm_stream_ps(x + index + 12, m4);
                    _mm_stream_ps(x + index + 16, m5);
                    _mm_stream_ps(x + index + 20, m6);
                    _mm_stream_ps(x + index + 24, m7);
                }

                for (unsigned long index(0) ; index < quad_start ; index++)
                {
                    x[index] += a;
                }

                for (unsigned long index(quad_end) ; index < size ; index++)
                {
                    x[index] += a;
                }

                _mm_sfence();
            }

            inline void sum(const double a, double * x, unsigned long size)
            {
                __m128d m1, m2, m3, m4, m5, m6, m7,  m8;
                double __attribute__((aligned(16))) a_data;
                a_data= a;
                m8 = _mm_load_pd1(&a_data);

                unsigned long x_address = reinterpret_cast<unsigned long>(x);
                unsigned long x_offset = x_address % 16;

                unsigned long z_offset(x_offset / 8);

                unsigned long quad_start = z_offset;
                unsigned long quad_end(size - ((size - quad_start) % 14));

                if (size < 24)
                {
                    quad_end = 0;
                    quad_start = 0;
                }

                for (unsigned long index(quad_start) ; index < quad_end ; index += 14)
                {
                    m1 = _mm_load_pd(x + index);
                    m2 = _mm_load_pd(x + index + 2);
                    m3 = _mm_load_pd(x + index + 4);
                    m4 = _mm_load_pd(x + index + 6);
                    m5 = _mm_load_pd(x + index + 8);
                    m6 = _mm_load_pd(x + index + 10);
                    m7 = _mm_load_pd(x + index + 12);

                    m1 = _mm_add_pd(m1, m8);
                    m2 = _mm_add_pd(m2, m8);
                    m3 = _mm_add_pd(m3, m8);
                    m4 = _mm_add_pd(m4, m8);
                    m5 = _mm_add_pd(m5, m8);
                    m6 = _mm_add_pd(m6, m8);
                    m7 = _mm_add_pd(m7, m8);

                    _mm_stream_pd(x + index, m1);
                    _mm_stream_pd(x + index + 2, m2);
                    _mm_stream_pd(x + index + 4, m3);
                    _mm_stream_pd(x + index + 6, m4);
                    _mm_stream_pd(x + index + 8, m5);
                    _mm_stream_pd(x + index + 10, m6);
                    _mm_stream_pd(x + index + 12, m7);
                }

                for (unsigned long index(0) ; index < quad_start ; index++)
                {
                    x[index] += a;
                }

                for (unsigned long index = quad_end ; index < size ; index++)
                {
                    x[index] += a;
                }

                _mm_sfence();
            }
        }
    }
}

using namespace honei;

DenseVectorContinuousBase<float> & Sum<tags::CPU::SSE>::value(DenseVectorContinuousBase<float> & a,
        const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When adding DenseVectorContinuousBase<float> to DenseVectorContinuousBase<float> with SSE:");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    intern::sse::sum(a.elements(), b.elements(), a.size());

    return a;
}

DenseVectorContinuousBase<double> & Sum<tags::CPU::SSE>::value(DenseVectorContinuousBase<double> & a,
        const DenseVectorContinuousBase<double> & b)
{
    CONTEXT("When adding DenseVectorContinuousBase<double> to DenseVectorContinuousBase<double> with SSE:");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    intern::sse::sum(a.elements(), b.elements(), a.size());

    return a;
}

DenseMatrix<float> & Sum<tags::CPU::SSE>::value(DenseMatrix<float> & a, const DenseMatrix<float> & b)
{
    CONTEXT("When adding DenseMatrix<float> to DenseMatrix<float> with SSE:");

    if (a.columns() != b.columns())
    {
        throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
    }

    if (a.rows() != b.rows())
    {
        throw MatrixRowsDoNotMatch(b.rows(), a.rows());
    }

    intern::sse::sum(a.elements(), b.elements(), a.rows() * a.columns());

    return a;
}

DenseMatrix<double> & Sum<tags::CPU::SSE>::value(DenseMatrix<double> & a, const DenseMatrix<double> & b)
{
    CONTEXT("When adding DenseMatrix<double> to DenseMatrix<double> with SSE:");

    if (a.columns() != b.columns())
    {
        throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
    }

    if (a.rows() != b.rows())
    {
        throw MatrixRowsDoNotMatch(b.rows(), a.rows());
    }

    intern::sse::sum(a.elements(), b.elements(), a.rows() * a.columns());

    return a;
}

DenseVectorContinuousBase<float> & Sum<tags::CPU::SSE>::value(DenseVectorContinuousBase<float> & x, const float a)
{
    CONTEXT("When adding DenseVectorContinuousBase<float> and float with SSE:");

    intern::sse::sum(a, x.elements(), x.size());

    return x;
}

DenseVectorContinuousBase<double> & Sum<tags::CPU::SSE>::value(DenseVectorContinuousBase<double> & x, const double a)
{
    CONTEXT("When adding DenseVectorContinuousBase<double> and double with SSE:");

    intern::sse::sum(a, x.elements(), x.size());

    return x;
}

DenseMatrix<float> & Sum<tags::CPU::SSE>::value(DenseMatrix<float> & x, const float a)
{
    CONTEXT("When adding DenseMatrix<float> and float with SSE:");

    intern::sse::sum(a, x.elements(), x.rows() * x.columns());

    return x;
}

DenseMatrix<double> & Sum<tags::CPU::SSE>::value(DenseMatrix<double> & x, const double a)
{
    CONTEXT("When adding DenseMatrix<double> and double with SSE:");

    intern::sse::sum(a, x.elements(), x.rows() * x.columns());

    return x;
}

