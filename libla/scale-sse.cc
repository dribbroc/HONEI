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

#include <libla/scale.hh>

#include <xmmintrin.h>
#include <emmintrin.h>

namespace honei
{
    namespace intern
    {
        namespace sse
        {
            inline void scale(const float a, float * x, unsigned long size)
            {
                __m128 m1, m2, m3, m4, m5, m6, m7, m8;
                float __attribute__((aligned(16))) a_data;
                a_data= a;
                m8 = _mm_load_ps1(&a_data);

                unsigned long x_address = (unsigned long)x;
                unsigned long x_offset = x_address % 16;

                unsigned long z_offset(x_offset / 4);
                z_offset = (4 - z_offset) % 4;

                unsigned long quad_start = z_offset;
                unsigned long quad_end(size - ((size - quad_start) % 28));
                if (size < 32)
                {
                    quad_end = 0;
                    quad_start = 0;
                }

                for (unsigned long index = quad_start ; index < quad_end ; index += 28)
                {
                    m1 = _mm_load_ps(x + index);
                    m2 = _mm_load_ps(x + index + 4);
                    m3 = _mm_load_ps(x + index + 8);
                    m4 = _mm_load_ps(x + index + 12);
                    m5 = _mm_load_ps(x + index + 16);
                    m6 = _mm_load_ps(x + index + 20);
                    m7 = _mm_load_ps(x + index + 24);

                    m1 = _mm_mul_ps(m1, m8);
                    m2 = _mm_mul_ps(m2, m8);
                    m3 = _mm_mul_ps(m3, m8);
                    m4 = _mm_mul_ps(m4, m8);
                    m5 = _mm_mul_ps(m5, m8);
                    m6 = _mm_mul_ps(m6, m8);
                    m7 = _mm_mul_ps(m7, m8);

                    _mm_stream_ps(x + index, m1);
                    _mm_stream_ps(x + index + 4, m2);
                    _mm_stream_ps(x + index + 8, m3);
                    _mm_stream_ps(x + index + 12, m4);
                    _mm_stream_ps(x + index + 16, m5);
                    _mm_stream_ps(x + index + 20, m6);
                    _mm_stream_ps(x + index + 24, m7);
                }

                for (unsigned long index = 0 ; index < quad_start ; index++)
                {
                    x[index] *= a;
                }
                for (unsigned long index = quad_end ; index < size ; index++)
                {
                    x[index] *= a;
                }
                _mm_sfence();
            }

            inline void scale(const double a, double * x, unsigned long size)
            {
                __m128d m1, m2, m3, m4, m5, m6, m7,  m8;
                double __attribute__((aligned(16))) a_data;
                a_data= a;
                m8 = _mm_load_pd1(&a_data);

                unsigned long x_address = (unsigned long)x;
                unsigned long x_offset = x_address % 16;

                unsigned long z_offset(x_offset / 8);

                unsigned long quad_start = z_offset;
                unsigned long quad_end(size - ((size - quad_start) % 14));
                if (size < 16)
                {
                    quad_end = 0;
                    quad_start = 0;
                }

                for (unsigned long index = quad_start ; index < quad_end ; index += 14)
                {
                    m1 = _mm_load_pd(x + index);
                    m2 = _mm_load_pd(x + index + 2);
                    m3 = _mm_load_pd(x + index + 4);
                    m4 = _mm_load_pd(x + index + 6);
                    m5 = _mm_load_pd(x + index + 8);
                    m6 = _mm_load_pd(x + index + 10);
                    m7 = _mm_load_pd(x + index + 12);

                    m1 = _mm_mul_pd(m1, m8);
                    m2 = _mm_mul_pd(m2, m8);
                    m3 = _mm_mul_pd(m3, m8);
                    m4 = _mm_mul_pd(m4, m8);
                    m5 = _mm_mul_pd(m5, m8);
                    m6 = _mm_mul_pd(m6, m8);
                    m7 = _mm_mul_pd(m7, m8);

                    _mm_stream_pd(x + index, m1);
                    _mm_stream_pd(x + index + 2, m2);
                    _mm_stream_pd(x + index + 4, m3);
                    _mm_stream_pd(x + index + 6, m4);
                    _mm_stream_pd(x + index + 8, m5);
                    _mm_stream_pd(x + index + 10, m6);
                    _mm_stream_pd(x + index + 12, m7);
                }

                for (unsigned long index = 0 ; index < quad_start ; index++)
                {
                    x[index] *= a;
                }
                for (unsigned long index = quad_end ; index < size ; index++)
                {
                    x[index] *= a;
                }
                _mm_sfence();
            }
        }
    }
}

using namespace honei;

DenseVectorContinuousBase<float> & Scale<tags::CPU::SSE>::value(const float a, DenseVectorContinuousBase<float> & x)
{
    CONTEXT("When scaling DenseVectorContinuousBase<float> by float with SSE:");

    intern::sse::scale(a, x.elements(), x.size());

    return x;
}

DenseVectorContinuousBase<double> & Scale<tags::CPU::SSE>::value(const double a, DenseVectorContinuousBase<double> & x)
{
    CONTEXT("When scaling DenseVectorContinuousBase<double> by double with SSE:");

    intern::sse::scale(a, x.elements(), x.size());

    return x;
}

DenseMatrix<float> & Scale<tags::CPU::SSE>::value(const float a, DenseMatrix<float> & x)
{
    CONTEXT("When scaling DenseMatrix<float> by float with SSE:");

    intern::sse::scale(a, x.elements(), x.rows() * x.columns());

    return x;
}

DenseMatrix<double> & Scale<tags::CPU::SSE>::value(const double a, DenseMatrix<double> & x)
{
    CONTEXT("When scaling DenseMatrix<double> by double with SSE:");

    intern::sse::scale(a, x.elements(), x.rows() * x.columns());

    return x;
}

SparseVector<float> & Scale<tags::CPU::SSE>::value(const float a, SparseVector<float> & x)
{
    CONTEXT("When scaling SparseVector<float> by float with SSE:");

    intern::sse::scale(a, x.elements(), x.used_elements());

    return x;
}

SparseVector<double> & Scale<tags::CPU::SSE>::value(const double a, SparseVector<double> & x)
{
    CONTEXT("When scaling SparseVector<double> by double with SSE:");

    intern::sse::scale(a, x.elements(), x.used_elements());
    return x;
}

SparseMatrix<float> & Scale<tags::CPU::SSE>::value(const float a, SparseMatrix<float> & x)
{
    CONTEXT("When scaling SparseMatrix<float> by float with SSE:");

    for (SparseMatrix<float>::RowIterator l(x.begin_non_zero_rows()),
            l_end(x.end_non_zero_rows()) ; l != l_end ; ++l)
    {
        intern::sse::scale(a, l->elements(), l->used_elements());
    }

    return x;
}

SparseMatrix<double> & Scale<tags::CPU::SSE>::value(const double a, SparseMatrix<double> & x)
{
    CONTEXT("When scaling SparseMatrix<double> by double with SSE:");

    for (SparseMatrix<double>::RowIterator l(x.begin_non_zero_rows()),
            l_end(x.end_non_zero_rows()) ; l != l_end ; ++l)
    {
        intern::sse::scale(a, l->elements(), l->used_elements());
    }

    return x;
}
