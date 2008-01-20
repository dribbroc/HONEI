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

#include <honei/libla/reduction.hh>

#include <xmmintrin.h>
#include <emmintrin.h>

namespace honei
{
    namespace intern
    {
        namespace sse
        {
            inline float reduction_sum(const float * a, unsigned long size)
            {
                float __attribute__((aligned(16))) result(0);
                union sse4
                {
                    __m128 m;
                    float f[4];
                } m1, m2, m3, m4, m5, m6, m7, m8;

                unsigned long a_address = reinterpret_cast<unsigned long>(a);
                unsigned long a_offset = a_address % 16;

                unsigned long x_offset(a_offset / 4);
                x_offset = (4 - x_offset) % 4;

                unsigned long quad_start = x_offset;
                unsigned long quad_end(size - ((size - quad_start) % 28));

                if (size < 32)
                {
                    quad_end = 0;
                    quad_start = 0;
                }

                m8.m = _mm_setzero_ps();

                for (unsigned long index(quad_start) ; index < quad_end ; index += 28)
                {
                    m1.m = _mm_load_ps(a + index);
                    m2.m = _mm_load_ps(a + index + 4);
                    m3.m = _mm_load_ps(a + index + 8);
                    m4.m = _mm_load_ps(a + index + 12);
                    m5.m = _mm_load_ps(a + index + 16);
                    m6.m = _mm_load_ps(a + index + 20);
                    m7.m = _mm_load_ps(a + index + 24);

                    m8.m = _mm_add_ps(m1.m, m8.m);
                    m8.m = _mm_add_ps(m2.m, m8.m);
                    m8.m = _mm_add_ps(m3.m, m8.m);
                    m8.m = _mm_add_ps(m4.m, m8.m);
                    m8.m = _mm_add_ps(m5.m, m8.m);
                    m8.m = _mm_add_ps(m6.m, m8.m);
                    m8.m = _mm_add_ps(m7.m, m8.m);
                }

                result += m8.f[0];
                result += m8.f[1];
                result += m8.f[2];
                result += m8.f[3];

                for (unsigned long index(0) ; index < quad_start ; index++)
                {
                    result += a[index];
                }

                for (unsigned long index(quad_end) ; index < size ; index++)
                {
                    result += a[index];
                }

                _mm_sfence();

                return result;
            }

            inline double reduction_sum(double * a, unsigned long size)
            {
                double __attribute__((aligned(16))) result(0);
                union sse2
                {
                    __m128d m;
                    double d[2];
                } m1, m2, m3, m4, m5, m6, m7, m8;

                unsigned long a_address = (unsigned long)a;
                unsigned long a_offset = a_address % 16;

                unsigned long x_offset(a_offset / 8);

                unsigned long quad_start = x_offset;
                unsigned long quad_end(size - ((size - quad_start) % 14));

                if (size < 20)
                {
                    quad_start = 0;
                    quad_end = 0;
                }

                m8.m = _mm_setzero_pd();

                for (unsigned long index(quad_start) ; index < quad_end ; index += 14)
                {
                    m1.m = _mm_load_pd(a + index);
                    m2.m = _mm_load_pd(a + index + 2);
                    m3.m = _mm_load_pd(a + index + 4);
                    m4.m = _mm_load_pd(a + index + 6);
                    m5.m = _mm_load_pd(a + index + 8);
                    m6.m = _mm_load_pd(a + index + 10);
                    m7.m = _mm_load_pd(a + index + 12);

                    m8.m = _mm_add_pd(m1.m, m8.m);
                    m8.m = _mm_add_pd(m2.m, m8.m);
                    m8.m = _mm_add_pd(m3.m, m8.m);
                    m8.m = _mm_add_pd(m4.m, m8.m);
                    m8.m = _mm_add_pd(m5.m, m8.m);
                    m8.m = _mm_add_pd(m6.m, m8.m);
                    m8.m = _mm_add_pd(m7.m, m8.m);
                }

                result += m8.d[0];
                result += m8.d[1];

                for (unsigned long index(0) ; index < quad_start ; index++)
                {
                    result += a[index];
                }

                for (unsigned long index = quad_end ; index < size ; index++)
                {
                    result += a[index];
                }

                _mm_sfence();

                return result;
            }
        }
    }
}

using namespace honei;

float Reduction<rt_sum, tags::CPU::SSE>::value(const DenseVectorContinuousBase<float> & a)
{
    CONTEXT("When reducing DenseVectorContinuousBase<float> to sum with SSE:");

    return intern::sse::reduction_sum(a.elements(), a.size());
}

double Reduction<rt_sum, tags::CPU::SSE>::value(const DenseVectorContinuousBase<double> & a)
{
    CONTEXT("When reducing DenseVectorContinuousBase<double> to sum with SSE:");

    return intern::sse::reduction_sum(a.elements(), a.size());
}

DenseVector<float> Reduction<rt_sum, tags::CPU::SSE>::value(const DenseMatrix<float> & a)
{
    CONTEXT("When reducing DenseMatrix<float> to sum with SSE:");

    DenseVector<float> result(a.rows());

    for (unsigned long i(0) ; i < a.rows() ; ++i)
    {
        result[i] = intern::sse::reduction_sum(a[i].elements(), a[i].size());
    }

    return result;
}

DenseVector<double> Reduction<rt_sum, tags::CPU::SSE>::value(const DenseMatrix<double> & a)
{
    CONTEXT("When reducing DenseMatrix<double> to sum with SSE:");

    DenseVector<double> result(a.rows());

    for (unsigned long i(0) ; i < a.rows() ; ++i)
    {
        result[i] = intern::sse::reduction_sum(a[i].elements(), a[i].size());
    }

    return result;
}

float Reduction<rt_sum, tags::CPU::SSE>::value(const SparseVector<float> & a)
{
    CONTEXT("When reducing SparseVector<float> to sum with SSE:");

    return intern::sse::reduction_sum(a.elements(), a.used_elements());
}

double Reduction<rt_sum, tags::CPU::SSE>::value(const SparseVector<double> & a)
{
    CONTEXT("When reducing SparseVector<double> to sum with SSE:");

    return intern::sse::reduction_sum(a.elements(), a.used_elements());
}

DenseVector<float> Reduction<rt_sum, tags::CPU::SSE>::value(const SparseMatrix<float> & a)
{
    CONTEXT("When reducing SparseMatrix<float> to sum with SSE:");

    DenseVector<float> result(a.rows());

    for (unsigned long i(0) ; i < a.rows() ; ++i)
    {
        result[i] = intern::sse::reduction_sum(a[i].elements(), a[i].used_elements());
    }

    return result;
}

DenseVector<double> Reduction<rt_sum, tags::CPU::SSE>::value(const SparseMatrix<double> & a)
{
    CONTEXT("When reducing SparseMatrix<double> to sum with SSE:");

    DenseVector<double> result(a.rows());

    for (unsigned long i(0) ; i < a.rows() ; ++i)
    {
        result[i] = intern::sse::reduction_sum(a[i].elements(), a[i].used_elements());
    }

    return result;
}

