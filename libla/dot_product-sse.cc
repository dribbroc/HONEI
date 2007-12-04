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

#include <libla/dot_product.hh>

#include <xmmintrin.h>
#include <emmintrin.h>

namespace honei
{
    namespace intern
    {
        namespace sse
        {
            float dot_product(const float * a, float * b, unsigned long size)
            {
                float __attribute__((aligned(16))) result(0);
                __m128 m1, m2, m8, m4, m5;

                unsigned long a_address = (unsigned long)a;
                unsigned long a_offset = a_address % 16;
                unsigned long b_address = (unsigned long)b;
                unsigned long b_offset = b_address % 16;

                unsigned long x_offset(a_offset / 4);
                x_offset = (4 - x_offset) % 4;

                unsigned long quad_start = x_offset;
                unsigned long quad_end(size - ((size - quad_start) % 4));
                if (size < 16)
                {
                    quad_end = 0;
                    quad_start = 0;
                }

                m8 = _mm_setzero_ps();
                m4 = _mm_setzero_ps();
                m5 = _mm_setzero_ps();

                /// \todo unroll loops
                if (a_offset == b_offset)
                {
                    for (unsigned long index(quad_start) ; index < quad_end ; index += 4) 
                    {
                        m1 = _mm_load_ps(a + index);
                        m2 = _mm_load_ps(b + index);

                        m1 = _mm_mul_ps(m1, m2);
                        m8 = _mm_add_ps(m1, m8);
                    }
                }
                else
                {
                    for (unsigned long index(quad_start) ; index < quad_end ; index += 4) 
                    {
                        m1 = _mm_load_ps(a + index);
                        m2 = _mm_loadu_ps(b + index);

                        m1 = _mm_mul_ps(m1, m2);
                        m8 = _mm_add_ps(m1, m8);
                    }
                }

                //move the 4 m8 floats to the lower m4 field and sum them up in m5
                m5 = _mm_add_ss(m5, m8);
                m4 = _mm_shuffle_ps(m8, m8, 0x1);
                m5 = _mm_add_ss(m5, m4);
                m4 = _mm_shuffle_ps(m8, m8, 0x2);
                m5 = _mm_add_ss(m5, m4);
                m4 = _mm_shuffle_ps(m8, m8, 0x3);
                m5 = _mm_add_ss(m5, m4);
                _mm_store_ss(&result, m5);

                for (unsigned long index(0) ; index < quad_start ; index++)
                {
                    result += a[index] * b[index];
                }
                for (unsigned long index(quad_end) ; index < size ; index++)
                {
                    result += a[index] * b[index];
                }
                return result;
            }

            double dot_product(double * a, double * b, unsigned long size)
            {
                double __attribute__((aligned(16))) result(0);
                union sse2
                {
                    __m128d m;
                    double d[2];
                } m1, m2, m8;

                unsigned long a_address = (unsigned long)a;
                unsigned long a_offset = a_address % 16;
                unsigned long b_address = (unsigned long)b;
                unsigned long b_offset = b_address % 16;

                unsigned long x_offset(a_offset / 8);

                unsigned long quad_start = x_offset;
                unsigned long quad_end(size - ((size - quad_start) % 2));
                if (size < 16)
                {
                    quad_start = 0;
                    quad_end = 0;
                }
                m8.m = _mm_setzero_pd();

                /// \todo unroll loops
                if(a_offset == b_offset)
                {
                    for (unsigned long index = quad_start ; index < quad_end ; index += 2) 
                    {
                        m1.m = _mm_load_pd(a + index);
                        m2.m = _mm_load_pd(b + index);

                        m1.m = _mm_mul_pd(m1.m, m2.m);
                        m8.m = _mm_add_pd(m1.m, m8.m);
                    }
                }
                else
                {
                    for (unsigned long index = quad_start ; index < quad_end ; index += 2) 
                    {
                        m1.m = _mm_load_pd(a + index);
                        m2.m = _mm_loadu_pd(b + index);

                        m1.m = _mm_mul_pd(m1.m, m2.m);
                        m8.m = _mm_add_pd(m1.m, m8.m);
                    }
                }
                result += m8.d[0];
                result += m8.d[1];

                for (unsigned long index(0) ; index < quad_start ; index++)
                {
                    result += a[index] * b[index];
                }
                for (unsigned long index = quad_end ; index < size ; index++)
                {
                    result += a[index] * b[index];
                }

                return result;
            }
        }
    }
}

using namespace honei;


float DotProduct<tags::CPU::SSE>::value(const DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When calculating dot-product of DenseVectorContinuousBase<float> with DenseVectorContinuousBase<float> with SSE:");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    return intern::sse::dot_product(a.elements(), b.elements(), a.size());
}

double DotProduct<tags::CPU::SSE>::value(const DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b)
{
    CONTEXT("When calculating dot-product of DenseVectorContinuousBase<double> with DenseVectorContinuousBase<double> with SSE:");
    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    return intern::sse::dot_product(a.elements(), b.elements(), a.size());
}

