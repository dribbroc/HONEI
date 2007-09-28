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

using namespace honei;

DenseVector<float> & Scale<tags::CPU::SSE>::value(const float a, DenseVector<float> & x)
{
    CONTEXT("When scaling DenseVector by float with SSE:");

    __m128 m1, m2, m3, m4, m5, m6, m7, m8;
    float __attribute__((aligned(16))) x1_data[4];
    float __attribute__((aligned(16))) x2_data[4];
    float __attribute__((aligned(16))) x3_data[4];
    float __attribute__((aligned(16))) x4_data[4];
    float __attribute__((aligned(16))) x5_data[4];
    float __attribute__((aligned(16))) x6_data[4];
    float __attribute__((aligned(16))) x7_data[4];
    float __attribute__((aligned(16))) a_data;
    a_data= a;
    m8 = _mm_load_ps1(&a_data);

    unsigned long quad_end(x.size() - (x.size() % 28));
    for (unsigned long index = 0 ; index < quad_end ; index += 28) 
    {
        for (int i = 0 ; i < 4 ; ++i)
        {
            x1_data[i] = (x.elements())[index + i];
            x2_data[i] = (x.elements())[index + i + 4];
            x3_data[i] = (x.elements())[index + i + 8];
            x4_data[i] = (x.elements())[index + i + 12];
            x5_data[i] = (x.elements())[index + i + 16];
            x6_data[i] = (x.elements())[index + i + 20];
            x7_data[i] = (x.elements())[index + i + 24];
        }
        m1 = _mm_load_ps(x1_data);
        m2 = _mm_load_ps(x2_data);
        m3 = _mm_load_ps(x3_data);
        m4 = _mm_load_ps(x4_data);
        m5 = _mm_load_ps(x5_data);
        m6 = _mm_load_ps(x6_data);
        m7 = _mm_load_ps(x7_data);

        m1 = _mm_mul_ps(m1, m8);
        m2 = _mm_mul_ps(m2, m8);
        m3 = _mm_mul_ps(m3, m8);
        m4 = _mm_mul_ps(m4, m8);
        m5 = _mm_mul_ps(m5, m8);
        m6 = _mm_mul_ps(m6, m8);
        m7 = _mm_mul_ps(m7, m8);

        _mm_stream_ps(x1_data, m1);
        _mm_stream_ps(x2_data, m2);
        _mm_stream_ps(x3_data, m3);
        _mm_stream_ps(x4_data, m4);
        _mm_stream_ps(x5_data, m5);
        _mm_stream_ps(x6_data, m6);
        _mm_stream_ps(x7_data, m7);

        for (int i = 0 ; i < 4 ; ++i)
        {
            (x.elements())[index + i] = x1_data[i];
            (x.elements())[index + i + 4] = x2_data[i];
            (x.elements())[index + i + 8] = x3_data[i];
            (x.elements())[index + i + 12] = x4_data[i];
            (x.elements())[index + i + 16] = x5_data[i];
            (x.elements())[index + i + 20] = x6_data[i];
            (x.elements())[index + i + 24] = x7_data[i];
        }
    }

    for (unsigned long index = quad_end ; index < x.size() ; index++)
    {
        x.elements()[index] *= a;
    }
    return x;
}

DenseVector<double> & Scale<tags::CPU::SSE>::value(const double a, DenseVector<double> & x)
{
    CONTEXT("When scaling DenseVector by double with SSE:");

    __m128d m1, m2, m3, m4, m5, m6, m7, m8;
    double __attribute__((aligned(16))) x1_data[2];
    double __attribute__((aligned(16))) x2_data[2];
    double __attribute__((aligned(16))) x3_data[2];
    double __attribute__((aligned(16))) x4_data[2];
    double __attribute__((aligned(16))) x5_data[2];
    double __attribute__((aligned(16))) x6_data[2];
    double __attribute__((aligned(16))) x7_data[2];
    double __attribute__((aligned(16))) a_data;
    a_data= a;
    m8 = _mm_load_pd1(&a_data);

    unsigned long quad_end(x.size() - (x.size() % 14));
    for (unsigned long index = 0 ; index < quad_end ; index += 14) 
    {
        for (int i = 0 ; i < 2 ; ++i)
        {
            x1_data[i] = (x.elements())[index + i];
            x2_data[i] = (x.elements())[index + i + 2];
            x3_data[i] = (x.elements())[index + i + 4];
            x4_data[i] = (x.elements())[index + i + 6];
            x5_data[i] = (x.elements())[index + i + 8];
            x6_data[i] = (x.elements())[index + i + 10];
            x7_data[i] = (x.elements())[index + i + 12];
        }
        m1 = _mm_load_pd(x1_data);
        m2 = _mm_load_pd(x2_data);
        m3 = _mm_load_pd(x3_data);
        m4 = _mm_load_pd(x4_data);
        m5 = _mm_load_pd(x5_data);
        m6 = _mm_load_pd(x6_data);
        m7 = _mm_load_pd(x7_data);

        m1 = _mm_mul_pd(m1, m8);
        m2 = _mm_mul_pd(m2, m8);
        m3 = _mm_mul_pd(m3, m8);
        m4 = _mm_mul_pd(m4, m8);
        m5 = _mm_mul_pd(m5, m8);
        m6 = _mm_mul_pd(m6, m8);
        m7 = _mm_mul_pd(m7, m8);

        _mm_stream_pd(x1_data, m1);
        _mm_stream_pd(x2_data, m2);
        _mm_stream_pd(x3_data, m3);
        _mm_stream_pd(x4_data, m4);
        _mm_stream_pd(x5_data, m5);
        _mm_stream_pd(x6_data, m6);
        _mm_stream_pd(x7_data, m7);

        for (int i = 0 ; i < 2 ; ++i)
        {
            (x.elements())[index + i] = x1_data[i];
            (x.elements())[index + i + 2] = x2_data[i];
            (x.elements())[index + i + 4] = x3_data[i];
            (x.elements())[index + i + 6] = x4_data[i];
            (x.elements())[index + i + 8] = x5_data[i];
            (x.elements())[index + i + 10] = x6_data[i];
            (x.elements())[index + i + 12] = x7_data[i];
        }
    }

    for (unsigned long index = quad_end ; index < x.size() ; index++)
    {
        x.elements()[index] *= a;
    }
    return x;
}
