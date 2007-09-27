
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

#include <libla/scaled_sum.hh>

#include <xmmintrin.h>

using namespace honei;

DenseVector<float> & ScaledSum<tags::CPU::SSE>::value(DenseVector<float> & x, const DenseVector<float> & y, float b)
{
    CONTEXT("When calculating ScaledSum form DenseVector<float> with SSE:");
    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());

    __m128 m1, m2, m3, m4, m5, m6, m7;
    float __attribute__((aligned(16))) x1_data[4];
    float __attribute__((aligned(16))) y1_data[4];
    float __attribute__((aligned(16))) x2_data[4];
    float __attribute__((aligned(16))) y2_data[4];
    float __attribute__((aligned(16))) x3_data[4];
    float __attribute__((aligned(16))) y3_data[4];
    float __attribute__((aligned(16))) b_data;
    b_data = b;

    unsigned long quad_end(x.size() - (y.size() % 12));

    for (unsigned long index = 0 ; index < quad_end ; index += 12)
    {
        for (int i = 0 ; i < 4 ; ++i)
        {
            x1_data[i] = x.elements()[index + i];
            x2_data[i] = x.elements()[index + i + 4];
            x3_data[i] = x.elements()[index + i + 8];
            y1_data[i] = y.elements()[index + i];
            y2_data[i] = y.elements()[index + i + 4];
            y3_data[i] = y.elements()[index + i + 8];
        }

        m1 = _mm_load_ps(x1_data);
        m2 = _mm_load_ps(y1_data);
        m3 = _mm_load_ps(x2_data);
        m4 = _mm_load_ps(y2_data);
        m5 = _mm_load_ps(x3_data);
        m6 = _mm_load_ps(y3_data);
        m7 = _mm_load1_ps(&b_data);

        m2 = _mm_mul_ps(m2, m7);
        m4 = _mm_mul_ps(m4, m7);
        m6 = _mm_mul_ps(m6, m7);

        m1 = _mm_add_ps(m1, m2);
        m3 = _mm_add_ps(m3, m4);
        m5 = _mm_add_ps(m5, m6);

        _mm_store_ps(x1_data, m1);
        _mm_store_ps(x2_data, m3);
        _mm_store_ps(x3_data, m5);
        for (int i = 0 ; i < 4 ; ++i)

        {
            x.elements()[index + i] = x1_data[i];
            x.elements()[index + i + 4] = x2_data[i];
            x.elements()[index + i + 8] = x3_data[i];
        }
    }
    for (unsigned long index = quad_end ; index < x.size() ; index++)
    {
        x.elements()[index] += y.elements()[index] * b;
    }
    return x; 
}
