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

using namespace honei;

DenseVector<float> & Sum<tags::CPU::SSE>::value(DenseVector<float> & a, const DenseVector<float> & b)
{
    CONTEXT("When adding DenseVector to DenseVector with SSE:");


    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());


    __m128 m1, m2, m3, m4, m5, m6, m7, m8;
    float __attribute__((aligned(16))) a1_data[4];
    float __attribute__((aligned(16))) b1_data[4];
    float __attribute__((aligned(16))) a2_data[4];
    float __attribute__((aligned(16))) b2_data[4];
    float __attribute__((aligned(16))) a3_data[4];
    float __attribute__((aligned(16))) b3_data[4];
    float __attribute__((aligned(16))) a4_data[4];
    float __attribute__((aligned(16))) b4_data[4];

    unsigned long quad_end(a.size() - (b.size() % 16));

    for (unsigned long index = 0 ; index < quad_end ; index += 16) 
    {
        for (int i = 0 ; i < 4 ; ++i)
        {
            a1_data[i] = (a.elements())[index + i];
            a2_data[i] = (a.elements())[index + i + 4];
            a3_data[i] = (a.elements())[index + i + 8];
            a4_data[i] = (a.elements())[index + i + 12];
            b1_data[i] = (b.elements())[index + i];
            b2_data[i] = (b.elements())[index + i + 4];
            b3_data[i] = (b.elements())[index + i + 8];
            b4_data[i] = (b.elements())[index + i + 12];
        }
        m1 = _mm_load_ps(a1_data);
        m2 = _mm_load_ps(b1_data);
        m3 = _mm_load_ps(a2_data);
        m4 = _mm_load_ps(b2_data);
        m5 = _mm_load_ps(a3_data);
        m6 = _mm_load_ps(b3_data);
        m7 = _mm_load_ps(a4_data);
        m8 = _mm_load_ps(b4_data);

        m1 = _mm_add_ps(m1, m2);
        m3 = _mm_add_ps(m3, m4);
        m5 = _mm_add_ps(m5, m6);
        m7 = _mm_add_ps(m7, m8);

        _mm_store_ps(a1_data, m1);
        _mm_store_ps(a2_data, m3);
        _mm_store_ps(a3_data, m5);
        _mm_store_ps(a4_data, m7);

        for (int i = 0 ; i < 4 ; ++i)
        {
            (a.elements())[index + i] = a1_data[i];
            (a.elements())[index + i + 4] = a2_data[i];
            (a.elements())[index + i + 8] = a3_data[i];
            (a.elements())[index + i + 12] = a4_data[i];
        }
    }

    for (unsigned long index = quad_end ; index < a.size() ; index++)
    {
        a.elements()[index] += b.elements()[index];
    }
    return a;
}

DenseVector<double> & Sum<tags::CPU::SSE>::value(DenseVector<double> & a, const DenseVector<double> & b)
{
    CONTEXT("When adding DenseVector to DenseVector with SSE:");


    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());


    __m128 m1, m2, m3, m4, m5, m6, m7, m8;
    float __attribute__((aligned(16))) a1_data[2];
    float __attribute__((aligned(16))) b1_data[2];
    float __attribute__((aligned(16))) a2_data[2];
    float __attribute__((aligned(16))) b2_data[2];
    float __attribute__((aligned(16))) a3_data[2];
    float __attribute__((aligned(16))) b3_data[2];
    float __attribute__((aligned(16))) a4_data[2];
    float __attribute__((aligned(16))) b4_data[2];

    unsigned long quad_end(a.size() - (b.size() % 8));

    for (unsigned long index = 0 ; index < quad_end ; index += 8) 
    {
        for (int i = 0 ; i < 2 ; ++i)
        {
            a1_data[i] = (a.elements())[index + i];
            a2_data[i] = (a.elements())[index + i + 2];
            a3_data[i] = (a.elements())[index + i + 4];
            a4_data[i] = (a.elements())[index + i + 6];
            b1_data[i] = (b.elements())[index + i];
            b2_data[i] = (b.elements())[index + i + 2];
            b3_data[i] = (b.elements())[index + i + 4];
            b4_data[i] = (b.elements())[index + i + 6];
        }
        m1 = _mm_load_ps(a1_data);
        m2 = _mm_load_ps(b1_data);
        m3 = _mm_load_ps(a2_data);
        m4 = _mm_load_ps(b2_data);
        m5 = _mm_load_ps(a3_data);
        m6 = _mm_load_ps(b3_data);
        m7 = _mm_load_ps(a4_data);
        m8 = _mm_load_ps(b4_data);

        m1 = _mm_add_ps(m1, m2);
        m3 = _mm_add_ps(m3, m4);
        m5 = _mm_add_ps(m5, m6);
        m7 = _mm_add_ps(m7, m8);

        _mm_store_ps(a1_data, m1);
        _mm_store_ps(a2_data, m3);
        _mm_store_ps(a3_data, m5);
        _mm_store_ps(a4_data, m7);

        for (int i = 0 ; i < 2 ; ++i)
        {
            (a.elements())[index + i] = a1_data[i];
            (a.elements())[index + i + 2] = a2_data[i];
            (a.elements())[index + i + 4] = a3_data[i];
            (a.elements())[index + i + 6] = a4_data[i];
        }
    }

    for (unsigned long index = quad_end ; index < a.size() ; index++)
    {
        a.elements()[index] += b.elements()[index];
    }
    return a;
}
