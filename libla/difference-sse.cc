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

#include <libla/difference.hh>

#include <xmmintrin.h>
#include <emmintrin.h>

using namespace honei;

DenseVector<float> & Difference<tags::CPU::SSE>::value(DenseVector<float> & a, const DenseVector<float> & b)
{
    CONTEXT("When subtracting DenseVector<float> to DenseVector<float> with SSE:");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());


    __m128 m1, m2;

    unsigned long quad_end(a.size() - (a.size() % 4));

    for (unsigned long index(0) ; index < quad_end ; index += 4) 
    {
        m1 = _mm_load_ps(a.elements() + index);
        m2 = _mm_load_ps(b.elements() + index);

        m1 = _mm_sub_ps(m1, m2);

        _mm_stream_ps(a.elements() + index, m1);
    }

    for (unsigned long index = quad_end ; index < a.size() ; index++)
    {
        a.elements()[index] -= b.elements()[index];
    }
    return a;
}

DenseVector<double> & Difference<tags::CPU::SSE>::value(DenseVector<double> & a, const DenseVector<double> & b)
{
    CONTEXT("When subtacting DenseVector<double> to DenseVector<double> with SSE:");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());


    __m128d m1, m2;

    unsigned long quad_end(a.size() - (a.size() % 2));

    for (unsigned long index(0) ; index < quad_end ; index += 2) 
    {
        m1 = _mm_load_pd(a.elements() + index);
        m2 = _mm_load_pd(b.elements() + index);

        m1 = _mm_sub_pd(m1, m2);

        _mm_stream_pd(a.elements() + index, m1);
    }

    for (unsigned long index = quad_end ; index < a.size() ; index++)
    {
        a.elements()[index] -= b.elements()[index];
    }
    return a;
}

DenseVectorContinuousBase<float> & Difference<tags::CPU::SSE>::value(DenseVectorContinuousBase<float> & a, const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When subtracting DenseVectorContinuousBase<float> to DenseVectorContinuousBase<float> with SSE:");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    __m128 m1, m2;

    unsigned long a_address = (unsigned long)a.elements();
    unsigned long a_offset = a_address % 8;

    unsigned long x_offset((4-a_offset)%4);

    unsigned long quad_start = x_offset;
    unsigned long quad_end(a.size() - ((a.size()-quad_start) % 4));

    for (unsigned long index = quad_start ; index < quad_end ; index += 4) 
    {
        m1 = _mm_load_ps(a.elements() + index);
        m2 = _mm_loadu_ps(b.elements() + index);

        m1 = _mm_sub_ps(m1, m2);

        _mm_stream_ps(a.elements() + index, m1);
    }

    for (unsigned long index = 0 ; index < quad_start ; index++)
    {
        a.elements()[index] -= b.elements()[index];
    }
    for (unsigned long index = quad_end ; index < a.size() ; index++)
    {
        a.elements()[index] -= b.elements()[index];
    }
    return a;
}

DenseVectorContinuousBase<double> & Difference<tags::CPU::SSE>::value(DenseVectorContinuousBase<double> & a, const DenseVectorContinuousBase<double> & b)
{
    CONTEXT("When subtacting DenseVectorContinuousBase<double> to DenseVectorContinuousBase<double> with SSE:");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    __m128d m1,m2;

    unsigned long a_address = (unsigned long)a.elements();
    unsigned long a_offset = a_address % 16;

    unsigned long x_offset((2 - a_offset) % 2);

    unsigned long quad_start = x_offset;
    unsigned long quad_end(a.size() - ((a.size() - quad_start) % 2));

    for (unsigned long index = quad_start ; index < quad_end ; index += 2) 
    {
        m1 = _mm_load_pd(a.elements() + index);
        m2 = _mm_loadu_pd(b.elements() + index);

        m1 = _mm_sub_pd(m1, m2);

        _mm_stream_pd(a.elements() + index, m1);
    }

    for (unsigned long index = 0 ; index < quad_start ; index++)
    {
        a.elements()[index] -= b.elements()[index];
    }
    for (unsigned long index = quad_end ; index < a.size() ; index++)
    {
        a.elements()[index] -= b.elements()[index];
    }
    return a;
}
