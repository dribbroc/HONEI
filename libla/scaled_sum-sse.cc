
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
#include <emmintrin.h>

using namespace honei;

DenseVector<float> & ScaledSum<tags::CPU::SSE>::value(DenseVector<float> & x, const DenseVector<float> & y, float b)
{
    CONTEXT("When calculating ScaledSum form DenseVector<float> with SSE:");
    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());

    __m128 m1, m2, m8;
    float __attribute__((aligned(16))) b_data;
    b_data = b;
    m8 = _mm_load1_ps(&b_data);

    unsigned long quad_end(x.size() - (y.size() % 4));

    for (unsigned long index = 0 ; index < quad_end ; index += 4)
    {

        m1 = _mm_load_ps(x.elements() + index);
        m2 = _mm_load_ps(y.elements() + index);

        m2 = _mm_mul_ps(m2, m8);

        m1 = _mm_add_ps(m1, m2);

        _mm_stream_ps(x.elements() + index, m1);
    }
    for (unsigned long index = quad_end ; index < x.size() ; index++)
    {
        x.elements()[index] += y.elements()[index] * b;
    }
    return x; 
}

DenseVector<double> & ScaledSum<tags::CPU::SSE>::value(DenseVector<double> & x, const DenseVector<double> & y, double b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorBase<double> with SSE:");
    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());

    __m128d m1, m2, m8;
    double __attribute__((aligned(16))) b_data;
    b_data = b;
    m8 = _mm_load1_pd(&b_data);

    unsigned long quad_end(x.size() - (y.size() % 2));

    for (unsigned long index = 0 ; index < quad_end ; index += 2)
    {

        m1 = _mm_load_pd(x.elements() + index);
        m2 = _mm_load_pd(y.elements() + index);

        m2 = _mm_mul_pd(m2, m8);

        m1 = _mm_add_pd(m1, m2);

        _mm_stream_pd(x.elements() + index, m1);
    }
    for (unsigned long index = quad_end ; index < x.size() ; index++)
    {
        x.elements()[index] += y.elements()[index] * b;
    }
    return x; 
}

DenseVectorContinuousBase<float> & ScaledSum<tags::CPU::SSE>::value(DenseVectorContinuousBase<float> & x, const DenseVectorContinuousBase<float> & y, float b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<float> with SSE:");
    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());

    __m128 m1, m2, m8;
    float __attribute__((aligned(16))) b_data;
    b_data = b;
    m8 = _mm_load1_ps(&b_data);

    unsigned long x_address = (unsigned long)x.elements();
    unsigned long x_offset = x_address % 16;

    unsigned long z_offset(x_offset / 4);
    z_offset = (4 - z_offset) % 4;

    unsigned long quad_start = z_offset;
    unsigned long quad_end(x.size() - ((x.size()-quad_start) % 4));
    if (quad_end < 4 || quad_end > x.size())
    {
        quad_end = x.size();
        quad_start = x.size();
    }

    for (unsigned long index = quad_start ; index < quad_end ; index += 4)
    {

        m1 = _mm_load_ps(x.elements() + index);
        m2 = _mm_loadu_ps(y.elements() + index);

        m2 = _mm_mul_ps(m2, m8);

        m1 = _mm_add_ps(m1, m2);

        _mm_stream_ps(x.elements() + index, m1);
    }
    for (unsigned long index(0) ; index < quad_start ; index++)
    {
        x.elements()[index] += y.elements()[index] * b;
    }
    for (unsigned long index = quad_end ; index < x.size() ; index++)
    {
        x.elements()[index] += y.elements()[index] * b;
    }
    return x; 
}

DenseVectorContinuousBase<double> & ScaledSum<tags::CPU::SSE>::value(DenseVectorContinuousBase<double> & x, const DenseVectorContinuousBase<double> & y, double b)
{
    CONTEXT("When calculating ScaledSum form DenseVectoriContinuousBase<double> with SSE:");
    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());

    __m128d m1, m2, m8;
    double __attribute__((aligned(16))) b_data;
    b_data = b;
    m8 = _mm_load1_pd(&b_data);

    unsigned long x_address = (unsigned long)x.elements();
    unsigned long x_offset = x_address % 16;

    unsigned long z_offset(x_offset / 8);

    unsigned long quad_start = z_offset;
    unsigned long quad_end(x.size() - ((x.size()-quad_start) % 2));

    for (unsigned long index = quad_start ; index < quad_end ; index += 2)
    {

        m1 = _mm_load_pd(x.elements() + index);
        m2 = _mm_loadu_pd(y.elements() + index);

        m2 = _mm_mul_pd(m2, m8);

        m1 = _mm_add_pd(m1, m2);

        _mm_stream_pd(x.elements() + index, m1);
    }
    for (unsigned long index(0) ; index < quad_start ; index++)
    {
        x.elements()[index] += y.elements()[index] * b;
    }
    for (unsigned long index = quad_end ; index < x.size() ; index++)
    {
        x.elements()[index] += y.elements()[index] * b;
    }
    return x;
}
