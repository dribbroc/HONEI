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
    CONTEXT("When scaling DenseVector<float> by float with SSE:");

    __m128 m1, m8;
    float __attribute__((aligned(16))) a_data;
    a_data= a;
    m8 = _mm_load_ps1(&a_data);

    for (unsigned long index = 0 ; index < x.size() ; index += 4) 
    {
        m1 = _mm_load_ps(x.elements() + index);

        m1 = _mm_mul_ps(m1, m8);

        _mm_stream_ps(x.elements() + index, m1);
    }

    return x;
}

DenseVector<double> & Scale<tags::CPU::SSE>::value(const double a, DenseVector<double> & x)
{
    CONTEXT("When scaling DenseVector<double> by double with SSE:");

    __m128d m1, m8;
    double __attribute__((aligned(16))) a_data;
    a_data= a;
    m8 = _mm_load_pd1(&a_data);

    for (unsigned long index = 0 ; index < x.size() ; index += 2) 
    {
        m1 = _mm_load_pd(x.elements() + index);

        m1 = _mm_mul_pd(m1, m8);

        _mm_stream_pd(x.elements() + index, m1);
    }

    return x;
}

DenseMatrix<float> & Scale<tags::CPU::SSE>::value(const float a, DenseMatrix<float> & x)
{
    CONTEXT("When scaling DenseMatrix<float> by float with SSE:");

    __m128 m1, m8;
    float __attribute__((aligned(16))) a_data;
    a_data= a;
    m8 = _mm_load_ps1(&a_data);

    for (unsigned long index = 0 ; index < x.rows() * x.columns() ; index += 4) 
    {
        m1 = _mm_load_ps(x.elements() + index);

        m1 = _mm_mul_ps(m1, m8);

        _mm_stream_ps(x.elements() + index, m1);
    }

    return x;
}

DenseMatrix<double> & Scale<tags::CPU::SSE>::value(const double a, DenseMatrix<double> & x)
{
    CONTEXT("When scaling DenseMatrix<double> by double with SSE:");

    __m128d m1, m8;
    double __attribute__((aligned(16))) a_data;
    a_data= a;
    m8 = _mm_load_pd1(&a_data);

    for (unsigned long index = 0 ; index < x.rows() * x.columns() ; index += 2) 
    {
        m1 = _mm_load_pd(x.elements() + index);

        m1 = _mm_mul_pd(m1, m8);

        _mm_stream_pd(x.elements() + index, m1);
    }

    return x;
}

DenseVectorContinuousBase<float> & Scale<tags::CPU::SSE>::value(const float a, DenseVectorContinuousBase<float> & x)
{
    CONTEXT("When scaling DenseVectorContinuousBase<float> by float with SSE:");

    __m128 m1, m8;
    float __attribute__((aligned(16))) a_data;
    a_data= a;
    m8 = _mm_load_ps1(&a_data);

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

        m1 = _mm_mul_ps(m1, m8);

        _mm_stream_ps(x.elements() + index, m1);
    }

    for (unsigned long index = 0 ; index < quad_start ; index++)
    {
        x.elements()[index] *= a;
    }
    for (unsigned long index = quad_end ; index < x.size() ; index++)
    {
        x.elements()[index] *= a;
    }
    return x;
}

DenseVectorContinuousBase<double> & Scale<tags::CPU::SSE>::value(const double a, DenseVectorContinuousBase<double> & x)
{
    CONTEXT("When scaling DenseVectorContinuousBase<double> by double with SSE:");

    __m128d m1, m8;
    double __attribute__((aligned(16))) a_data;
    a_data= a;
    m8 = _mm_load_pd1(&a_data);

    unsigned long x_address = (unsigned long)x.elements();
    unsigned long x_offset = x_address % 16;

    unsigned long z_offset(x_offset / 8);

    unsigned long quad_start = z_offset;
    unsigned long quad_end(x.size() - ((x.size()-quad_start) % 2));

    for (unsigned long index = quad_start ; index < quad_end ; index += 2) 
    {
        m1 = _mm_load_pd(x.elements() + index);

        m1 = _mm_mul_pd(m1, m8);

        _mm_stream_pd(x.elements() + index, m1);
    }

    for (unsigned long index = 0 ; index < quad_start ; index++)
    {
        x.elements()[index] *= a;
    }
    for (unsigned long index = quad_end ; index < x.size() ; index++)
    {
        x.elements()[index] *= a;
    }
    return x;
}
