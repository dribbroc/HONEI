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

using namespace honei;

float DotProduct<tags::CPU::SSE>::value(const DenseVector<float> & a, const DenseVector<float> & b)
{
    CONTEXT("When calculating dot-product of DenseVector<float> with DenseVector<float> with SSE:");


    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    float result(0);
    union sse4
    {
        __m128 m;
        float f[4];
    } m1, m2, m8;

    unsigned long quad_end(a.size() - (b.size() % 4));

    for (unsigned long index = 0 ; index < quad_end ; index += 4) 
    {
        m1.m = _mm_load_ps(&a.elements()[index]);
        m2.m = _mm_load_ps(&b.elements()[index]);

        m1.m = _mm_mul_ps(m1.m, m2.m);
        m8.m = _mm_add_ps(m1.m, m8.m);
    }
    result = m8.f[0] + m8.f[1] + m8.f[2] + m8.f[3];

    for (unsigned long index = quad_end ; index < a.size() ; index++)
    {
        result += a.elements()[index] * b.elements()[index];
    }
    return result;
}

double DotProduct<tags::CPU::SSE>::value(const DenseVector<double> & a, const DenseVector<double> & b)
{
    CONTEXT("When calculating dot-product of DenseVector<double> with DenseVector<double> with SSE:");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    double result(0);
    union sse2
    {
        __m128d m;
        double d[2];
    } m1, m2, m8;

    unsigned long quad_end(a.size() - (b.size() % 2));

    for (unsigned long index = 0 ; index < quad_end ; index += 2) 
    {
        m1.m = _mm_load_pd(&a.elements()[index]);
        m2.m = _mm_load_pd(&b.elements()[index]);

        m1.m = _mm_mul_pd(m1.m, m2.m);
        m8.m = _mm_add_pd(m1.m, m8.m);
    }
    result += m8.d[0];
    result += m8.d[1];

    for (unsigned long index = quad_end ; index < a.size() ; index++)
    {
        result += a.elements()[index] * b.elements()[index];
    }
    return result;

}
