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

#include <libla/element_product.hh>

#include <xmmintrin.h>
#include <emmintrin.h>

using namespace honei;

DenseVector<float> & ElementProduct<tags::CPU::SSE>::value(DenseVector<float> & a, const DenseVector<float> & b)
{
    CONTEXT("When multiplyign DenseVector<float> with DenseVector<float> elementwise with SSE:");


    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());


    __m128 m1, m2;

    unsigned long quad_end(a.size() - (b.size() % 4));

    for (unsigned long index = 0 ; index < quad_end ; index += 4) 
    {
        m1 = _mm_load_ps(a.elements() + index);
        m2 = _mm_load_ps(b.elements() + index);

        m1 = _mm_mul_ps(m1, m2);

        _mm_stream_ps(a.elements() + index, m1);

    }

    for (unsigned long index = quad_end ; index < a.size() ; index++)
    {
        a.elements()[index] *= b.elements()[index];
    }
    return a;
}

DenseVector<double> & ElementProduct<tags::CPU::SSE>::value(DenseVector<double> & a, const DenseVector<double> & b)
{
    CONTEXT("When multiplyign DenseVector<double> with DenseVector<double> elementwise with SSE:");


    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());


    __m128d m1, m2;

    unsigned long quad_end(a.size() - (b.size() % 2));

    for (unsigned long index = 0 ; index < quad_end ; index += 2) 
    {
        m1 = _mm_load_pd(a.elements() + index);
        m2 = _mm_load_pd(b.elements() + index);

        m1 = _mm_mul_pd(m1, m2);

        _mm_stream_pd(a.elements() + index, m1);
    }

    for (unsigned long index = quad_end ; index < a.size() ; index++)
    {
        a.elements()[index] *= b.elements()[index];
    }
    return a;
}
