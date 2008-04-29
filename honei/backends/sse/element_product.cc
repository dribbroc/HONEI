/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/util/attributes.hh>

#include <xmmintrin.h>
#include <emmintrin.h>

namespace honei
{
    namespace sse
    {
        void element_product(float * a, const float * b, unsigned long size)
        {
            __m128 m1, m2, m3, m4, m5, m6, m7, m8;

            unsigned long a_address = reinterpret_cast<unsigned long>(a);
            unsigned long a_offset = a_address % 16;
            unsigned long b_address = reinterpret_cast<unsigned long>(b);
            unsigned long b_offset = b_address % 16;

            unsigned long x_offset(a_offset / 4);
            x_offset = (4 - x_offset) % 4;

            unsigned long quad_start = x_offset;
            unsigned long quad_end(size - ((size - quad_start) % 16));

            if (size < 24)
            {
                quad_end = 0;
                quad_start = 0;
            }

            if (a_offset == b_offset)
            {
                for (unsigned long index(quad_start) ; index < quad_end ; index += 16)
                {
                    m1 = _mm_load_ps(a + index);
                    m3 = _mm_load_ps(a + index + 4);
                    m5 = _mm_load_ps(a + index + 8);
                    m7 = _mm_load_ps(a + index + 12);
                    m2 = _mm_load_ps(b + index);
                    m4 = _mm_load_ps(b + index + 4);
                    m6 = _mm_load_ps(b + index + 8);
                    m8 = _mm_load_ps(b + index + 12);

                    m1 = _mm_mul_ps(m1, m2);
                    m3 = _mm_mul_ps(m3, m4);
                    m5 = _mm_mul_ps(m5, m6);
                    m7 = _mm_mul_ps(m7, m8);

                    _mm_store_ps(a + index, m1);
                    _mm_store_ps(a + index + 4, m3);
                    _mm_store_ps(a + index + 8, m5);
                    _mm_store_ps(a + index + 12, m7);
                }
            }
            else
            {
                for (unsigned long index(quad_start) ; index < quad_end ; index += 16)
                {
                    m1 = _mm_load_ps(a + index);
                    m3 = _mm_load_ps(a + index + 4);
                    m5 = _mm_load_ps(a + index + 8);
                    m7 = _mm_load_ps(a + index + 12);
                    m2 = _mm_loadu_ps(b + index);
                    m4 = _mm_loadu_ps(b + index + 4);
                    m6 = _mm_loadu_ps(b + index + 8);
                    m8 = _mm_loadu_ps(b + index + 12);

                    m1 = _mm_mul_ps(m1, m2);
                    m3 = _mm_mul_ps(m3, m4);
                    m5 = _mm_mul_ps(m5, m6);
                    m7 = _mm_mul_ps(m7, m8);

                    _mm_store_ps(a + index, m1);
                    _mm_store_ps(a + index + 4, m3);
                    _mm_store_ps(a + index + 8, m5);
                    _mm_store_ps(a + index + 12, m7);
                }
            }

            for (unsigned long index(0) ; index < quad_start ; index++)
            {
                a[index] *= b[index];
            }

            for (unsigned long index(quad_end) ; index < size ; index++)
            {
                a[index] *= b[index];
            }
        }

        void element_product(double * a, const double * b, unsigned long size)
        {
            __m128d m1, m2, m3, m4, m5, m6, m7, m8;

            unsigned long a_address = reinterpret_cast<unsigned long>(a);
            unsigned long a_offset = a_address % 16;
            unsigned long b_address = reinterpret_cast<unsigned long>(b);
            unsigned long b_offset = b_address % 16;

            unsigned long x_offset(a_offset / 8);

            unsigned long quad_start = x_offset;
            unsigned long quad_end(size - ((size - quad_start) % 8));

            if (size < 16)
            {
                quad_end = 0;
                quad_start = 0;
            }

            if (a_offset == b_offset)
            {
                for (unsigned long index(quad_start) ; index < quad_end ; index += 8)
                {
                    m1 = _mm_load_pd(a + index);
                    m3 = _mm_load_pd(a + index + 2);
                    m5 = _mm_load_pd(a + index + 4);
                    m7 = _mm_load_pd(a + index + 6);
                    m2 = _mm_load_pd(b + index);
                    m4 = _mm_load_pd(b + index + 2);
                    m6 = _mm_load_pd(b + index + 4);
                    m8 = _mm_load_pd(b + index + 6);

                    m1 = _mm_mul_pd(m1, m2);
                    m3 = _mm_mul_pd(m3, m4);
                    m5 = _mm_mul_pd(m5, m6);
                    m7 = _mm_mul_pd(m7, m8);

                    _mm_store_pd(a + index, m1);
                    _mm_store_pd(a + index + 2, m3);
                    _mm_store_pd(a + index + 4, m5);
                    _mm_store_pd(a + index + 6, m7);
                }
            }
            else
            {
                for (unsigned long index(quad_start) ; index < quad_end ; index += 8)
                {
                    m1 = _mm_load_pd(a + index);
                    m3 = _mm_load_pd(a + index + 2);
                    m5 = _mm_load_pd(a + index + 4);
                    m7 = _mm_load_pd(a + index + 6);
                    m2 = _mm_loadu_pd(b + index);
                    m4 = _mm_loadu_pd(b + index + 2);
                    m6 = _mm_loadu_pd(b + index + 4);
                    m8 = _mm_loadu_pd(b + index + 6);

                    m1 = _mm_mul_pd(m1, m2);
                    m3 = _mm_mul_pd(m3, m4);
                    m5 = _mm_mul_pd(m5, m6);
                    m7 = _mm_mul_pd(m7, m8);

                    _mm_store_pd(a + index, m1);
                    _mm_store_pd(a + index + 2, m3);
                    _mm_store_pd(a + index + 4, m5);
                    _mm_store_pd(a + index + 6, m7);
                }
            }

            for (unsigned long index(0) ; index < quad_start ; index++)
            {
                a[index] *= b[index];
            }

            for (unsigned long index(quad_end) ; index < size ; index++)
            {
                a[index] *= b[index];
            }
        }
    }
}

