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
        void product_dm(float * x, float * y, float b, unsigned long size)
        {
            _mm_prefetch(y, _MM_HINT_T0);
            _mm_prefetch(x, _MM_HINT_T0);

            __m128 m1, m2, m3, m4, m5, m6, m8;
            float HONEI_ALIGNED(16) b_data;
            b_data = b;
            m8 = _mm_load1_ps(&b_data);

            unsigned long x_address = reinterpret_cast<unsigned long>(x);
            unsigned long x_offset = x_address % 16;
            unsigned long y_address = reinterpret_cast<unsigned long>(y);
            unsigned long y_offset = y_address % 16;

            unsigned long z_offset(x_offset / 4);
            z_offset = (4 - z_offset) % 4;

            unsigned long quad_start = z_offset;
            unsigned long quad_end(size - ((size - quad_start) % 12));

            if (size < 16)
            {
                quad_end = 0;
                quad_start = 0;
            }

            if (x_offset == y_offset)
            {
                for (unsigned long index(quad_start) ; index < quad_end ; index += 12)
                {

                    m1 = _mm_load_ps(x + index);
                    m3 = _mm_load_ps(x + index + 4);
                    m5 = _mm_load_ps(x + index + 8);
                    m2 = _mm_load_ps(y + index);
                    m4 = _mm_load_ps(y + index + 4);
                    m6 = _mm_load_ps(y + index + 8);

                    m2 = _mm_mul_ps(m2, m8);
                    m4 = _mm_mul_ps(m4, m8);
                    m6 = _mm_mul_ps(m6, m8);

                    m1 = _mm_add_ps(m1, m2);
                    m3 = _mm_add_ps(m3, m4);
                    m5 = _mm_add_ps(m5, m6);

                    _mm_store_ps(x + index, m1);
                    _mm_store_ps(x + index + 4, m3);
                    _mm_store_ps(x + index + 8, m5);
                }
            }
            else
            {
                for (unsigned long index(quad_start) ; index < quad_end ; index += 12)
                {

                    m1 = _mm_load_ps(x + index);
                    m3 = _mm_load_ps(x + index + 4);
                    m5 = _mm_load_ps(x + index + 8);
                    m2 = _mm_loadu_ps(y + index);
                    m4 = _mm_loadu_ps(y + index + 4);
                    m6 = _mm_loadu_ps(y + index + 8);

                    m2 = _mm_mul_ps(m2, m8);
                    m4 = _mm_mul_ps(m4, m8);
                    m6 = _mm_mul_ps(m6, m8);

                    m1 = _mm_add_ps(m1, m2);
                    m3 = _mm_add_ps(m3, m4);
                    m5 = _mm_add_ps(m5, m6);

                    _mm_store_ps(x + index, m1);
                    _mm_store_ps(x + index + 4, m3);
                    _mm_store_ps(x + index + 8, m5);
                }
            }

            for (unsigned long index(0) ; index < quad_start ; index++)
            {
                x[index] += y[index] * b;
            }

            for (unsigned long index(quad_end) ; index < size ; index++)
            {
                x[index] += y[index] * b;
            }
        }

        void product_dm(double * x, double * y, double b, unsigned long size)
        {
            _mm_prefetch(y, _MM_HINT_T0);
            _mm_prefetch(x, _MM_HINT_T0);

            __m128d m1, m2, m3, m4, m5, m6, m8;
            double HONEI_ALIGNED(16) b_data;
            b_data = b;
            m8 = _mm_load1_pd(&b_data);

            unsigned long x_address = (unsigned long)x;
            unsigned long x_offset = x_address % 16;
            unsigned long y_address = (unsigned long)y;
            unsigned long y_offset = y_address % 16;

            unsigned long z_offset(x_offset / 8);

            unsigned long quad_start = z_offset;
            unsigned long quad_end(size - ((size - quad_start) % 6));
            if (size < 12)
            {
                quad_end = 0;
                quad_start = 0;
            }

            if (x_offset == y_offset)
            {
                for (unsigned long index = quad_start ; index < quad_end ; index += 6)
                {

                    m1 = _mm_load_pd(x + index);
                    m3 = _mm_load_pd(x + index + 2);
                    m5 = _mm_load_pd(x + index + 4);
                    m2 = _mm_load_pd(y + index);
                    m4 = _mm_load_pd(y + index + 2);
                    m6 = _mm_load_pd(y + index + 4);

                    m2 = _mm_mul_pd(m2, m8);
                    m4 = _mm_mul_pd(m4, m8);
                    m6 = _mm_mul_pd(m6, m8);

                    m1 = _mm_add_pd(m1, m2);
                    m3 = _mm_add_pd(m3, m4);
                    m5 = _mm_add_pd(m5, m6);

                    _mm_store_pd(x + index, m1);
                    _mm_store_pd(x + index + 2, m3);
                    _mm_store_pd(x + index + 4, m5);
                }
            }
            else
            {
                for (unsigned long index = quad_start ; index < quad_end ; index += 6)
                {

                    m1 = _mm_load_pd(x + index);
                    m3 = _mm_load_pd(x + index + 2);
                    m5 = _mm_load_pd(x + index + 4);
                    m2 = _mm_loadu_pd(y + index);
                    m4 = _mm_loadu_pd(y + index + 2);
                    m6 = _mm_loadu_pd(y + index + 4);

                    m2 = _mm_mul_pd(m2, m8);
                    m4 = _mm_mul_pd(m4, m8);
                    m6 = _mm_mul_pd(m6, m8);

                    m1 = _mm_add_pd(m1, m2);
                    m3 = _mm_add_pd(m3, m4);
                    m5 = _mm_add_pd(m5, m6);

                    _mm_store_pd(x + index, m1);
                    _mm_store_pd(x + index + 2, m3);
                    _mm_store_pd(x + index + 4, m5);
                }
            }
            for (unsigned long index(0) ; index < quad_start ; index++)
            {
                x[index] += y[index] * b;
            }
            for (unsigned long index = quad_end ; index < size ; index++)
            {
                x[index] += y[index] * b;
            }
        }

        void product_dm_nx2(float * result, const float * a, const float * b, unsigned long size)
        {
            /*
               float HONEI_ALIGNED(16) result1(0);
               float HONEI_ALIGNED(16) result2(0);

               union sse4
               {
               __m128 m;
               float f[4];
               } m1, m2, m3, m4, m5, m6, m8;
               m8.m = _mm_setzero_ps();

               unsigned long a_address = (unsigned long)a;
               unsigned long a_offset = a_address % 16;

               unsigned long x_offset(a_offset / 4);
               x_offset = (4 - x_offset) % 4;

               unsigned long quad_start = x_offset;
               unsigned long quad_end(size - ((size - quad_start) % 2));

               if (size < 24)
               {
               quad_end = 0;
               quad_start = 0;
               }

               for (unsigned long index(0) ; index < quad_start ; ++index)
               {
               result1 += a[index] * b[index * 2];
               result2 += a[index] * b[index * 2 + 1];
               }
               for (unsigned long index(quad_start) ; index < quad_end ; index += 2)
               {
               m1.m = _mm_set_ps(a[index +1 ], a[index + 1], a[index], a[index]);
               m2.m = _mm_load_ps(b + index * 2);
               m2.m = _mm_mul_ps(m1.m, m2.m);
               m8.m = _mm_add_ps(m8.m, m2.m);
               }
               result1 += m8.f[0];
               result2 += m8.f[1];
               result1 += m8.f[2];
               result2 += m8.f[3];
               for (unsigned long index(quad_end) ; index < size ; ++index)
               {
               result1 += a[index] * b[index * 2];
               result2 += a[index] * b[index * 2 + 1];
               }
               result[0] = result1;
               result[1] = result2;*/
            float result1(0);
            float result2(0);

            for (unsigned long row(0) ; row < size ; ++row)
            {
                result1 += a[row] * b[row * 2];
                result2 += a[row] * b[row * 2 + 1];
            }
            result[0] = result1;
            result[1] = result2;
        }

        void product_dm_nx2(double * result, const double * a, const double * b, unsigned long size)
        {
            double result1(0);
            double result2(0);

            for (unsigned long row(0) ; row < size ; ++row)
            {
                result1 += a[row] * b[row * 2];
                result2 += a[row] * b[row * 2 + 1];
            }
            result[0] = result1;
            result[1] = result2;
        }

        void product_bmdv(float * x, const float * y, const float * z, unsigned long size)
        {
            __m128 m1, m2, m3, m4, m5, m6;

            unsigned long x_address = reinterpret_cast<unsigned long>(x);
            unsigned long x_offset = x_address % 16;
            unsigned long z_address = reinterpret_cast<unsigned long>(z);
            unsigned long z_offset = z_address % 16;

            unsigned long w_offset(x_offset / 4);
            w_offset = (4 - w_offset) % 4;

            unsigned long quad_start = w_offset;
            unsigned long quad_end(size - ((size - quad_start) % 8));

            if (size < 16)
            {
                quad_end = 0;
                quad_start = 0;
            }

            if (x_offset == z_offset)
            {
                for (unsigned long index(quad_start) ; index < quad_end ; index += 8)
                {
                    m1 = _mm_load_ps(x + index);
                    m4 = _mm_load_ps(x + index + 4);
                    m2 = _mm_load_ps(y + index);
                    m5 = _mm_load_ps(y + index + 4);
                    m3 = _mm_load_ps(z + index);
                    m6 = _mm_load_ps(z + index + 4);

                    m2 = _mm_mul_ps(m2, m3);
                    m5 = _mm_mul_ps(m5, m6);

                    m1 = _mm_add_ps(m1, m2);
                    m4 = _mm_add_ps(m4, m5);

                    _mm_store_ps(x + index, m1);
                    _mm_store_ps(x + index + 4, m4);
                }
            }
            else
            {
                for (unsigned long index(quad_start) ; index < quad_end ; index += 8)
                {
                    m1 = _mm_load_ps(x + index);
                    m4 = _mm_load_ps(x + index + 4);
                    m2 = _mm_load_ps(y + index);
                    m5 = _mm_load_ps(y + index + 4);
                    m3 = _mm_loadu_ps(z + index);
                    m6 = _mm_loadu_ps(z + index + 4);

                    m2 = _mm_mul_ps(m2, m3);
                    m5 = _mm_mul_ps(m5, m6);

                    m1 = _mm_add_ps(m1, m2);
                    m4 = _mm_add_ps(m4, m5);

                    _mm_store_ps(x + index, m1);
                    _mm_store_ps(x + index + 4, m4);
                }
            }

            for (unsigned long index(0) ; index < quad_start ; index++)
            {
                x[index] += y[index] * z[index];
            }

            for (unsigned long index(quad_end) ; index < size ; index++)
            {
                x[index] += y[index] * z[index];
            }
        }

        void product_bmdv(double * x, const double * y, const double * z, unsigned long size)
        {
            __m128d m1, m2, m3, m4, m5, m6;

            unsigned long x_address = reinterpret_cast<unsigned long>(x);
            unsigned long x_offset = x_address % 16;
            unsigned long z_address = reinterpret_cast<unsigned long>(z);
            unsigned long z_offset = z_address % 16;

            unsigned long w_offset(x_offset / 8);

            unsigned long quad_start = w_offset;
            unsigned long quad_end(size - ((size - quad_start) % 4));

            if (size < 16)
            {
                quad_end = 0;
                quad_start = 0;
            }

            if (x_offset == z_offset)
            {
                for (unsigned long index(quad_start) ; index < quad_end ; index += 4)
                {
                    m1 = _mm_load_pd(x + index);
                    m4 = _mm_load_pd(x + index + 2);
                    m2 = _mm_load_pd(y + index);
                    m5 = _mm_load_pd(y + index + 2);
                    m3 = _mm_load_pd(z + index);
                    m6 = _mm_load_pd(z + index + 2);

                    m2 = _mm_mul_pd(m2, m3);
                    m5 = _mm_mul_pd(m5, m6);

                    m1 = _mm_add_pd(m1, m2);
                    m4 = _mm_add_pd(m4, m5);

                    _mm_store_pd(x + index, m1);
                    _mm_store_pd(x + index + 2, m4);
                }
            }
            else
            {
                for (unsigned long index(quad_start) ; index < quad_end ; index += 4)
                {
                    m1 = _mm_load_pd(x + index);
                    m4 = _mm_load_pd(x + index + 2);
                    m2 = _mm_load_pd(y + index);
                    m5 = _mm_load_pd(y + index + 2);
                    m3 = _mm_loadu_pd(z + index);
                    m6 = _mm_loadu_pd(z + index + 2);

                    m2 = _mm_mul_pd(m2, m3);
                    m5 = _mm_mul_pd(m5, m6);

                    m1 = _mm_add_pd(m1, m2);
                    m4 = _mm_add_pd(m4, m5);

                    _mm_store_pd(x + index, m1);
                    _mm_store_pd(x + index + 2, m4);
                }
            }

            for (unsigned long index(0) ; index < quad_start ; index++)
            {
                x[index] += y[index] * z[index];
            }

            for (unsigned long index(quad_end) ; index < size ; index++)
            {
                x[index] += y[index] * z[index];
            }
        }
    }
}

