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
            _mm_prefetch((char *)y, _MM_HINT_T0);
            _mm_prefetch((char *)x, _MM_HINT_T0);

            __m128 m1, m2, m3, m4, m5, m6, m8;
            float HONEI_ALIGNED(16) b_data;
            b_data = b;
            m8 = _mm_load1_ps(&b_data);

            unsigned long x_address((unsigned long)x);
            unsigned long x_offset(x_address % 16);
            unsigned long y_address((unsigned long)y);
            unsigned long y_offset(y_address % 16);

            unsigned long z_offset(x_offset / 4);
            z_offset = (4 - z_offset) % 4;

            unsigned long quad_start(z_offset);
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

            for (unsigned long index(0) ; index < quad_start ; ++index)
            {
                x[index] += y[index] * b;
            }

            for (unsigned long index(quad_end) ; index < size ; ++index)
            {
                x[index] += y[index] * b;
            }
        }

        void product_dm(double * x, double * y, double b, unsigned long size)
        {
            _mm_prefetch((char *)y, _MM_HINT_T0);
            _mm_prefetch((char *)x, _MM_HINT_T0);

            __m128d m1, m2, m3, m4, m5, m6, m8;
            double HONEI_ALIGNED(16) b_data;
            b_data = b;
            m8 = _mm_load1_pd(&b_data);

            unsigned long x_address((unsigned long)x);
            unsigned long x_offset(x_address % 16);
            unsigned long y_address((unsigned long)y);
            unsigned long y_offset(y_address % 16);

            unsigned long z_offset(x_offset / 8);

            unsigned long quad_start(z_offset);
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
            for (unsigned long index(0) ; index < quad_start ; ++index)
            {
                x[index] += y[index] * b;
            }
            for (unsigned long index = quad_end ; index < size ; ++index)
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

            unsigned long x_address((unsigned long)x);
            unsigned long x_offset(x_address % 16);
            unsigned long z_address((unsigned long)z);
            unsigned long z_offset(z_address % 16);

            unsigned long w_offset(x_offset / 4);
            w_offset = (4 - w_offset) % 4;

            unsigned long quad_start(w_offset);
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

            for (unsigned long index(0) ; index < quad_start ; ++index)
            {
                x[index] += y[index] * z[index];
            }

            for (unsigned long index(quad_end) ; index < size ; ++index)
            {
                x[index] += y[index] * z[index];
            }
        }

        void product_bmdv(double * x, const double * y, const double * z, unsigned long size)
        {
            __m128d m1, m2, m3, m4, m5, m6;

            unsigned long x_address((unsigned long)x);
            unsigned long x_offset(x_address % 16);
            unsigned long z_address((unsigned long)z);
            unsigned long z_offset(z_address % 16);

            unsigned long w_offset(x_offset / 8);

            unsigned long quad_start(w_offset);
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

            for (unsigned long index(0) ; index < quad_start ; ++index)
            {
                x[index] += y[index] * z[index];
            }

            for (unsigned long index(quad_end) ; index < size ; ++index)
            {
                x[index] += y[index] * z[index];
            }
        }

        void product_bmdv_q1(float * ll, float * ld, float * lu,
                float * dl, float * dd, float * du,
                float * ul, float * ud, float * uu,
                float * b, float * result,
                unsigned long size, unsigned long m)
        {
            __m128 m1, m2, m3;

            unsigned long index(0);

            //index 0
            result[index] = dd[index] * b[index]
                + du[index] * b[index + 1]
                + ul[index] * b[index + m - 1]
                + ud[index] * b[index + m]
                + uu[index] * b[index + m + 1];

            //index in [1, root_n -1[
            index = 1;
            unsigned long quad_start(index + 4 - (index % 4));
            unsigned long quad_end(m - 1 - ((m - 1 - quad_start) % 4));

            for (unsigned long si(index) ; si < quad_start ; ++si)
            {
                result[si] = dd[si] * b[si]
                    + dl[si] * b[si - 1]
                    + du[si] * b[si + 1]
                    + ul[si] * b[si + m - 1]
                    + ud[si] * b[si + m]
                    + uu[si] * b[si + m + 1];
            }

            for (unsigned long si(quad_start) ; si < quad_end ; si+=4)
            {
                //result[si] = dd[si] * b[si]
                m2 = _mm_loadu_ps(b + si);
                m3 = _mm_load_ps(dd + si);
                m1 = _mm_mul_ps(m3, m2);

                //+ dl[si] * b[si - 1]
                m2 = _mm_loadu_ps(b + si - 1);
                m3 = _mm_load_ps(dl + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ du[si] * b[si + 1]
                m2 = _mm_loadu_ps(b + si + 1);
                m3 = _mm_load_ps(du + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ ul[si] * b[si + m - 1]
                m2 = _mm_loadu_ps(b + si + m - 1);
                m3 = _mm_load_ps(ul + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ ud[si] * b[si + m]
                m2 = _mm_loadu_ps(b + si + m);
                m3 = _mm_load_ps(ud + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ uu[si] * b[si + m + 1];
                m2 = _mm_loadu_ps(b + si + m + 1);
                m3 = _mm_load_ps(uu + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                _mm_store_ps(result + si, m1);
            }

            for (unsigned long si(quad_end) ; si < m - 1 ; ++si)
            {
                result[si] = dd[si] * b[si]
                    + dl[si] * b[si - 1]
                    + du[si] * b[si + 1]
                    + ul[si] * b[si + m - 1]
                    + ud[si] * b[si + m]
                    + uu[si] * b[si + m + 1];
            }


            //index root_n -1
            index = m - 1;
            result[index] = dd[index] * b[index]
                + lu[index] * b[index - m + 1]
                + dl[index] * b[index - 1]
                + du[index] * b[index + 1]
                + ul[index] * b[index + m - 1]
                + ud[index] * b[index + m]
                + uu[index] * b[index + m + 1];

            //index root_n
            index = m;
            result[index] = dd[index] * b[index]
                + lu[index] * b[index - m + 1]
                + ld[index] * b[index - m]
                + dl[index] * b[index - 1]
                + du[index] * b[index + 1]
                + ul[index] * b[index + m - 1]
                + ud[index] * b[index + m]
                + uu[index] * b[index + m + 1];

            //index in [root_n + 1, n - (root_n + 1)[
            index = m + 1;
            quad_start = index + 4 - (index % 4);
            quad_end = size - m - 1 - ((size - m - 1 - quad_start) % 4);

            for (unsigned long si(index) ; si < quad_start ; ++si)
            {
                result[si] = dd[si] * b[si]
                    + ll[si] * b[si - m - 1]
                    + lu[si] * b[si - m + 1]
                    + ld[si] * b[si - m]
                    + dl[si] * b[si - 1]
                    + du[si] * b[si + 1]
                    + ul[si] * b[si + m - 1]
                    + ud[si] * b[si + m]
                    + uu[si] * b[si + m + 1];
            }

            for (unsigned long si(quad_start) ; si < quad_end ; si+=4)
            {
                //result[si] = dd[si] * b[si]
                m2 = _mm_loadu_ps(b + si);
                m3 = _mm_load_ps(dd + si);
                m1 = _mm_mul_ps(m3, m2);

                //+ ll[si] * b[si - m - 1]
                m2 = _mm_loadu_ps(b + si -m - 1);
                m3 = _mm_load_ps(ll + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ lu[si] * b[si - m + 1]
                m2 = _mm_loadu_ps(b + si - m + 1);
                m3 = _mm_load_ps(lu + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ ld[si] * b[si - m]
                m2 = _mm_loadu_ps(b + si - m);
                m3 = _mm_load_ps(ld + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ dl[si] * b[si - 1]
                m2 = _mm_loadu_ps(b + si - 1);
                m3 = _mm_load_ps(dl + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ du[si] * b[si + 1]
                m2 = _mm_loadu_ps(b + si + 1);
                m3 = _mm_load_ps(du + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ ul[si] * b[si + m - 1]
                m2 = _mm_loadu_ps(b + si + m - 1);
                m3 = _mm_load_ps(ul + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ ud[si] * b[si + m]
                m2 = _mm_loadu_ps(b + si + m);
                m3 = _mm_load_ps(ud + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ uu[si] * b[si + m + 1];
                m2 = _mm_loadu_ps(b + si + m + 1);
                m3 = _mm_load_ps(uu + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                _mm_store_ps(result + si, m1);
            }

            for (unsigned long si(quad_end) ; si < size - m - 1 ; ++si)
            {
                result[si] = dd[si] * b[si]
                    + ll[si] * b[si - m - 1]
                    + lu[si] * b[si - m + 1]
                    + ld[si] * b[si - m]
                    + dl[si] * b[si - 1]
                    + du[si] * b[si + 1]
                    + ul[si] * b[si + m - 1]
                    + ud[si] * b[si + m]
                    + uu[si] * b[si + m + 1];
            }

            //index n - (root_n + 1)
            index = size - m - 1;
            result[index] = dd[index] * b[index]
                + ll[index] * b[index - m - 1]
                + lu[index] * b[index - m + 1]
                + ld[index] * b[index - m]
                + dl[index] * b[index - 1]
                + du[index] * b[index + 1]
                + ul[index] * b[index + m - 1]
                + ud[index] * b[index + m];

            //index n - root_n
            index = size - m;
            result[index] = dd[index] * b[index]
                + ll[index] * b[index - m - 1]
                + lu[index] * b[index - m + 1]
                + ld[index] * b[index - m]
                + dl[index] * b[index - 1]
                + du[index] * b[index + 1]
                + ul[index] * b[index + m - 1];

            //index in [n - root_n + 1, n -1[
            index = size - m + 1;
            quad_start = index + 4 - (index % 4);
            quad_end = size - 1 - ((size - 1 - quad_start) % 4);

            for (unsigned long si(index) ; si < quad_start ; ++si)
            {
                result[si] = dd[si] * b[si]
                    + ll[si] * b[si - m - 1]
                    + lu[si] * b[si - m + 1]
                    + ld[si] * b[si - m]
                    + dl[si] * b[si - 1]
                    + du[si] * b[si + 1];
            }

            for (unsigned long si(quad_start) ; si < quad_end ; si+=4)
            {
                //result[si] = dd[si] * b[si]
                m2 = _mm_loadu_ps(b + si);
                m3 = _mm_load_ps(dd + si);
                m1 = _mm_mul_ps(m3, m2);

                //+ ll[si] * b[si - m - 1]
                m2 = _mm_loadu_ps(b + si -m - 1);
                m3 = _mm_load_ps(ll + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ lu[si] * b[si - m + 1]
                m2 = _mm_loadu_ps(b + si - m + 1);
                m3 = _mm_load_ps(lu + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ ld[si] * b[si - m]
                m2 = _mm_loadu_ps(b + si - m);
                m3 = _mm_load_ps(ld + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ dl[si] * b[si - 1]
                m2 = _mm_loadu_ps(b + si - 1);
                m3 = _mm_load_ps(dl + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ du[si] * b[si + 1]
                m2 = _mm_loadu_ps(b + si + 1);
                m3 = _mm_load_ps(du + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                _mm_store_ps(result + si, m1);
            }

            for (unsigned long si(quad_end) ; si < size - m - 1 ; ++si)
            {
                result[si] = dd[si] * b[si]
                    + ll[si] * b[si - m - 1]
                    + lu[si] * b[si - m + 1]
                    + ld[si] * b[si - m]
                    + dl[si] * b[si - 1]
                    + du[si] * b[si + 1];
            }

            //index n - 1
            index = size - 1;
            result[index] = dd[index] * b[index]
                + ll[index] * b[index - m - 1]
                + lu[index] * b[index - m + 1]
                + ld[index] * b[index - m]
                + dl[index] * b[index - 1];
        }

        void product_bmdv_q1(double * ll, double * ld, double * lu,
                double * dl, double * dd, double * du,
                double * ul, double * ud, double * uu,
                double * b, double * result,
                unsigned long size, unsigned long m)
        {
            __m128d m1, m2, m3;

            unsigned long index(0);

            //index 0
            result[index] = dd[index] * b[index]
                + du[index] * b[index + 1]
                + ul[index] * b[index + m - 1]
                + ud[index] * b[index + m]
                + uu[index] * b[index + m + 1];

            //index in [1, root_n -1[
            index = 1;
            unsigned long quad_start(index + (index % 2));
            unsigned long quad_end(m - 1 - ((m - 1 - quad_start) % 2));

            for (unsigned long si(index) ; si < quad_start ; ++si)
            {
                result[si] = dd[si] * b[si]
                    + dl[si] * b[si - 1]
                    + du[si] * b[si + 1]
                    + ul[si] * b[si + m - 1]
                    + ud[si] * b[si + m]
                    + uu[si] * b[si + m + 1];
            }

            for (unsigned long si(quad_start) ; si < quad_end ; si+=2)
            {
                //result[si] = dd[si] * b[si]
                m2 = _mm_loadu_pd(b + si);
                m3 = _mm_load_pd(dd + si);
                m1 = _mm_mul_pd(m3, m2);

                //+ dl[si] * b[si - 1]
                m2 = _mm_loadu_pd(b + si - 1);
                m3 = _mm_load_pd(dl + si);
                m3 = _mm_mul_pd(m3, m2);
                m1 = _mm_add_pd(m1, m3);

                //+ du[si] * b[si + 1]
                m2 = _mm_loadu_pd(b + si + 1);
                m3 = _mm_load_pd(du + si);
                m3 = _mm_mul_pd(m3, m2);
                m1 = _mm_add_pd(m1, m3);

                //+ ul[si] * b[si + m - 1]
                m2 = _mm_loadu_pd(b + si + m - 1);
                m3 = _mm_load_pd(ul + si);
                m3 = _mm_mul_pd(m3, m2);
                m1 = _mm_add_pd(m1, m3);

                //+ ud[si] * b[si + m]
                m2 = _mm_loadu_pd(b + si + m);
                m3 = _mm_load_pd(ud + si);
                m3 = _mm_mul_pd(m3, m2);
                m1 = _mm_add_pd(m1, m3);

                //+ uu[si] * b[si + m + 1];
                m2 = _mm_loadu_pd(b + si + m + 1);
                m3 = _mm_load_pd(uu + si);
                m3 = _mm_mul_pd(m3, m2);
                m1 = _mm_add_pd(m1, m3);

                _mm_store_pd(result + si, m1);
            }

            for (unsigned long si(quad_end) ; si < m - 1 ; ++si)
            {
                result[si] = dd[si] * b[si]
                    + dl[si] * b[si - 1]
                    + du[si] * b[si + 1]
                    + ul[si] * b[si + m - 1]
                    + ud[si] * b[si + m]
                    + uu[si] * b[si + m + 1];
            }


            //index root_n -1
            index = m - 1;
            result[index] = dd[index] * b[index]
                + lu[index] * b[index - m + 1]
                + dl[index] * b[index - 1]
                + du[index] * b[index + 1]
                + ul[index] * b[index + m - 1]
                + ud[index] * b[index + m]
                + uu[index] * b[index + m + 1];

            //index root_n
            index = m;
            result[index] = dd[index] * b[index]
                + lu[index] * b[index - m + 1]
                + ld[index] * b[index - m]
                + dl[index] * b[index - 1]
                + du[index] * b[index + 1]
                + ul[index] * b[index + m - 1]
                + ud[index] * b[index + m]
                + uu[index] * b[index + m + 1];

            //index in [root_n + 1, n - (root_n + 1)[
            index = m + 1;
            quad_start = index + (index % 2);
            quad_end = size - m - 1 - ((size - m - 1 - quad_start) % 2);

            for (unsigned long si(index) ; si < quad_start ; ++si)
            {
                result[si] = dd[si] * b[si]
                    + ll[si] * b[si - m - 1]
                    + lu[si] * b[si - m + 1]
                    + ld[si] * b[si - m]
                    + dl[si] * b[si - 1]
                    + du[si] * b[si + 1]
                    + ul[si] * b[si + m - 1]
                    + ud[si] * b[si + m]
                    + uu[si] * b[si + m + 1];
            }

            for (unsigned long si(quad_start) ; si < quad_end ; si+=2)
            {
                //result[si] = dd[si] * b[si]
                m2 = _mm_loadu_pd(b + si);
                m3 = _mm_load_pd(dd + si);
                m1 = _mm_mul_pd(m3, m2);

                //+ ll[si] * b[si - m - 1]
                m2 = _mm_loadu_pd(b + si -m - 1);
                m3 = _mm_load_pd(ll + si);
                m3 = _mm_mul_pd(m3, m2);
                m1 = _mm_add_pd(m1, m3);

                //+ lu[si] * b[si - m + 1]
                m2 = _mm_loadu_pd(b + si - m + 1);
                m3 = _mm_load_pd(lu + si);
                m3 = _mm_mul_pd(m3, m2);
                m1 = _mm_add_pd(m1, m3);

                //+ ld[si] * b[si - m]
                m2 = _mm_loadu_pd(b + si - m);
                m3 = _mm_load_pd(ld + si);
                m3 = _mm_mul_pd(m3, m2);
                m1 = _mm_add_pd(m1, m3);

                //+ dl[si] * b[si - 1]
                m2 = _mm_loadu_pd(b + si - 1);
                m3 = _mm_load_pd(dl + si);
                m3 = _mm_mul_pd(m3, m2);
                m1 = _mm_add_pd(m1, m3);

                //+ du[si] * b[si + 1]
                m2 = _mm_loadu_pd(b + si + 1);
                m3 = _mm_load_pd(du + si);
                m3 = _mm_mul_pd(m3, m2);
                m1 = _mm_add_pd(m1, m3);

                //+ ul[si] * b[si + m - 1]
                m2 = _mm_loadu_pd(b + si + m - 1);
                m3 = _mm_load_pd(ul + si);
                m3 = _mm_mul_pd(m3, m2);
                m1 = _mm_add_pd(m1, m3);

                //+ ud[si] * b[si + m]
                m2 = _mm_loadu_pd(b + si + m);
                m3 = _mm_load_pd(ud + si);
                m3 = _mm_mul_pd(m3, m2);
                m1 = _mm_add_pd(m1, m3);

                //+ uu[si] * b[si + m + 1];
                m2 = _mm_loadu_pd(b + si + m + 1);
                m3 = _mm_load_pd(uu + si);
                m3 = _mm_mul_pd(m3, m2);
                m1 = _mm_add_pd(m1, m3);

                _mm_store_pd(result + si, m1);
            }

            for (unsigned long si(quad_end) ; si < size - m - 1 ; ++si)
            {
                result[si] = dd[si] * b[si]
                    + ll[si] * b[si - m - 1]
                    + lu[si] * b[si - m + 1]
                    + ld[si] * b[si - m]
                    + dl[si] * b[si - 1]
                    + du[si] * b[si + 1]
                    + ul[si] * b[si + m - 1]
                    + ud[si] * b[si + m]
                    + uu[si] * b[si + m + 1];
            }

            //index n - (root_n + 1)
            index = size - m - 1;
            result[index] = dd[index] * b[index]
                + ll[index] * b[index - m - 1]
                + lu[index] * b[index - m + 1]
                + ld[index] * b[index - m]
                + dl[index] * b[index - 1]
                + du[index] * b[index + 1]
                + ul[index] * b[index + m - 1]
                + ud[index] * b[index + m];

            //index n - root_n
            index = size - m;
            result[index] = dd[index] * b[index]
                + ll[index] * b[index - m - 1]
                + lu[index] * b[index - m + 1]
                + ld[index] * b[index - m]
                + dl[index] * b[index - 1]
                + du[index] * b[index + 1]
                + ul[index] * b[index + m - 1];

            //index in [n - root_n + 1, n -1[
            index = size - m + 1;
            quad_start = index + (index % 2);
            quad_end = size - 1 - ((size - 1 - quad_start) % 2);

            for (unsigned long si(index) ; si < quad_start ; ++si)
            {
                result[si] = dd[si] * b[si]
                    + ll[si] * b[si - m - 1]
                    + lu[si] * b[si - m + 1]
                    + ld[si] * b[si - m]
                    + dl[si] * b[si - 1]
                    + du[si] * b[si + 1];
            }

            for (unsigned long si(quad_start) ; si < quad_end ; si+=2)
            {
                //result[si] = dd[si] * b[si]
                m2 = _mm_loadu_pd(b + si);
                m3 = _mm_load_pd(dd + si);
                m1 = _mm_mul_pd(m3, m2);

                //+ ll[si] * b[si - m - 1]
                m2 = _mm_loadu_pd(b + si -m - 1);
                m3 = _mm_load_pd(ll + si);
                m3 = _mm_mul_pd(m3, m2);
                m1 = _mm_add_pd(m1, m3);

                //+ lu[si] * b[si - m + 1]
                m2 = _mm_loadu_pd(b + si - m + 1);
                m3 = _mm_load_pd(lu + si);
                m3 = _mm_mul_pd(m3, m2);
                m1 = _mm_add_pd(m1, m3);

                //+ ld[si] * b[si - m]
                m2 = _mm_loadu_pd(b + si - m);
                m3 = _mm_load_pd(ld + si);
                m3 = _mm_mul_pd(m3, m2);
                m1 = _mm_add_pd(m1, m3);

                //+ dl[si] * b[si - 1]
                m2 = _mm_loadu_pd(b + si - 1);
                m3 = _mm_load_pd(dl + si);
                m3 = _mm_mul_pd(m3, m2);
                m1 = _mm_add_pd(m1, m3);

                //+ du[si] * b[si + 1]
                m2 = _mm_loadu_pd(b + si + 1);
                m3 = _mm_load_pd(du + si);
                m3 = _mm_mul_pd(m3, m2);
                m1 = _mm_add_pd(m1, m3);

                _mm_store_pd(result + si, m1);
            }

            for (unsigned long si(quad_end) ; si < size - m - 1 ; ++si)
            {
                result[si] = dd[si] * b[si]
                    + ll[si] * b[si - m - 1]
                    + lu[si] * b[si - m + 1]
                    + ld[si] * b[si - m]
                    + dl[si] * b[si - 1]
                    + du[si] * b[si + 1];
            }

            //index n - 1
            index = size - 1;
            result[index] = dd[index] * b[index]
                + ll[index] * b[index - m - 1]
                + lu[index] * b[index - m + 1]
                + ld[index] * b[index - m]
                + dl[index] * b[index - 1];
        }

        void product_smell_dv(float * result, const unsigned long * Aj, const float * Ax, const unsigned long * Arl, const float * b,
                unsigned long stride, unsigned long /*rows*/, unsigned long /*num_cols_per_row*/,
                unsigned long row_start, unsigned long row_end, const unsigned long threads)
        {
            if (threads % 4 != 0)
            {
                for (unsigned long row(row_start) ; row < row_end ; ++row)
                {
                    const unsigned long * tAj(Aj);
                    const float * tAx(Ax);
                    float sum(0);
                    tAj += row * threads;
                    tAx += row * threads;

                    const unsigned long max(Arl[row]);
                    for(unsigned long n = 0; n < max ; n++)
                    {
                        for (unsigned long thread(0) ; thread < threads ; ++thread)
                        {
                            const float A_ij = *(tAx + thread);

                            const unsigned long col = *(tAj + thread);
                            sum += A_ij * b[col];
                        }

                        tAj += stride;
                        tAx += stride;
                    }
                    result[row] = sum;
                }
            }
            else
            {
                for (unsigned long row(row_start) ; row < row_end ; ++row)
                {
                    const unsigned long * tAj(Aj);
                    const float * tAx(Ax);
                    tAj += row * threads;
                    tAx += row * threads;

                    const unsigned long max(Arl[row]);
                    __m128 A_ij, b_v;
                    union sse4
                    {
                        __m128 m;
                        float f[4];
                    } sum_v;
                    sum_v.m = _mm_setzero_ps();

                    for(unsigned long n = 0; n < max ; n++)
                    {
                        for (unsigned long thread(0) ; thread < threads ; thread+=4)
                        {
                            A_ij = _mm_load_ps(tAx + thread);
                            const unsigned long col0 = *(tAj + thread);
                            const unsigned long col1 = *(tAj + thread + 1);
                            const unsigned long col2 = *(tAj + thread + 2);
                            const unsigned long col3 = *(tAj + thread + 3);
                            b_v = _mm_set_ps(*(b+col3), *(b+col2), *(b+col1), *(b+col0));
                            b_v = _mm_mul_ps(A_ij, b_v);
                            sum_v.m = _mm_add_ps(b_v, sum_v.m);
                        }

                        tAj += stride;
                        tAx += stride;
                    }
                    float sum = sum_v.f[0] + sum_v.f[1] + sum_v.f[2] + sum_v.f[3];

                    result[row] = sum;
                }
            }
        }

        void product_smell_dv(double * result, const unsigned long * Aj, const double * Ax, const unsigned long * Arl, const double * b,
                unsigned long stride, unsigned long /*rows*/, unsigned long /*num_cols_per_row*/,
                unsigned long row_start, unsigned long row_end, const unsigned long threads)
        {
            if (threads % 2 != 0)
            {
                for (unsigned long row(row_start) ; row < row_end ; ++row)
                {
                    const unsigned long * tAj(Aj);
                    const double * tAx(Ax);
                    double sum(0);
                    tAj += row * threads;
                    tAx += row * threads;

                    const unsigned long max(Arl[row]);
                    for(unsigned long n = 0; n < max ; n++)
                    {
                        for (unsigned long thread(0) ; thread < threads ; ++thread)
                        {
                            const double A_ij = *(tAx + thread);

                            const unsigned long col = *(tAj + thread);
                            sum += A_ij * b[col];
                        }

                        tAj += stride;
                        tAx += stride;
                    }
                    result[row] = sum;
                }
            }
            else
            {
                for (unsigned long row(row_start) ; row < row_end ; ++row)
                {
                    const unsigned long * tAj(Aj);
                    const double * tAx(Ax);
                    tAj += row * threads;
                    tAx += row * threads;

                    const unsigned long max(Arl[row]);
                    __m128d A_ij, b_v;
                    union sse2
                    {
                        __m128d m;
                        double d[2];
                    } sum_v;
                    sum_v.m = _mm_setzero_pd();

                    for(unsigned long n = 0; n < max ; n++)
                    {
                        for (unsigned long thread(0) ; thread < threads ; thread+=2)
                        {
                            A_ij = _mm_load_pd(tAx + thread);
                            const unsigned long col0 = *(tAj + thread);
                            const unsigned long col1 = *(tAj + thread + 1);
                            b_v = _mm_set_pd(*(b+col1), *(b+col0));
                            b_v = _mm_mul_pd(A_ij, b_v);
                            sum_v.m = _mm_add_pd(b_v, sum_v.m);
                        }

                        tAj += stride;
                        tAx += stride;
                    }
                    double sum = sum_v.d[0] + sum_v.d[1];

                    result[row] = sum;
                }
            }
        }

        void product_csr_dv(float * result, const unsigned long * Aj, const float * Ax, const unsigned long * Ar, const float * b,
                unsigned long blocksize, unsigned long row_start, unsigned long row_end)
        {
            if (blocksize != 4)
            {
                switch(blocksize)
                {
                    case 1:
                        for (unsigned long row(row_start) ; row < row_end ; ++row)
                        {
                            float sum(0);
                            const unsigned long end(Ar[row+1]);
                            for (unsigned long i(Ar[row]) ; i < end ; ++i)
                            {
                                sum += Ax[i] * b[Aj[i]];
                            }
                            result[row] = sum;
                        }
                        break;
                    case 2:
                        for (unsigned long row(row_start) ; row < row_end ; ++row)
                        {
                            float sum(0);
                            const unsigned long end(Ar[row+1]);
                            for (unsigned long i(Ar[row]) ; i < end ; ++i)
                            {
                                const unsigned long col(Aj[i]);
                                sum += Ax[2*i] * b[col];
                                sum += Ax[2*i+1] * b[col+1];
                            }
                            result[row] = sum;
                        }
                        break;
                    /*case 4:
                        for (unsigned long row(row_start) ; row < row_end ; ++row)
                        {
                            float sum(0);
                            const unsigned long end(Ar[row+1]);
                            for (unsigned long i(Ar[row]) ; i < end ; ++i)
                            {
                                const unsigned long col(Aj[i]);
                                sum += Ax[2*i] * b[col];
                                sum += Ax[2*i+1] * b[col+1];
                                sum += Ax[2*i+2] * b[col+2];
                                sum += Ax[2*i+3] * b[col+3];
                            }
                            result[row] = sum;
                        }
                        break;*/
                    case 8:
                        for (unsigned long row(row_start) ; row < row_end ; ++row)
                        {
                            float sum(0);
                            const unsigned long end(Ar[row+1]);
                            for (unsigned long i(Ar[row]) ; i < end ; ++i)
                            {
                                const unsigned long col(Aj[i]);
                                sum += Ax[2*i] * b[col];
                                sum += Ax[2*i+1] * b[col+1];
                                sum += Ax[2*i+2] * b[col+2];
                                sum += Ax[2*i+3] * b[col+3];
                                sum += Ax[2*i+4] * b[col+4];
                                sum += Ax[2*i+5] * b[col+5];
                                sum += Ax[2*i+6] * b[col+6];
                                sum += Ax[2*i+7] * b[col+7];
                            }
                            result[row] = sum;
                        }
                        break;
                    default:
                        for (unsigned long row(row_start) ; row < row_end ; ++row)
                        {
                            float sum(0);
                            const unsigned long end(Ar[row+1]);
                            for (unsigned long i(Ar[row]) ; i < end ; ++i)
                            {
                                for (unsigned long blocki(0) ; blocki < blocksize ; ++blocki)
                                {
                                    sum += Ax[(i*blocksize)+blocki] * b[Aj[i] + blocki];
                                }
                            }
                            result[row] = sum;
                        }
                        break;
                }
            }
            else
            {
                __m128 A_ij, b_v;
                union sse4
                {
                    __m128 m;
                    float f[4];
                } sum_v;

                for (unsigned long row(row_start) ; row < row_end ; ++row)
                {
                    sum_v.m = _mm_setzero_ps();

                    const unsigned long end(Ar[row+1]);
                    for (unsigned long i(Ar[row]) ; i < end ; i++)
                    {
                        A_ij = _mm_load_ps(Ax + (i * blocksize));
                        const unsigned long col = Aj[i];
                        b_v = _mm_loadu_ps(b+col);
                        b_v = _mm_mul_ps(A_ij, b_v);
                        sum_v.m = _mm_add_ps(b_v, sum_v.m);
                    }
                    float sum = sum_v.f[0] + sum_v.f[1] + sum_v.f[2] + sum_v.f[3];

                    result[row] = sum;
                }
            }
        }

        void product_csr_dv(double * result, const unsigned long * Aj, const double * Ax, const unsigned long * Ar, const double * b,
                unsigned long blocksize, unsigned long row_start, unsigned long row_end)
        {
            if (blocksize != 2)
            {
                switch(blocksize)
                {
                    case 1:
                        for (unsigned long row(row_start) ; row < row_end ; ++row)
                        {
                            double sum(0);
                            const unsigned long end(Ar[row+1]);
                            for (unsigned long i(Ar[row]) ; i < end ; ++i)
                            {
                                sum += Ax[i] * b[Aj[i]];
                            }
                            result[row] = sum;
                        }
                        break;
                    /*case 2:
                        for (unsigned long row(row_start) ; row < row_end ; ++row)
                        {
                            double sum(0);
                            const unsigned long end(Ar[row+1]);
                            for (unsigned long i(Ar[row]) ; i < end ; ++i)
                            {
                                const unsigned long col(Aj[i]);
                                sum += Ax[2*i] * b[col];
                                sum += Ax[2*i+1] * b[col+1];
                            }
                            result[row] = sum;
                        }
                        break;*/
                    case 4:
                        for (unsigned long row(row_start) ; row < row_end ; ++row)
                        {
                            double sum(0);
                            const unsigned long end(Ar[row+1]);
                            for (unsigned long i(Ar[row]) ; i < end ; ++i)
                            {
                                const unsigned long col(Aj[i]);
                                sum += Ax[2*i] * b[col];
                                sum += Ax[2*i+1] * b[col+1];
                                sum += Ax[2*i+2] * b[col+2];
                                sum += Ax[2*i+3] * b[col+3];
                            }
                            result[row] = sum;
                        }
                        break;
                    case 8:
                        for (unsigned long row(row_start) ; row < row_end ; ++row)
                        {
                            double sum(0);
                            const unsigned long end(Ar[row+1]);
                            for (unsigned long i(Ar[row]) ; i < end ; ++i)
                            {
                                const unsigned long col(Aj[i]);
                                sum += Ax[2*i] * b[col];
                                sum += Ax[2*i+1] * b[col+1];
                                sum += Ax[2*i+2] * b[col+2];
                                sum += Ax[2*i+3] * b[col+3];
                                sum += Ax[2*i+4] * b[col+4];
                                sum += Ax[2*i+5] * b[col+5];
                                sum += Ax[2*i+6] * b[col+6];
                                sum += Ax[2*i+7] * b[col+7];
                            }
                            result[row] = sum;
                        }
                        break;
                    default:
                        for (unsigned long row(row_start) ; row < row_end ; ++row)
                        {
                            double sum(0);
                            const unsigned long end(Ar[row+1]);
                            for (unsigned long i(Ar[row]) ; i < end ; ++i)
                            {
                                for (unsigned long blocki(0) ; blocki < blocksize ; ++blocki)
                                {
                                    sum += Ax[(i*blocksize)+blocki] * b[Aj[i] + blocki];
                                }
                            }
                            result[row] = sum;
                        }
                        break;
                }
            }
            else
            {
                __m128d A_ij, b_v;
                union sse4
                {
                    __m128d m;
                    double d[2];
                } sum_v;

                for (unsigned long row(row_start) ; row < row_end ; ++row)
                {
                    sum_v.m = _mm_setzero_pd();

                    const unsigned long end(Ar[row+1]);
                    for (unsigned long i(Ar[row]) ; i < end ; i++)
                    {
                        A_ij = _mm_load_pd(Ax + (i * blocksize));
                        const unsigned long col = Aj[i];
                        b_v = _mm_loadu_pd(b+col);
                        b_v = _mm_mul_pd(A_ij, b_v);
                        sum_v.m = _mm_add_pd(b_v, sum_v.m);
                    }
                    double sum = sum_v.d[0] + sum_v.d[1];

                    result[row] = sum;
                }
            }
        }
    }
}

