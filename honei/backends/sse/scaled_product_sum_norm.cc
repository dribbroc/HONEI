/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2009 Markus Geveler <apryde@gmx.de>
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
#include <iostream>

namespace honei {

    namespace sse {

        float scaled_product_sum_norm(unsigned long size, float a, float * y, float b, float * A_x)
        {

            union sse4
            {
                __m128 m;
                float f[4];
            } m1, m2, m3, m4, m5;

            float result(0);

            //compute ||ay + bA_x|| (SSE)
            m1.m = _mm_set_ps1(a);
            m2.m = _mm_set_ps1(b);

            m5.m = _mm_setzero_ps();

            unsigned long quad_end(size - (size % 4));
            for (unsigned long i(0) ; i < quad_end ; i += 4)
            {
                m3.m = _mm_load_ps(A_x + i);
                m4.m = _mm_load_ps(y + i);

                m4.m = _mm_mul_ps(m1.m, m4.m);
                m3.m = _mm_mul_ps(m2.m, m3.m);

                m3.m = _mm_add_ps(m3.m, m4.m);
                m3.m = _mm_mul_ps(m3.m, m3.m);

                m5.m = _mm_add_ps(m5.m, m3.m);
            }

            result += m5.f[0];
            result += m5.f[1];
            result += m5.f[2];
            result += m5.f[3];

            //compute ||ay + bA_x|| (FPU)
            for (unsigned long i(quad_end) ; i < size ; ++i)
            {
                result += (a * y[i] + b * A_x[i]) * (a * y[i] + b * A_x[i]);
            }

            return result;
        }

        double scaled_product_sum_norm(unsigned long size, double a, double * y, double b, double * A_x)
        {

            union sse2
            {
                __m128d m;
                double f[2];
            } m1, m2, m3, m4, m5;

            double result(0);

            //compute ||ay + bA_x|| (SSE)
            m1.m = _mm_set_pd1(a);
            m2.m = _mm_set_pd1(b);

            m5.m = _mm_setzero_pd();

            unsigned long quad_end(size - (size % 2));
            for (unsigned long i(0) ; i < quad_end ; i += 2)
            {
                m3.m = _mm_load_pd(A_x + i);
                m4.m = _mm_load_pd(y + i);

                m4.m = _mm_mul_pd(m1.m, m4.m);
                m3.m = _mm_mul_pd(m2.m, m3.m);

                m3.m = _mm_add_pd(m3.m, m4.m);
                m3.m = _mm_mul_pd(m3.m, m3.m);

                m5.m = _mm_add_pd(m5.m, m3.m);
            }

            result += m5.f[0];
            result += m5.f[1];

            //compute ||ay + bA_x|| (FPU)
            for (unsigned long i(quad_end) ; i < size ; ++i)
            {
                result += (a * y[i] + b * A_x[i]) * (a * y[i] + b * A_x[i]);
            }

            return result;
        }

        float scaled_product_sum_norm(unsigned long size, unsigned long m, float a, float * y, float b, float * ll, float * ld, float * lu, float * dl, float * dd, float * du, float * ul, float * ud, float * uu, float * x)
        {

            __m128 m1, m2, m3, m4, m5, m6, m7;
            m4 = _mm_set_ps1(a);
            m5 = _mm_set_ps1(b);

            union sse4
            {
                __m128 m;
                float f[4];
            } m8;

            float result(0);
            float temp(0);

            unsigned long index(0);

            //index 0
            temp = (a * y[index] + b * (dd[index] * x[index]
                + du[index] * x[index + 1]
                + ul[index] * x[index + m - 1]
                + ud[index] * x[index + m]
                + uu[index] * x[index + m + 1]));
            result += temp * temp;

            //index in [1, root_n -1[
            index = 1;
            unsigned long quad_start(index + 4 - (index % 4));
            unsigned long quad_end(m - 1 - ((m - 1 - quad_start) % 4));

            for (unsigned long si(index) ; si < quad_start ; ++si)
            {
                temp = (a * y[si] + b * (dd[si] * x[si]
                    + dl[si] * x[si - 1]
                    + du[si] * x[si + 1]
                    + ul[si] * x[si + m - 1]
                    + ud[si] * x[si + m]
                    + uu[si] * x[si + m + 1]));
                result += temp * temp;
            }

            for (unsigned long si(quad_start) ; si < quad_end ; si+=4)
            {
                m6 = _mm_load_ps(y + si);
                m6 = _mm_mul_ps(m6, m4);

                //result[si] = dd[si] * b[si]
                m2 = _mm_loadu_ps(x + si);
                m3 = _mm_load_ps(dd + si);
                m1 = _mm_mul_ps(m3, m2);

                //+ dl[si] * b[si - 1]
                m2 = _mm_loadu_ps(x + si - 1);
                m3 = _mm_load_ps(dl + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ du[si] * b[si + 1]
                m2 = _mm_loadu_ps(x + si + 1);
                m3 = _mm_load_ps(du + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ ul[si] * b[si + m - 1]
                m2 = _mm_loadu_ps(x + si + m - 1);
                m3 = _mm_load_ps(ul + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ ud[si] * b[si + m]
                m2 = _mm_loadu_ps(x + si + m);
                m3 = _mm_load_ps(ud + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ uu[si] * b[si + m + 1];
                m2 = _mm_loadu_ps(x + si + m + 1);
                m3 = _mm_load_ps(uu + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                m1 = _mm_mul_ps(m1, m5);
                m1 = _mm_add_ps(m1, m6);
                m1 = _mm_mul_ps(m1, m1);

                m8.m = _mm_add_ps(m8.m, m1);
            }

            for (unsigned long si(quad_end) ; si < m - 1 ; ++si)
            {
                temp = (a * y[si] + b * (dd[si] * x[si]
                    + dl[si] * x[si - 1]
                    + du[si] * x[si + 1]
                    + ul[si] * x[si + m - 1]
                    + ud[si] * x[si + m]
                    + uu[si] * x[si + m + 1]));
                result += temp * temp;
            }

            //index root_n -1
            index = m - 1;
            temp = (a * y[index] + b * (dd[index] * x[index]
                + lu[index] * x[index - m + 1]
                + dl[index] * x[index - 1]
                + du[index] * x[index + 1]
                + ul[index] * x[index + m - 1]
                + ud[index] * x[index + m]
                + uu[index] * x[index + m + 1]));
            result += temp * temp;

            //index root_n
            index = m;
            temp = (a * y[index] + b * (dd[index] * x[index]
                + lu[index] * x[index - m + 1]
                + ld[index] * x[index - m]
                + dl[index] * x[index - 1]
                + du[index] * x[index + 1]
                + ul[index] * x[index + m - 1]
                + ud[index] * x[index + m]
                + uu[index] * x[index + m + 1]));
            result += temp * temp;

            //index in [root_n + 1, n - (root_n + 1)[
            index = m + 1;
            quad_start = index + 4 - (index % 4);
            quad_end = size - m - 1 - ((size - m - 1 - quad_start) % 4);

            for (unsigned long si(index) ; si < quad_start ; ++si)
            {
                temp = (a * y[si] + b * (dd[si] * x[si]
                    + ll[si] * x[si - m - 1]
                    + lu[si] * x[si - m + 1]
                    + ld[si] * x[si - m]
                    + dl[si] * x[si - 1]
                    + du[si] * x[si + 1]
                    + ul[si] * x[si + m - 1]
                    + ud[si] * x[si + m]
                    + uu[si] * x[si + m + 1]));
                result += temp * temp;
            }

            for (unsigned long si(quad_start) ; si < quad_end ; si+=4)
            {
                m6 = _mm_load_ps(y + si);
                m6 = _mm_mul_ps(m6, m4);

                //result[si] = dd[si] * b[si]
                m2 = _mm_loadu_ps(x + si);
                m3 = _mm_load_ps(dd + si);
                m1 = _mm_mul_ps(m3, m2);

                //+ ll[si] * b[si - m - 1]
                m2 = _mm_loadu_ps(x + si -m - 1);
                m3 = _mm_load_ps(ll + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ lu[si] * b[si - m + 1]
                m2 = _mm_loadu_ps(x + si - m + 1);
                m3 = _mm_load_ps(lu + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ ld[si] * b[si - m]
                m2 = _mm_loadu_ps(x + si - m);
                m3 = _mm_load_ps(ld + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ dl[si] * b[si - 1]
                m2 = _mm_loadu_ps(x + si - 1);
                m3 = _mm_load_ps(dl + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ du[si] * b[si + 1]
                m2 = _mm_loadu_ps(x + si + 1);
                m3 = _mm_load_ps(du + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ ul[si] * b[si + m - 1]
                m2 = _mm_loadu_ps(x + si + m - 1);
                m3 = _mm_load_ps(ul + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ ud[si] * b[si + m]
                m2 = _mm_loadu_ps(x + si + m);
                m3 = _mm_load_ps(ud + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ uu[si] * b[si + m + 1];
                m2 = _mm_loadu_ps(x + si + m + 1);
                m3 = _mm_load_ps(uu + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                m1 = _mm_mul_ps(m1, m5);
                m1 = _mm_add_ps(m1, m6);
                m1 = _mm_mul_ps(m1, m1);

                m8.m = _mm_add_ps(m8.m, m1);
            }

            for (unsigned long si(quad_end) ; si < size - m - 1 ; ++si)
            {
                temp = (a * y[si] + b * (dd[si] * x[si]
                    + ll[si] * x[si - m - 1]
                    + lu[si] * x[si - m + 1]
                    + ld[si] * x[si - m]
                    + dl[si] * x[si - 1]
                    + du[si] * x[si + 1]
                    + ul[si] * x[si + m - 1]
                    + ud[si] * x[si + m]
                    + uu[si] * x[si + m + 1]));
                result += temp * temp;
            }

            //index n - (root_n + 1)
            index = size - m - 1;
            temp = (a * y[index] + b * (dd[index] * x[index]
                + ll[index] * x[index - m - 1]
                + lu[index] * x[index - m + 1]
                + ld[index] * x[index - m]
                + dl[index] * x[index - 1]
                + du[index] * x[index + 1]
                + ul[index] * x[index + m - 1]
                + ud[index] * x[index + m]));
            result += temp * temp;

            //index n - root_n
            index = size - m;
            temp = (a * y[index] + b * (dd[index] * x[index]
                + ll[index] * x[index - m - 1]
                + lu[index] * x[index - m + 1]
                + ld[index] * x[index - m]
                + dl[index] * x[index - 1]
                + du[index] * x[index + 1]
                + ul[index] * x[index + m - 1]));
            result += temp * temp;

            //index in [n - root_n + 1, n -1[
            index = size - m + 1;
            quad_start = index + 4 - (index % 4);
            quad_end = size - 1 - ((size - 1 - quad_start) % 4);

            for (unsigned long si(index) ; si < quad_start ; ++si)
            {
                temp = (a * y[si] + b * (dd[si] * x[si]
                    + ll[si] * x[si - m - 1]
                    + lu[si] * x[si - m + 1]
                    + ld[si] * x[si - m]
                    + dl[si] * x[si - 1]
                    + du[si] * x[si + 1]));
            }

            for (unsigned long si(quad_start) ; si < quad_end ; si+=4)
            {
                m6 = _mm_load_ps(y + si);
                m6 = _mm_mul_ps(m6, m4);

                //result[si] = dd[si] * b[si]
                m2 = _mm_loadu_ps(x + si);
                m3 = _mm_load_ps(dd + si);
                m1 = _mm_mul_ps(m3, m2);

                //+ ll[si] * b[si - m - 1]
                m2 = _mm_loadu_ps(x + si -m - 1);
                m3 = _mm_load_ps(ll + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ lu[si] * b[si - m + 1]
                m2 = _mm_loadu_ps(x + si - m + 1);
                m3 = _mm_load_ps(lu + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ ld[si] * b[si - m]
                m2 = _mm_loadu_ps(x + si - m);
                m3 = _mm_load_ps(ld + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ dl[si] * b[si - 1]
                m2 = _mm_loadu_ps(x + si - 1);
                m3 = _mm_load_ps(dl + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                //+ du[si] * b[si + 1]
                m2 = _mm_loadu_ps(x + si + 1);
                m3 = _mm_load_ps(du + si);
                m3 = _mm_mul_ps(m3, m2);
                m1 = _mm_add_ps(m1, m3);

                m1 = _mm_mul_ps(m1, m5);
                m1 = _mm_add_ps(m1, m6);
                m1 = _mm_mul_ps(m1, m1);

                m8.m = _mm_add_ps(m8.m, m1);
            }

            for (unsigned long si(quad_end) ; si < size - m - 1 ; ++si)
            {
                temp = (a * y[si] + b * (dd[si] * x[si]
                            + ll[si] * x[si - m - 1]
                            + lu[si] * x[si - m + 1]
                            + ld[si] * x[si - m]
                            + dl[si] * x[si - 1]
                            + du[si] * x[si + 1]));
                result += temp * temp;
            }

            //index n - 1
            index = size - 1;
            temp = (a * y[index] + b * (dd[index] * x[index]
                + ll[index] * x[index - m - 1]
                + lu[index] * x[index - m + 1]
                + ld[index] * x[index - m]
                + dl[index] * x[index - 1]));
            result += temp * temp;

            result += m8.f[0];
            result += m8.f[1];
            result += m8.f[2];
            result += m8.f[3];

            return result;
        }
/*
        double scaled_product_sum_norm(unsigned long size, double a, double * y, double b, double * ll, double * ld, double * lu, double * dl, double * dd, double * du, double * ul, double * ud, double * uu, double * x)
        {

            union sse2
            {
                __m128d m;
                double f[2];
            } m1, m2, m3, m4, m5;

            double result(0);

            return result;
        }*/
    }
}
