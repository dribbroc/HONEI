/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2009 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
        void eq_dist_grid_dir_0(unsigned long begin, unsigned long end,
                float g, float e,
                float * h, float * u, float * v,
                float * f_eq_0)
        {
            unsigned long size(end - begin);

            unsigned long x_address((unsigned long)&h[begin]);
            unsigned long x_offset(x_address % 16);

            unsigned long z_offset(x_offset / 4);
            z_offset = (4 - z_offset) % 4;

            unsigned long quad_start(z_offset + begin);
            unsigned long quad_end(end - ((end - quad_start) % 4));

            if (size < 32)
            {
                quad_end = begin;
                quad_start = begin;
            }

            float e2(e);
            float e23(float(3.) * e2);
            float e26(float(6.) * e2);
            for (unsigned long index(begin) ; index < quad_start ; ++index)
            {
                float u2(u[index] * u[index]);
                float v2(v[index] * v[index]);
                float h2(h[index] * h[index]);

                float t1, t2;

                t1 = (float(5.) * g * h2) / e26;
                t2 = (float(2.) * h[index]) / e23 * (u2 + v2);
                f_eq_0[index] = h[index] - t1 - t2;
            }
            for (unsigned long index(quad_end) ; index < end ; ++index)
            {
                float u2(u[index] * u[index]);
                float v2(v[index] * v[index]);
                float h2(h[index] * h[index]);

                float t1, t2;

                t1 = (float(5.) * g * h2) / e26;
                t2 = (float(2.) * h[index]) / e23 * (u2 + v2);
                f_eq_0[index] = h[index] - t1 - t2;
            }

            __m128 gv = _mm_set1_ps(g);
            __m128 e2v = _mm_set1_ps(e);
            __m128 t1, t2, m1, m2;
            __m128 scal1, scal2;

            // f_eq_0
            for (unsigned long index(quad_start) ; index < quad_end ; index += 4)
            {

                scal1 = _mm_set1_ps(float(5));
                scal2 = _mm_set1_ps(float(6));
                m1 = _mm_mul_ps(scal1, gv);
                m2 = _mm_load_ps(h + index);
                m2 = _mm_mul_ps(m2, m2);
                m1 = _mm_mul_ps(m1, m2);
                m2 = _mm_mul_ps(e2v, scal2);
                t1 = _mm_div_ps(m1, m2);

                scal1 = _mm_set1_ps(float(2));
                scal2 = _mm_set1_ps(float(3));
                m1 = _mm_load_ps(h + index);
                m1 = _mm_mul_ps(scal1, m1);
                m2 = _mm_mul_ps(e2v, scal2);
                t2 = _mm_div_ps(m1, m2);
                m1 = _mm_load_ps(u + index);
                m1 = _mm_mul_ps(m1, m1);
                m2 = _mm_load_ps(v + index);
                m2 = _mm_mul_ps(m2, m2);
                m1 = _mm_add_ps(m1, m2);
                t2 = _mm_mul_ps(t2, m1);

                m1 = _mm_load_ps(h + index);
                m1 = _mm_sub_ps(m1, t1);
                m1 = _mm_sub_ps(m1, t2);
                _mm_store_ps(f_eq_0 + index, m1);
            }
        }

        void eq_dist_grid_dir_odd(unsigned long begin, unsigned long end,
                float g, float e,
                float * h, float * u, float * v,
                float * distribution_x, float * distribution_y,
                float * f_eq,
                unsigned long dir)
        {
            unsigned long size(end - begin);

            unsigned long x_address((unsigned long)&h[begin]);
            unsigned long x_offset(x_address % 16);

            unsigned long z_offset(x_offset / 4);
            z_offset = (4 - z_offset) % 4;

            unsigned long quad_start(z_offset + begin);
            unsigned long quad_end(end - ((end - quad_start) % 4));

            if (size < 32)
            {
                quad_end = begin;
                quad_start = begin;
            }

            float e2(e);
            float e42(float(2.) * e2 * e2);
            float e23(float(3.) * e2);
            float e26(float(6.) * e2);
            for (unsigned long index(begin) ; index < quad_start ; ++index)
            {
                float u2(u[index] * u[index]);
                float v2(v[index] * v[index]);
                float h2(h[index] * h[index]);

                float dxu, dyv;
                float t1, t2, t3, t4;
                dxu = distribution_x[dir] * u[index];
                dyv = distribution_y[dir] * v[index];
                t1 = (g * h2) / e26;
                t2 = (h[index] / e23) * (dxu + dyv);
                t3 = (h[index] / e42) * (dxu * dxu + float(2.) * dxu * dyv + dyv * dyv);
                t4 = (h[index] / e26) * (u2 + v2);
                f_eq[index] = t1 + t2 + t3 - t4;
            }
            for (unsigned long index(quad_end) ; index < end ; ++index)
            {
                float u2(u[index] * u[index]);
                float v2(v[index] * v[index]);
                float h2(h[index] * h[index]);

                float dxu, dyv;
                float t1, t2, t3, t4;
                dxu = distribution_x[dir] * u[index];
                dyv = distribution_y[dir] * v[index];
                t1 = (g * h2) / e26;
                t2 = (h[index] / e23) * (dxu + dyv);
                t3 = (h[index] / e42) * (dxu * dxu + float(2.) * dxu * dyv + dyv * dyv);
                t4 = (h[index] / e26) * (u2 + v2);
                f_eq[index] = t1 + t2 + t3 - t4;
            }

            __m128 gv = _mm_set1_ps(g);
            __m128 e2v = _mm_set1_ps(e);
            __m128 t1, t2, t3, t4, m1, m2, dxu, dyv;
            __m128 scal1;

            // f_eq_odd
            for (unsigned long index(quad_start) ; index < quad_end ; index += 4)
            {
                m1 = _mm_set1_ps(distribution_x[dir]);
                m2 = _mm_load_ps(u + index);
                dxu = _mm_mul_ps(m1, m2);
                m1 = _mm_set1_ps(distribution_y[dir]);
                m2 = _mm_load_ps(v + index);
                dyv = _mm_mul_ps(m1, m2);

                scal1 = _mm_set1_ps(float(6));
                m1 = _mm_load_ps(h + index);
                m1 = _mm_mul_ps(m1, m1);
                m1 = _mm_mul_ps(m1, gv);
                scal1 = _mm_mul_ps(e2v, scal1);
                t1 = _mm_div_ps(m1, scal1);

                m1 = _mm_load_ps(h + index);
                t4 = _mm_div_ps(m1, scal1);
                m1 = _mm_load_ps(u + index);
                m1 = _mm_mul_ps(m1, m1);
                m2 = _mm_load_ps(v + index);
                m2 = _mm_mul_ps(m2, m2);
                m1 = _mm_add_ps(m1, m2);
                t4 = _mm_mul_ps(t4, m1);

                scal1 = _mm_set1_ps(float(3));
                m1 = _mm_load_ps(h + index);
                m2 = _mm_mul_ps(e2v, scal1);
                m1 = _mm_div_ps(m1, m2);
                m2 = _mm_add_ps(dxu, dyv);
                t2 = _mm_mul_ps(m1, m2);

                scal1 = _mm_set1_ps(float(2));
                m1 = _mm_load_ps(h + index);
                m2 = _mm_mul_ps(e2v, e2v);
                m2 = _mm_mul_ps(m2, scal1);
                t3 = _mm_div_ps(m1, m2);

                m1 = _mm_mul_ps(dxu, dxu);
                m2 = _mm_mul_ps(scal1, dxu);
                m2 = _mm_mul_ps(m2, dyv);
                m1 = _mm_add_ps(m1, m2);
                m2 = _mm_mul_ps(dyv, dyv);
                m1 = _mm_add_ps(m1, m2);

                t3 = _mm_mul_ps(m1, t3);

                m1 = _mm_add_ps(t1, t2);
                m1 = _mm_add_ps(m1, t3);
                m1 = _mm_sub_ps(m1, t4);
                _mm_store_ps(f_eq + index, m1);
            }
        }

        void eq_dist_grid_dir_even(unsigned long begin, unsigned long end,
                float g, float e,
                float * h, float * u, float * v,
                float * distribution_x, float * distribution_y,
                float * f_eq,
                unsigned long dir)
        {
            unsigned long size(end - begin);

            unsigned long x_address((unsigned long)&h[begin]);
            unsigned long x_offset(x_address % 16);

            unsigned long z_offset(x_offset / 4);
            z_offset = (4 - z_offset) % 4;

            unsigned long quad_start(z_offset + begin);
            unsigned long quad_end(end - ((end - quad_start) % 4));

            if (size < 32)
            {
                quad_end = begin;
                quad_start = begin;
            }

            float e2(e);
            float e48(float(8.) * e2 * e2);
            float e212(float(12.) * e2);
            float e224(float(24.) * e2);
            for (unsigned long index(begin) ; index < quad_start ; ++index)
            {
                float u2(u[index] * u[index]);
                float v2(v[index] * v[index]);
                float h2(h[index] * h[index]);

                float dxu, dyv;
                float t1, t2, t3, t4;
                dxu = distribution_x[dir] * u[index];
                dyv = distribution_y[dir] * v[index];
                t1 = (g * h2) / e224;
                t2 = (h[index] / e212) * (dxu + dyv);
                t3 = (h[index] / e48) * (dxu * dxu + float(2.) * dxu * dyv + dyv * dyv);
                t4 = (h[index] / e224) * (u2 + v2);
                f_eq[index] =  t1 + t2 + t3 - t4;
            }
            for (unsigned long index(quad_end) ; index < end ; ++index)
            {
                float u2(u[index] * u[index]);
                float v2(v[index] * v[index]);
                float h2(h[index] * h[index]);

                float dxu, dyv;
                float t1, t2, t3, t4;
                dxu = distribution_x[dir] * u[index];
                dyv = distribution_y[dir] * v[index];
                t1 = (g * h2) / e224;
                t2 = (h[index] / e212) * (dxu + dyv);
                t3 = (h[index] / e48) * (dxu * dxu + float(2.) * dxu * dyv + dyv * dyv);
                t4 = (h[index] / e224) * (u2 + v2);
                f_eq[index] =  t1 + t2 + t3 - t4;
            }

            __m128 gv = _mm_set1_ps(g);
            __m128 e2v = _mm_set1_ps(e);
            __m128 t1, t2, t3, t4, m1, m2, dxu, dyv;
            __m128 scal1;

            // f_eq_even
            for (unsigned long index(quad_start) ; index < quad_end ; index += 4)
            {
                m1 = _mm_set1_ps(distribution_x[dir]);
                m2 = _mm_load_ps(u + index);
                dxu = _mm_mul_ps(m1, m2);
                m1 = _mm_set1_ps(distribution_y[dir]);
                m2 = _mm_load_ps(v + index);
                dyv = _mm_mul_ps(m1, m2);

                scal1 = _mm_set1_ps(float(24));
                m1 = _mm_load_ps(h + index);
                m1 = _mm_mul_ps(m1, m1);
                m1 = _mm_mul_ps(m1, gv);
                scal1 = _mm_mul_ps(e2v, scal1);
                t1 = _mm_div_ps(m1, scal1);

                m1 = _mm_load_ps(h + index);
                t4 = _mm_div_ps(m1, scal1);
                m1 = _mm_load_ps(u + index);
                m1 = _mm_mul_ps(m1, m1);
                m2 = _mm_load_ps(v + index);
                m2 = _mm_mul_ps(m2, m2);
                m1 = _mm_add_ps(m1, m2);
                t4 = _mm_mul_ps(t4, m1);

                scal1 = _mm_set1_ps(float(12));
                m1 = _mm_load_ps(h + index);
                m2 = _mm_mul_ps(e2v, scal1);
                m1 = _mm_div_ps(m1, m2);
                m2 = _mm_add_ps(dxu, dyv);
                t2 = _mm_mul_ps(m1, m2);

                scal1 = _mm_set1_ps(float(8));
                m1 = _mm_load_ps(h + index);
                m2 = _mm_mul_ps(e2v, e2v);
                m2 = _mm_mul_ps(m2, scal1);
                t3 = _mm_div_ps(m1, m2);

                m1 = _mm_mul_ps(dxu, dxu);
                scal1 = _mm_set1_ps(float(2));
                m2 = _mm_mul_ps(scal1, dxu);
                m2 = _mm_mul_ps(m2, dyv);
                m1 = _mm_add_ps(m1, m2);
                m2 = _mm_mul_ps(dyv, dyv);
                m1 = _mm_add_ps(m1, m2);

                t3 = _mm_mul_ps(m1, t3);

                m1 = _mm_add_ps(t1, t2);
                m1 = _mm_add_ps(m1, t3);
                m1 = _mm_sub_ps(m1, t4);
                _mm_store_ps(f_eq + index, m1);
            }
        }

        void eq_dist_grid_dir_0(unsigned long begin, unsigned long end,
                double g, double e,
                double * h, double * u, double * v,
                double * f_eq_0)
        {
            unsigned long size(end - begin);

            unsigned long x_address((unsigned long)&h[begin]);
            unsigned long x_offset(x_address % 16);

            unsigned long z_offset(x_offset / 8);

            unsigned long quad_start(z_offset + begin);
            unsigned long quad_end(end - ((end - quad_start) % 2));

            if (size < 32)
            {
                quad_end = begin;
                quad_start = begin;
            }

            double e2(e);
            double e23(double(3.) * e2);
            double e26(double(6.) * e2);
            for (unsigned long index(begin) ; index < quad_start ; ++index)
            {
                double u2(u[index] * u[index]);
                double v2(v[index] * v[index]);
                double h2(h[index] * h[index]);

                double t1, t2;

                t1 = (double(5.) * g * h2) / e26;
                t2 = (double(2.) * h[index]) / e23 * (u2 + v2);
                f_eq_0[index] = h[index] - t1 - t2;
            }
            for (unsigned long index(quad_end) ; index < end ; ++index)
            {
                double u2(u[index] * u[index]);
                double v2(v[index] * v[index]);
                double h2(h[index] * h[index]);

                double t1, t2;

                t1 = (double(5.) * g * h2) / e26;
                t2 = (double(2.) * h[index]) / e23 * (u2 + v2);
                f_eq_0[index] = h[index] - t1 - t2;
            }

            __m128d gv = _mm_set1_pd(g);
            __m128d e2v = _mm_set1_pd(e);
            __m128d t1, t2, m1, m2;
            __m128d scal1, scal2;

            // f_eq_0
            for (unsigned long index(quad_start) ; index < quad_end ; index += 2)
            {
                scal1 = _mm_set1_pd(double(5));
                scal2 = _mm_set1_pd(double(6));
                m1 = _mm_mul_pd(scal1, gv);
                m2 = _mm_load_pd(h + index);
                m2 = _mm_mul_pd(m2, m2);
                m1 = _mm_mul_pd(m1, m2);
                m2 = _mm_mul_pd(e2v, scal2);
                t1 = _mm_div_pd(m1, m2);

                scal1 = _mm_set1_pd(double(2));
                scal2 = _mm_set1_pd(double(3));
                m1 = _mm_load_pd(h + index);
                m1 = _mm_mul_pd(scal1, m1);
                m2 = _mm_mul_pd(e2v, scal2);
                t2 = _mm_div_pd(m1, m2);
                m1 = _mm_load_pd(u + index);
                m1 = _mm_mul_pd(m1, m1);
                m2 = _mm_load_pd(v + index);
                m2 = _mm_mul_pd(m2, m2);
                m1 = _mm_add_pd(m1, m2);
                t2 = _mm_mul_pd(t2, m1);

                m1 = _mm_load_pd(h + index);
                m1 = _mm_sub_pd(m1, t1);
                m1 = _mm_sub_pd(m1, t2);
                _mm_store_pd(f_eq_0 + index, m1);
            }
        }

        void eq_dist_grid_dir_odd(unsigned long begin, unsigned long end,
                double g, double e,
                double * h, double * u, double * v,
                double * distribution_x, double * distribution_y,
                double * f_eq,
                unsigned long dir)
        {
            unsigned long size(end - begin);

            unsigned long x_address((unsigned long)&h[begin]);
            unsigned long x_offset(x_address % 16);

            unsigned long z_offset(x_offset / 8);

            unsigned long quad_start(z_offset + begin);
            unsigned long quad_end(end - ((end - quad_start) % 2));

            if (size < 32)
            {
                quad_end = begin;
                quad_start = begin;
            }

            double e2(e);
            double e42(double(2.) * e2 * e2);
            double e23(double(3.) * e2);
            double e26(double(6.) * e2);
            for (unsigned long index(begin) ; index < quad_start ; ++index)
            {
                double u2(u[index] * u[index]);
                double v2(v[index] * v[index]);
                double h2(h[index] * h[index]);

                double dxu, dyv;
                double t1, t2, t3, t4;
                dxu = distribution_x[dir] * u[index];
                dyv = distribution_y[dir] * v[index];
                t1 = (g * h2) / e26;
                t2 = (h[index] / e23) * (dxu + dyv);
                t3 = (h[index] / e42) * (dxu * dxu + double(2.) * dxu * dyv + dyv * dyv);
                t4 = (h[index] / e26) * (u2 + v2);
                f_eq[index] = t1 + t2 + t3 - t4;
            }
            for (unsigned long index(quad_end) ; index < end ; ++index)
            {
                double u2(u[index] * u[index]);
                double v2(v[index] * v[index]);
                double h2(h[index] * h[index]);

                double dxu, dyv;
                double t1, t2, t3, t4;
                dxu = distribution_x[dir] * u[index];
                dyv = distribution_y[dir] * v[index];
                t1 = (g * h2) / e26;
                t2 = (h[index] / e23) * (dxu + dyv);
                t3 = (h[index] / e42) * (dxu * dxu + double(2.) * dxu * dyv + dyv * dyv);
                t4 = (h[index] / e26) * (u2 + v2);
                f_eq[index] = t1 + t2 + t3 - t4;
            }

            __m128d gv = _mm_set1_pd(g);
            __m128d e2v = _mm_set1_pd(e);
            __m128d t1, t2, t3, t4, m1, m2, dxu, dyv;
            __m128d scal1;

            // f_eq_odd
            for (unsigned long index(quad_start) ; index < quad_end ; index += 2)
            {
                m1 = _mm_set1_pd(distribution_x[dir]);
                m2 = _mm_load_pd(u + index);
                dxu = _mm_mul_pd(m1, m2);
                m1 = _mm_set1_pd(distribution_y[dir]);
                m2 = _mm_load_pd(v + index);
                dyv = _mm_mul_pd(m1, m2);

                scal1 = _mm_set1_pd(double(6));
                m1 = _mm_load_pd(h + index);
                m1 = _mm_mul_pd(m1, m1);
                m1 = _mm_mul_pd(m1, gv);
                scal1 = _mm_mul_pd(e2v, scal1);
                t1 = _mm_div_pd(m1, scal1);

                m1 = _mm_load_pd(h + index);
                t4 = _mm_div_pd(m1, scal1);
                m1 = _mm_load_pd(u + index);
                m1 = _mm_mul_pd(m1, m1);
                m2 = _mm_load_pd(v + index);
                m2 = _mm_mul_pd(m2, m2);
                m1 = _mm_add_pd(m1, m2);
                t4 = _mm_mul_pd(t4, m1);

                scal1 = _mm_set1_pd(double(3));
                m1 = _mm_load_pd(h + index);
                m2 = _mm_mul_pd(e2v, scal1);
                m1 = _mm_div_pd(m1, m2);
                m2 = _mm_add_pd(dxu, dyv);
                t2 = _mm_mul_pd(m1, m2);

                scal1 = _mm_set1_pd(double(2));
                m1 = _mm_load_pd(h + index);
                m2 = _mm_mul_pd(e2v, e2v);
                m2 = _mm_mul_pd(m2, scal1);
                t3 = _mm_div_pd(m1, m2);

                m1 = _mm_mul_pd(dxu, dxu);
                m2 = _mm_mul_pd(scal1, dxu);
                m2 = _mm_mul_pd(m2, dyv);
                m1 = _mm_add_pd(m1, m2);
                m2 = _mm_mul_pd(dyv, dyv);
                m1 = _mm_add_pd(m1, m2);

                t3 = _mm_mul_pd(m1, t3);

                m1 = _mm_add_pd(t1, t2);
                m1 = _mm_add_pd(m1, t3);
                m1 = _mm_sub_pd(m1, t4);
                _mm_store_pd(f_eq + index, m1);
            }
        }

        void eq_dist_grid_dir_even(unsigned long begin, unsigned long end,
                double g, double e,
                double * h, double * u, double * v,
                double * distribution_x, double * distribution_y,
                double * f_eq,
                unsigned long dir)
        {
            unsigned long size(end - begin);

            unsigned long x_address((unsigned long)&h[begin]);
            unsigned long x_offset(x_address % 16);

            unsigned long z_offset(x_offset / 8);

            unsigned long quad_start(z_offset + begin);
            unsigned long quad_end(end - ((end - quad_start) % 2));

            if (size < 32)
            {
                quad_end = begin;
                quad_start = begin;
            }

            double e2(e);
            double e48(double(8.) * e2 * e2);
            double e212(double(12.) * e2);
            double e224(double(24.) * e2);
            for (unsigned long index(begin) ; index < quad_start ; ++index)
            {
                double u2(u[index] * u[index]);
                double v2(v[index] * v[index]);
                double h2(h[index] * h[index]);

                double dxu, dyv;
                double t1, t2, t3, t4;
                dxu = distribution_x[dir] * u[index];
                dyv = distribution_y[dir] * v[index];
                t1 = (g * h2) / e224;
                t2 = (h[index] / e212) * (dxu + dyv);
                t3 = (h[index] / e48) * (dxu * dxu + double(2.) * dxu * dyv + dyv * dyv);
                t4 = (h[index] / e224) * (u2 + v2);
                f_eq[index] =  t1 + t2 + t3 - t4;
            }
            for (unsigned long index(quad_end) ; index < end ; ++index)
            {
                double u2(u[index] * u[index]);
                double v2(v[index] * v[index]);
                double h2(h[index] * h[index]);

                double dxu, dyv;
                double t1, t2, t3, t4;
                dxu = distribution_x[dir] * u[index];
                dyv = distribution_y[dir] * v[index];
                t1 = (g * h2) / e224;
                t2 = (h[index] / e212) * (dxu + dyv);
                t3 = (h[index] / e48) * (dxu * dxu + double(2.) * dxu * dyv + dyv * dyv);
                t4 = (h[index] / e224) * (u2 + v2);
                f_eq[index] =  t1 + t2 + t3 - t4;
            }

            __m128d gv = _mm_set1_pd(g);
            __m128d e2v = _mm_set1_pd(e);
            __m128d t1, t2, t3, t4, m1, m2, dxu, dyv;
            __m128d scal1;

            // f_eq_even
            for (unsigned long index(quad_start) ; index < quad_end ; index += 2)
            {
                m1 = _mm_set1_pd(distribution_x[dir]);
                m2 = _mm_load_pd(u + index);
                dxu = _mm_mul_pd(m1, m2);
                m1 = _mm_set1_pd(distribution_y[dir]);
                m2 = _mm_load_pd(v + index);
                dyv = _mm_mul_pd(m1, m2);

                scal1 = _mm_set1_pd(double(24));
                m1 = _mm_load_pd(h + index);
                m1 = _mm_mul_pd(m1, m1);
                m1 = _mm_mul_pd(m1, gv);
                scal1 = _mm_mul_pd(e2v, scal1);
                t1 = _mm_div_pd(m1, scal1);

                m1 = _mm_load_pd(h + index);
                t4 = _mm_div_pd(m1, scal1);
                m1 = _mm_load_pd(u + index);
                m1 = _mm_mul_pd(m1, m1);
                m2 = _mm_load_pd(v + index);
                m2 = _mm_mul_pd(m2, m2);
                m1 = _mm_add_pd(m1, m2);
                t4 = _mm_mul_pd(t4, m1);

                scal1 = _mm_set1_pd(double(12));
                m1 = _mm_load_pd(h + index);
                m2 = _mm_mul_pd(e2v, scal1);
                m1 = _mm_div_pd(m1, m2);
                m2 = _mm_add_pd(dxu, dyv);
                t2 = _mm_mul_pd(m1, m2);

                scal1 = _mm_set1_pd(double(8));
                m1 = _mm_load_pd(h + index);
                m2 = _mm_mul_pd(e2v, e2v);
                m2 = _mm_mul_pd(m2, scal1);
                t3 = _mm_div_pd(m1, m2);

                m1 = _mm_mul_pd(dxu, dxu);
                scal1 = _mm_set1_pd(double(2));
                m2 = _mm_mul_pd(scal1, dxu);
                m2 = _mm_mul_pd(m2, dyv);
                m1 = _mm_add_pd(m1, m2);
                m2 = _mm_mul_pd(dyv, dyv);
                m1 = _mm_add_pd(m1, m2);

                t3 = _mm_mul_pd(m1, t3);

                m1 = _mm_add_pd(t1, t2);
                m1 = _mm_add_pd(m1, t3);
                m1 = _mm_sub_pd(m1, t4);
                _mm_store_pd(f_eq + index, m1);
            }
        }
    }
}
