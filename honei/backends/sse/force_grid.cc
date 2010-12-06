/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
#include <honei/backends/sse/sse_mathfun.hh>
#include <cmath>
#include <limits>
#include <iostream>

#include <xmmintrin.h>
#include <emmintrin.h>

namespace honei
{
    namespace sse
    {
        void force_friction(const unsigned long start, unsigned long end, float * f_temp, const float * h,
                const float * u, const float * v, float force_multiplier, float dist, float manning_const_sq)
        {
            float divisor(float(1.)/float(3.));
            float factor(force_multiplier * dist * manning_const_sq);
            __m128 factorv = _mm_set1_ps(factor);
            float HONEI_ALIGNED(16) inf(std::numeric_limits<float>::infinity());
            __m128 infv = _mm_load1_ps(&inf);

            unsigned long size(end - start);

            unsigned long x_address((unsigned long)&f_temp[start]);
            unsigned long x_offset(x_address % 16);

            unsigned long z_offset(x_offset / 4);
            z_offset = (4 - z_offset) % 4;

            unsigned long quad_start(z_offset + start);
            unsigned long quad_end(end - ((end - quad_start) % 4));

            if (size < 16)
            {
                quad_end = start;
                quad_start = start;
            }

            for (unsigned long index(start) ; index < quad_start ; ++index)
            {
                float powh(pow(h[index], divisor));
                if ( (powh) > std::numeric_limits<float>::epsilon() || (powh) < float(-std::numeric_limits<float>::epsilon()) )
                {
                    f_temp[index] -= factor * u[index] * sqrt(u[index] * u[index] + v[index] * v[index]) / (powh);
                }
            }
            for (unsigned long index(quad_end) ; index < end ; ++index)
            {
                float powh(pow(h[index], divisor));
                if ( (powh) > std::numeric_limits<float>::epsilon() || (powh) < float(-std::numeric_limits<float>::epsilon()) )
                {
                    f_temp[index] -= factor * u[index] * sqrt(u[index] * u[index] + v[index] * v[index]) / (powh);
                }
            }

            __m128 m1, m2, m3, m4;
            for (unsigned long index(quad_start) ; index < quad_end ; index += 4)
            {
                // calc cubic root of h by exp(ln(h) / 3)
                m1 = _mm_load_ps(h + index);
                m1 = log_ps(m1);
                m2 = _mm_set1_ps(float(3.0));
                m1 = _mm_div_ps(m1, m2);
                m1 = exp_ps(m1);

                // sqrt(u[index] * u[index] + v[index] * v[index])
                m4 = _mm_load_ps(u + index);
                m2 = _mm_mul_ps(m4, m4);
                m3 = _mm_load_ps(v + index);
                m3 = _mm_mul_ps(m3, m3);
                m2 = _mm_add_ps(m2, m3);
                m2 = _mm_sqrt_ps(m2);

                // factor * u * sqrt() / powh
                m2 = _mm_mul_ps(m2, m4);
                m2 = _mm_mul_ps(factorv, m2);
                m2 = _mm_div_ps(m2, m1);

                // cmp with infinity
                // \todo Capture negative values of h resulting in negative infinity...
                m1 = _mm_cmpneq_ps(m2, infv);
                m2 = _mm_and_ps(m2, m1);
                //m1 = _mm_cmpnle_ps(m1, epsv);
                //m2 = _mm_and_ps(m2, m1);
                //m1 = _mm_cmple_ps(m2, mepsv);
                //m2 = _mm_and_ps(m2, m1);

                m1 = _mm_load_ps(f_temp + index);
                m1 = _mm_sub_ps(m1, m2);
                _mm_store_ps(f_temp + index, m1);
            }
        }

        void force_friction(const unsigned long start, const unsigned long end, double * f_temp, const double * h,
                const double * u, const double * v, double force_multiplier, double dist, double manning_const_sq)
        {
            double divisor(double(1.)/double(3.));
            double factor(force_multiplier * dist * manning_const_sq);
            //__m128d factorv = _mm_set1_pd(factor);
            //double HONEI_ALIGNED(16) inf(std::numeric_limits<double>::infinity());
            //__m128d infv = _mm_load1_pd(&inf);
            double HONEI_ALIGNED(16) eps(std::numeric_limits<double>::epsilon());
            __m128d epsv = _mm_load1_pd(&eps);

            unsigned long size(end - start);

            unsigned long x_address((unsigned long)&f_temp[start]);
            unsigned long x_offset(x_address % 16);

            unsigned long z_offset(x_offset / 8);

            unsigned long quad_start(z_offset + start);
            unsigned long quad_end(end - ((end - quad_start) % 2));

            if (size < 16)
            {
                quad_end = start;
                quad_start = start;
            }

            for (unsigned long index(start) ; index < quad_start ; ++index)
            {
                double powh(pow(h[index], divisor));
                if ( (powh) > std::numeric_limits<double>::epsilon() || (powh) < double(-std::numeric_limits<double>::epsilon()) )
                {
                    f_temp[index] -= factor * u[index] * sqrt(u[index] * u[index] + v[index] * v[index]) / (powh);
                }
            }
            for (unsigned long index(quad_end) ; index < end ; ++index)
            {
                double powh(pow(h[index], divisor));
                if ( (powh) > std::numeric_limits<double>::epsilon() || (powh) < double(-std::numeric_limits<double>::epsilon()) )
                {
                    f_temp[index] -= factor * u[index] * sqrt(u[index] * u[index] + v[index] * v[index]) / (powh);
                }
            }

            __m128d m1, m2, m3/*, m4*/;
            for (unsigned long index(quad_start) ; index < quad_end ; index += 2)
            {
                /*    // calc cubic root of h by exp(ln(h) / 3)
                      m1 = _mm_load_pd(h + index);
                      m1 = log_pd(m1);
                      m2 = _mm_set1_pd(double(3.0));
                      m1 = _mm_div_pd(m1, m2);
                      m1 = exp_pd(m1);

                // sqrt(u[index] * u[index] + v[index] * v[index])
                m4 = _mm_load_pd(u + index);
                m2 = _mm_mul_pd(m4, m4);
                m3 = _mm_load_pd(v + index);
                m3 = _mm_mul_pd(m3, m3);
                m2 = _mm_add_pd(m2, m3);
                m2 = _mm_sqrt_pd(m2);

                // factor * u * sqrt() / powh
                m2 = _mm_mul_pd(m2, m4);
                m2 = _mm_mul_pd(factorv, m2);
                m2 = _mm_div_pd(m2, m1);
                */
                double HONEI_ALIGNED(16) temp_pow[2];
                temp_pow[0]=(pow(h[index], divisor));
                temp_pow[1]=(pow(h[index + 1], divisor));
                double HONEI_ALIGNED(16) temp[2];
                temp[0] = (factor * u[index] * sqrt(u[index] * u[index] + v[index] * v[index]));
                temp[1] = (factor * u[index + 1] * sqrt(u[index + 1] * u[index + 1] + v[index + 1] * v[index + 1]));
                m1 = _mm_load_pd(temp_pow);
                m2 = _mm_load_pd(temp);
                m2 = _mm_div_pd(m2, m1);


                // cmp with infinity
                // \todo Capture negative values of h resulting in negative infinity...
                //m1 = _mm_cmpneq_pd(m2, infv);
                //m2 = _mm_and_pd(m2, m1);

                m3 = _mm_load_pd(h + index);
                m1 = _mm_cmpnle_pd(m3, epsv);
                m2 = _mm_and_pd(m2, m1);
                //m1 = _mm_cmple_pd(m2, mepsv);
                //m2 = _mm_and_pd(m2, m1);

                m1 = _mm_load_pd(f_temp + index);
                m1 = _mm_sub_pd(m1, m2);
                _mm_store_pd(f_temp + index, m1);
            }
        }
    }
}

