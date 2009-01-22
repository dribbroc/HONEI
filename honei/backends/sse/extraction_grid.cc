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
        void extraction_grid(unsigned long begin, unsigned long end,
                float * distribution_x, float * distribution_y,
                float * h, float * u, float * v,
                float * f_0, float * f_1, float * f_2,
                float * f_3, float * f_4, float * f_5,
                float * f_6, float * f_7, float * f_8)
        {
            //unsigned long size(end - begin);

            unsigned long x_address((unsigned long)&h[begin]);
            unsigned long x_offset(x_address % 16);

            unsigned long z_offset(x_offset / 4);
            z_offset = (4 - z_offset) % 4;

            unsigned long quad_start(z_offset + begin);
            unsigned long quad_end(end - ((end - quad_start) % 4));

            //if (size < 32)
            {
                quad_end = begin;
                quad_start = begin;
            }

            for (unsigned long index(begin) ; index < quad_start ; ++index)
            {
                h[index] = f_0[index] + f_1[index] + f_2[index] + f_3[index] + f_4[index] +
                    f_5[index] + f_6[index] + f_7[index] + f_8[index];

                u[index] = (distribution_x[0] * f_0[index] +
                        distribution_x[1] * f_1[index] +
                        distribution_x[2] * f_2[index] +
                        distribution_x[3] * f_3[index] +
                        distribution_x[4] * f_4[index] +
                        distribution_x[5] * f_5[index] +
                        distribution_x[6] * f_6[index] +
                        distribution_x[7] * f_7[index] +
                        distribution_x[8] * f_8[index]) / h[index];

                v[index] = (distribution_y[0] * f_0[index] +
                        distribution_y[1] * f_1[index] +
                        distribution_y[2] * f_2[index] +
                        distribution_y[3] * f_3[index] +
                        distribution_y[4] * f_4[index] +
                        distribution_y[5] * f_5[index] +
                        distribution_y[6] * f_6[index] +
                        distribution_y[7] * f_7[index] +
                        distribution_y[8] * f_8[index]) / h[index];
            }
            for (unsigned long index(quad_end) ; index < end ; ++index)
            {
                h[index] = f_0[index] + f_1[index] + f_2[index] + f_3[index] + f_4[index] +
                    f_5[index] + f_6[index] + f_7[index] + f_8[index];

                u[index] = (distribution_x[0] * f_0[index] +
                        distribution_x[1] * f_1[index] +
                        distribution_x[2] * f_2[index] +
                        distribution_x[3] * f_3[index] +
                        distribution_x[4] * f_4[index] +
                        distribution_x[5] * f_5[index] +
                        distribution_x[6] * f_6[index] +
                        distribution_x[7] * f_7[index] +
                        distribution_x[8] * f_8[index]) / h[index];

                v[index] = (distribution_y[0] * f_0[index] +
                        distribution_y[1] * f_1[index] +
                        distribution_y[2] * f_2[index] +
                        distribution_y[3] * f_3[index] +
                        distribution_y[4] * f_4[index] +
                        distribution_y[5] * f_5[index] +
                        distribution_y[6] * f_6[index] +
                        distribution_y[7] * f_7[index] +
                        distribution_y[8] * f_8[index]) / h[index];
            }


            /// \todo use sse intrinsics
            for (unsigned long index(quad_start) ; index < quad_end ; index += 4)
            {
                /*__m128 t1, t2, m1, m2;
                  __m128 scal1, scal2;

                  m1 = _mm_mul_ps(scal1, gv);
                  m2 = _mm_load_ps(h + index);
                  _mm_store_ps(f_eq_0 + index, m1);*/
            }
        }

        void extraction_grid(unsigned long begin, unsigned long end,
                double * distribution_x, double * distribution_y,
                double * h, double * u, double * v,
                double * f_0, double * f_1, double * f_2,
                double * f_3, double * f_4, double * f_5,
                double * f_6, double * f_7, double * f_8)
        {
            //unsigned long size(end - begin);

            unsigned long x_address((unsigned long)&h[begin]);
            unsigned long x_offset(x_address % 16);

            unsigned long z_offset(x_offset / 8);

            unsigned long quad_start(z_offset + begin);
            unsigned long quad_end(end - ((end - quad_start) % 2));

            //if (size < 32)
            {
                quad_end = begin;
                quad_start = begin;
            }

            for (unsigned long index(begin) ; index < quad_start ; ++index)
            {
                h[index] = f_0[index] + f_1[index] + f_2[index] + f_3[index] + f_4[index] +
                    f_5[index] + f_6[index] + f_7[index] + f_8[index];

                u[index] = (distribution_x[0] * f_0[index] +
                        distribution_x[1] * f_1[index] +
                        distribution_x[2] * f_2[index] +
                        distribution_x[3] * f_3[index] +
                        distribution_x[4] * f_4[index] +
                        distribution_x[5] * f_5[index] +
                        distribution_x[6] * f_6[index] +
                        distribution_x[7] * f_7[index] +
                        distribution_x[8] * f_8[index]) / h[index];

                v[index] = (distribution_y[0] * f_0[index] +
                        distribution_y[1] * f_1[index] +
                        distribution_y[2] * f_2[index] +
                        distribution_y[3] * f_3[index] +
                        distribution_y[4] * f_4[index] +
                        distribution_y[5] * f_5[index] +
                        distribution_y[6] * f_6[index] +
                        distribution_y[7] * f_7[index] +
                        distribution_y[8] * f_8[index]) / h[index];
            }
            for (unsigned long index(quad_end) ; index < end ; ++index)
            {
                h[index] = f_0[index] + f_1[index] + f_2[index] + f_3[index] + f_4[index] +
                    f_5[index] + f_6[index] + f_7[index] + f_8[index];

                u[index] = (distribution_x[0] * f_0[index] +
                        distribution_x[1] * f_1[index] +
                        distribution_x[2] * f_2[index] +
                        distribution_x[3] * f_3[index] +
                        distribution_x[4] * f_4[index] +
                        distribution_x[5] * f_5[index] +
                        distribution_x[6] * f_6[index] +
                        distribution_x[7] * f_7[index] +
                        distribution_x[8] * f_8[index]) / h[index];

                v[index] = (distribution_y[0] * f_0[index] +
                        distribution_y[1] * f_1[index] +
                        distribution_y[2] * f_2[index] +
                        distribution_y[3] * f_3[index] +
                        distribution_y[4] * f_4[index] +
                        distribution_y[5] * f_5[index] +
                        distribution_y[6] * f_6[index] +
                        distribution_y[7] * f_7[index] +
                        distribution_y[8] * f_8[index]) / h[index];
            }


            /// \todo use sse intrinsics
            for (unsigned long index(quad_start) ; index < quad_end ; index += 4)
            {
                /*__m128 t1, t2, m1, m2;
                  __m128 scal1, scal2;

                  m1 = _mm_mul_ps(scal1, gv);
                  m2 = _mm_load_ps(h + index);
                  _mm_store_ps(f_eq_0 + index, m1);*/
            }
        }
    }
}
