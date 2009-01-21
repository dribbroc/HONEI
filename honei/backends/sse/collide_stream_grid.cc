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
        void collide_stream_grid_dir_0(unsigned long begin, unsigned long end, float tau,
                float * f_temp_0, float * f_0, float * f_eq_0)
        {
            unsigned long size(end - begin);

            unsigned long x_address((unsigned long)&f_temp_0[begin]);
            unsigned long x_offset(x_address % 16);

            unsigned long z_offset(x_offset / 4);
            z_offset = (4 - z_offset) % 4;

            unsigned long quad_start(z_offset + begin);
            unsigned long quad_end(end - ((end - quad_start) % 4));

            if (size < 16)
            {
                quad_end = begin;
                quad_start = begin;
            }

            for (unsigned long index(begin) ; index < quad_start ; ++index)
            {
                f_temp_0[index] = f_0[index] - (f_0[index] - f_eq_0[index])/tau;
            }
            for (unsigned long index(quad_end) ; index < end ; ++index)
            {
                f_temp_0[index] = f_0[index] - (f_0[index] - f_eq_0[index])/tau;
            }

            __m128 tauv = _mm_set1_ps(tau);

            for (unsigned long index(quad_start) ; index < quad_end ; index += 4)
            {
                __m128 m1, m2;

                m1 = _mm_load_ps(f_0 + index);
                m2 = _mm_load_ps(f_eq_0 + index);
                m2 = _mm_sub_ps(m1, m2);
                m2 = _mm_div_ps(m2, tauv);
                m2 = _mm_sub_ps(m1, m2);
                _mm_store_ps(f_temp_0 + index, m2);
            }
        }

        void collide_stream_grid_dir_n(unsigned long end, float tau,
                unsigned long * dir, unsigned long * dir_index,
                float * f_temp, float * f, float * f_eq)
        {
            __m128 tauv = _mm_set1_ps(tau);
            for (unsigned long begin(0), half(0) ; begin < end; begin+=2, ++half)
            {
                //for (unsigned long i(dir_index[begin]), offset(0) ; i < dir_index[begin + 1] ; ++i, ++offset)
                //{
                //}
                unsigned long start(dir_index[begin]);
                unsigned long size(dir_index[begin + 1] - start);

                unsigned long x_address((unsigned long)&f[start]);
                unsigned long x_offset(x_address % 16);

                unsigned long z_offset(x_offset / 4);
                z_offset = (4 - z_offset) % 4;

                unsigned long quad_start(z_offset + start);
                unsigned long quad_end(dir_index[begin + 1] - ((dir_index[begin + 1] - quad_start) % 4));

                if (size < 16)
                {
                    quad_end = start;
                    quad_start = start;
                }

                unsigned long offset(0);
                for (unsigned long index(start) ; index < quad_start ; ++index, ++offset)
                {
                    f_temp[dir[half] + offset] = f[index] - (f[index] - f_eq[index])/tau;
                }
                offset = quad_end - dir_index[begin];
                for (unsigned long index(quad_end) ; index < dir_index[begin + 1] ; ++index, ++offset)
                {
                    f_temp[dir[half] + offset] = f[index] - (f[index] - f_eq[index])/tau;
                }

                offset = quad_start - dir_index[begin];
                for (unsigned long index(quad_start) ; index < quad_end ; index += 4, offset += 4)
                {
                    //f_temp[dir[half] + offset] = f[index] - (f[index] - f_eq[index])/tau;
                    __m128 m1, m2;

                    m1 = _mm_load_ps(f + index);
                    m2 = _mm_load_ps(f_eq + index);
                    m2 = _mm_sub_ps(m1, m2);
                    m2 = _mm_div_ps(m2, tauv);
                    m2 = _mm_sub_ps(m1, m2);
                    _mm_storeu_ps(f_temp + dir[half] + offset, m2);
                }
            }
        }

        void collide_stream_grid_dir_0(unsigned long begin, unsigned long end, double tau,
                double * f_temp_0, double * f_0, double * f_eq_0)
        {
            unsigned long size(end - begin);

            unsigned long x_address((unsigned long)&f_temp_0[begin]);
            unsigned long x_offset(x_address % 16);

            unsigned long z_offset(x_offset / 8);

            unsigned long quad_start(z_offset + begin);
            unsigned long quad_end(end - ((end - quad_start) % 2));

            if (size < 16)
            {
                quad_end = begin;
                quad_start = begin;
            }

            for (unsigned long index(begin) ; index < quad_start ; ++index)
            {
                f_temp_0[index] = f_0[index] - (f_0[index] - f_eq_0[index])/tau;
            }
            for (unsigned long index(quad_end) ; index < end ; ++index)
            {
                f_temp_0[index] = f_0[index] - (f_0[index] - f_eq_0[index])/tau;
            }

            __m128d tauv = _mm_set1_pd(tau);

            for (unsigned long index(quad_start) ; index < quad_end ; index += 2)
            {
                __m128d m1, m2;

                m1 = _mm_load_pd(f_0 + index);
                m2 = _mm_load_pd(f_eq_0 + index);
                m2 = _mm_sub_pd(m1, m2);
                m2 = _mm_div_pd(m2, tauv);
                m2 = _mm_sub_pd(m1, m2);
                _mm_store_pd(f_temp_0 + index, m2);
            }
        }

        void collide_stream_grid_dir_n(unsigned long end, double tau,
                unsigned long * dir, unsigned long * dir_index,
                double * f_temp, double * f, double * f_eq)
        {
            __m128d tauv = _mm_set1_pd(tau);
            for (unsigned long begin(0), half(0) ; begin < end; begin+=2, ++half)
            {
                //for (unsigned long i(dir_index[begin]), offset(0) ; i < dir_index[begin + 1] ; ++i, ++offset)
                //{
                //}
                unsigned long start(dir_index[begin]);
                unsigned long size(dir_index[begin + 1] - start);

                unsigned long x_address((unsigned long)&f[start]);
                unsigned long x_offset(x_address % 16);

                unsigned long z_offset(x_offset / 8);

                unsigned long quad_start(z_offset + start);
                unsigned long quad_end(dir_index[begin + 1] - ((dir_index[begin + 1] - quad_start) % 2));

                if (size < 16)
                {
                    quad_end = start;
                    quad_start = start;
                }

                unsigned long offset(0);
                for (unsigned long index(start) ; index < quad_start ; ++index, ++offset)
                {
                    f_temp[dir[half] + offset] = f[index] - (f[index] - f_eq[index])/tau;
                }
                offset = quad_end - dir_index[begin];
                for (unsigned long index(quad_end) ; index < dir_index[begin + 1] ; ++index, ++offset)
                {
                    f_temp[dir[half] + offset] = f[index] - (f[index] - f_eq[index])/tau;
                }

                offset = quad_start - dir_index[begin];
                for (unsigned long index(quad_start) ; index < quad_end ; index += 2, offset += 2)
                {
                    //f_temp[dir[half] + offset] = f[index] - (f[index] - f_eq[index])/tau;
                    __m128d m1, m2;

                    m1 = _mm_load_pd(f + index);
                    m2 = _mm_load_pd(f_eq + index);
                    m2 = _mm_sub_pd(m1, m2);
                    m2 = _mm_div_pd(m2, tauv);
                    m2 = _mm_sub_pd(m1, m2);
                    _mm_storeu_pd(f_temp + dir[half] + offset, m2);
                }
            }
        }
    }
}

