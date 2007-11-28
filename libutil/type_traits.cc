/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * This file is part of the Utility C++ library. LibUtil is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibUtil is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <libutil/type_traits.hh>

#if defined(__ALTIVEC__)
#  include <altivec.h>
#  include <cstring>
#elif defined(__SSE__)
#  include <xmmintrin.h>
#  include <emmintrin.h>
#  include <cstring>
#else /* Unoptimised */
#  include <cstring>
#endif

namespace honei
{
    namespace intern
    {
#if defined(__ALTIVEC__)

        template <> void
        PODTraits<float>::copy(const float * source, float * dest, std::size_t count)
        {
            std::memcpy(dest, source, sizeof(float) * count);
        }

        template <> void
        PODTraits<double>::copy(const double * source, double * dest, std::size_t count)
        {
            std::memcpy(dest, source, sizeof(double) * count);
        }

        template <> void
        PODTraits<float>::fill(float * dest, std::size_t count, const float & proto)
        {
            vector float vectorised_proto = { proto, proto, proto, proto };
            vector float * vectorised_dest(reinterpret_cast<vector float *>(ptr));
            unsigned offset((reinterpret_cast<unsigned long long>(ptr) & 0xF) / sizeof(float));

            for (unsigned i(0) ; i < std::min(n, offset) ; ++i)
            {
                dest[i] = v;
            }

            for (unsigned i(sizeof(vector float) - sizeof(float) * offset) ; i < (n - offset) * sizeof(float) ;
                    i += sizeof(vector float))
            {
                vec_st(vectorised_proto, i, vectorised_dest);
            }

            for (unsigned i((n - offset) & ~0x3) ; i < n ; ++i)
            {
                dest[i] = proto;
            }
        }

        template <> void
        PODTraits<double>::fill(double * dest, std::size_t count, const double & proto)
        {
            vector double vectorised_proto = { proto, proto, proto, proto };
            vector double * vectorised_dest(reinterpret_cast<vector double *>(ptr));
            unsigned offset((reinterpret_cast<unsigned long long>(ptr) & 0xF) / sizeof(double));

            for (unsigned i(0) ; i < std::min(n, offset) ; ++i)
            {
                dest[i] = v;
            }

            for (unsigned i(sizeof(vector double) - sizeof(double) * offset) ; i < (n - offset) * sizeof(double) ;
                    i += sizeof(vector double))
            {
                vec_st(vectorised_proto, i, vectorised_dest);
            }

            for (unsigned i((n - offset) & ~0x1) ; i < n ; ++i)
            {
                dest[i] = proto;
            }
        }

#elif defined(__SSE__)

        template <> void
        PODTraits<float>::copy(const float * source, float * dest, std::size_t count)
        {
            __m128 m1;
            unsigned long dest_address = (unsigned long)dest;
            unsigned long dest_offset = dest_address % 16;

            unsigned long x_offset(dest_offset / 4);
            x_offset = (4 - x_offset) % 4;

            unsigned long quad_start = x_offset;
            unsigned long quad_end(count - ((count - quad_start) % 4));
            if (quad_end < 4 || quad_end > count)
            {
                quad_end = count;
                quad_start = count;
            }

            for (unsigned long index = quad_start ; index < quad_end ; index += 4)
            {
                m1 = _mm_loadu_ps(source + index);
                _mm_stream_ps(dest + index, m1);
            }

            for (unsigned long index (0) ; index < quad_start ; index++)
            {
                dest[index] = source[index];
            }
            for (unsigned long index = quad_end ; index < count ; index++)
            {
                dest[index] = source[index];
            }
        }

        template <> void
        PODTraits<double>::copy(const double * source, double * dest, std::size_t count)
        {
            __m128d m1;
            unsigned long dest_address = (unsigned long)dest;
            unsigned long dest_offset = dest_address % 16;

            unsigned long x_offset(dest_offset / 8);

            unsigned long quad_start = x_offset;
            unsigned long quad_end(count - ((count - quad_start) % 2));

            for (unsigned long index = quad_start ; index < quad_end ; index += 2)
            {
                m1 = _mm_loadu_pd(source + index);
                _mm_stream_pd(dest + index, m1);
            }

            for (unsigned long index (0) ; index < quad_start ; index++)
            {
                dest[index] = source[index];
            }
            for (unsigned long index = quad_end ; index < count ; index++)
            {
                dest[index] = source[index];
            }
        }

        template <> void
        PODTraits<float>::fill(float * dest, std::size_t count, const float & v)
        {
            __m128 m1;
            float __attribute__((aligned(16))) v_data;
            v_data= v;
            m1 = _mm_load_ps1(&v_data);

            unsigned long dest_address = (unsigned long)dest;
            unsigned long dest_offset = dest_address % 16;

            unsigned long x_offset(dest_offset / 4);
            x_offset = (4 - x_offset) % 4;

            unsigned long quad_start = x_offset;
            unsigned long quad_end(count - ((count - quad_start) % 4));
            if (quad_end < 4 || quad_end > count)
            {
                quad_end = count;
                quad_start = count;
            }

            for (unsigned long index = quad_start ; index < quad_end ; index += 4)
            {
                _mm_stream_ps(dest + index, m1);
            }

            for (unsigned long index (0) ; index < quad_start ; index++)
            {
                dest[index] = v;
            }
            for (unsigned long index = quad_end ; index < count ; index++)
            {
                dest[index] = v;
            }
        }

        template <> void
        PODTraits<double>::fill(double * dest, std::size_t count, const double & v)
        {
            __m128d m1;
            double __attribute__((aligned(16))) v_data;
            v_data= v;
            m1 = _mm_load_pd1(&v_data);

            unsigned long dest_address = (unsigned long)dest;
            unsigned long dest_offset = dest_address % 16;

            unsigned long x_offset(dest_offset / 8);

            unsigned long quad_start = x_offset;
            unsigned long quad_end(count - ((count - quad_start) % 2));

            for (unsigned long index = quad_start ; index < quad_end ; index += 2)
            {
                _mm_stream_pd(dest + index, m1);
            }

            for (unsigned long index (0) ; index < quad_start ; index++)
            {
                dest[index] = v;
            }
            for (unsigned long index = quad_end ; index < count ; index++)
            {
                dest[index] = v;
            }
        }

#else /* Unoptimised */

        template <> void
        PODTraits<float>::copy(const float * source, float * dest, std::size_t count)
        {
            std::memcpy(dest, source, sizeof(float) * count);
        }

        template <> void
        PODTraits<double>::copy(const double * source, double * dest, std::size_t count)
        {
            std::memcpy(dest, source, sizeof(double) * count);
        }

        template <> void
        PODTraits<float>::fill(float * dest, std::size_t count, const float & proto)
        {
            for (unsigned i(0) ; i < count ; ++i)
            {
                dest[i] = proto;
            }
        }

        template <> void
        PODTraits<double>::fill(double * dest, std::size_t count, const double & proto)
        {
            for (unsigned i(0) ; i < count ; ++i)
            {
                dest[i] = proto;
            }
        }

#endif

    }
}
