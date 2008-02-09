/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/attributes.hh>
#include <honei/libutil/type_traits.hh>

#if defined(__ALTIVEC__)
#  include <altivec.h>
#  include <algorithm>
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
            /// \todo load data, mask and and it and store using vec_stl
            vector float vectorised_proto = { proto, proto, proto, proto };
            std::size_t offset((reinterpret_cast<unsigned long long>(dest) & 0xF) / sizeof(float));

            for (unsigned i(0) ; i < std::min(count, offset) ; ++i)
            {
                dest[i] = proto;
            }

            for (unsigned i(std::min(count, offset) * sizeof(float)) ; i < (count - offset) * sizeof(float) ;
                    i += sizeof(vector float))
            {
                vec_stl(vectorised_proto, i, dest);
            }

            for (unsigned i((count - offset) & ~0x3) ; i < count ; ++i)
            {
                dest[i] = proto;
            }
        }

        template <> void
        PODTraits<double>::fill(double * dest, std::size_t count, const double & proto)
        {
            union
            {
                double d[2];
                vector unsigned int vui;
            } vectorised_proto = { { proto, proto } };
            vector unsigned int * vectorised_dest(reinterpret_cast<vector unsigned int *>(dest));
            std::size_t offset((reinterpret_cast<unsigned long long>(dest) & 0xF) / sizeof(double));

            /// \todo load data, mask and and it and store using vec_stl
            for (unsigned i(0) ; i < std::min(count, offset) ; ++i)
            {
                dest[i] = proto;
            }

            for (unsigned i(std::min(count, offset) * sizeof(double)) ; i < (count - offset) * sizeof(double) ;
                    i += 2 * sizeof(double))
            {
                vec_stl(vectorised_proto.vui, i, vectorised_dest);
            }

            for (unsigned i((count - offset) & ~0x3) ; i < count ; ++i)
            {
                dest[i] = proto;
            }
        }

        void PODConversionTraits<float>::convert(double * copy, const float * orig, std::size_t count)
        {
            const float * s(orig);
            for (double * d(copy), * d_end(copy + count) ; d != d_end ; ++d, ++s)
            {
                *d = *s;
            }
        }

        void PODConversionTraits<double>::convert(float * copy, const double * orig, std::size_t count)
        {
            const double * s(orig);
            for (float * d(copy), * d_end(copy + count) ; d != d_end ; ++d, ++s)
            {
                *d = *s;
            }
        }

#elif defined(__SSE__)

        template <> void
        PODTraits<float>::copy(const float * source, float * dest, std::size_t count)
        {
            __m128 m1, m2, m3, m4, m5, m6, m7, m8;
            unsigned long dest_address = (unsigned long)dest;
            unsigned long dest_offset = dest_address % 16;

            unsigned long x_offset(dest_offset / 4);
            x_offset = (4 - x_offset) % 4;

            unsigned long quad_start = x_offset;
            unsigned long quad_end(count - ((count - quad_start) % 32));
            if (count < 40)
            {
                quad_end = 0;
                quad_start = 0;
            }

            for (unsigned long index = quad_start ; index < quad_end ; index += 32)
            {
                m1 = _mm_loadu_ps(source + index);
                m2 = _mm_loadu_ps(source + index + 4);
                m3 = _mm_loadu_ps(source + index + 8);
                m4 = _mm_loadu_ps(source + index + 12);
                m5 = _mm_loadu_ps(source + index + 16);
                m6 = _mm_loadu_ps(source + index + 20);
                m7 = _mm_loadu_ps(source + index + 24);
                m8 = _mm_loadu_ps(source + index + 28);
                _mm_store_ps(dest + index, m1);
                _mm_store_ps(dest + index + 4, m2);
                _mm_store_ps(dest + index + 8, m3);
                _mm_store_ps(dest + index + 12, m4);
                _mm_store_ps(dest + index + 16, m5);
                _mm_store_ps(dest + index + 20, m6);
                _mm_store_ps(dest + index + 24, m7);
                _mm_store_ps(dest + index + 28, m8);
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
            __m128d m1, m2, m3, m4, m5, m6, m7, m8;
            unsigned long dest_address = (unsigned long)dest;
            unsigned long dest_offset = dest_address % 16;

            unsigned long x_offset(dest_offset / 8);

            unsigned long quad_start = x_offset;
            unsigned long quad_end(count - ((count - quad_start) % 16));
            if (count < 24)
            {
                quad_end = 0;
                quad_start = 0;
            }

            for (unsigned long index = quad_start ; index < quad_end ; index += 16)
            {
                m1 = _mm_loadu_pd(source + index);
                m2 = _mm_loadu_pd(source + index + 2);
                m3 = _mm_loadu_pd(source + index + 4);
                m4 = _mm_loadu_pd(source + index + 6);
                m5 = _mm_loadu_pd(source + index + 8);
                m6 = _mm_loadu_pd(source + index + 10);
                m7 = _mm_loadu_pd(source + index + 12);
                m8 = _mm_loadu_pd(source + index + 14);
                _mm_store_pd(dest + index, m1);
                _mm_store_pd(dest + index + 2, m2);
                _mm_store_pd(dest + index + 4, m3);
                _mm_store_pd(dest + index + 6, m4);
                _mm_store_pd(dest + index + 8, m5);
                _mm_store_pd(dest + index + 10, m6);
                _mm_store_pd(dest + index + 12, m7);
                _mm_store_pd(dest + index + 14, m8);
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
            float HONEI_ATTRIBUTE(aligned(16)) v_data;
            v_data= v;
            m1 = _mm_load1_ps(&v_data);

            unsigned long dest_address = (unsigned long)dest;
            unsigned long dest_offset = dest_address % 16;

            unsigned long x_offset(dest_offset / 4);
            x_offset = (4 - x_offset) % 4;

            unsigned long quad_start = x_offset;
            unsigned long quad_end(count - ((count - quad_start) % 32));
            if (count < 40)
            {
                quad_end = 0;
                quad_start = 0;
            }

            for (unsigned long index = quad_start ; index < quad_end ; index += 32)
            {
                _mm_store_ps(dest + index, m1);
                _mm_store_ps(dest + index + 4, m1);
                _mm_store_ps(dest + index + 8, m1);
                _mm_store_ps(dest + index + 12, m1);
                _mm_store_ps(dest + index + 16, m1);
                _mm_store_ps(dest + index + 20, m1);
                _mm_store_ps(dest + index + 24, m1);
                _mm_store_ps(dest + index + 28, m1);
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
            double HONEI_ATTRIBUTE(aligned(16)) v_data;
            v_data= v;
            m1 = _mm_load1_pd(&v_data);

            unsigned long dest_address = (unsigned long)dest;
            unsigned long dest_offset = dest_address % 16;

            unsigned long x_offset(dest_offset / 8);

            unsigned long quad_start = x_offset;
            unsigned long quad_end(count - ((count - quad_start) % 16));
            if (count < 24)
            {
                quad_end = 0;
                quad_start = 0;
            }

            for (unsigned long index = quad_start ; index < quad_end ; index += 16)
            {
                _mm_store_pd(dest + index, m1);
                _mm_store_pd(dest + index + 2, m1);
                _mm_store_pd(dest + index + 4, m1);
                _mm_store_pd(dest + index + 6, m1);
                _mm_store_pd(dest + index + 8, m1);
                _mm_store_pd(dest + index + 10, m1);
                _mm_store_pd(dest + index + 12, m1);
                _mm_store_pd(dest + index + 14, m1);
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

        void PODConversionTraits<float>::convert(double * copy, const float * orig, std::size_t count)
        {
            __m128 m1, m2, m3, m4;
            __m128d m5, m6, m7, m8;

            unsigned long orig_address = (unsigned long)orig;
            unsigned long orig_offset = orig_address % 16;

            unsigned long x_offset(orig_offset / 4);
            x_offset = (4 - x_offset) % 4;

            unsigned long quad_start = x_offset;
            unsigned long quad_end(count - 4 - ((count - quad_start) % 8));
            if (count < 24)
            {
                quad_end = 0;
                quad_start = 0;
            }

            for (unsigned long index(quad_start) ; index < quad_end ; index += 8)
            {
                m1 = _mm_load_ps(orig + index);
                m2 = _mm_load_ps(orig + index + 2);
                m3 = _mm_load_ps(orig + index + 4);
                m4 = _mm_load_ps(orig + index + 6);

                m5 = _mm_cvtps_pd(m1);
                m6 = _mm_cvtps_pd(m2);
                m7 = _mm_cvtps_pd(m3);
                m8 = _mm_cvtps_pd(m4);

                _mm_store_pd(copy + index, m5);
                _mm_store_pd(copy + index + 2, m6);
                _mm_store_pd(copy + index + 4, m7);
                _mm_store_pd(copy + index + 6, m8);
            }

            for (unsigned long index(0) ; index < quad_start ; index++)
            {
                copy[index] = orig[index];
            }
            for (unsigned long index(quad_end) ; index < count ; index++)
            {
                copy[index] = orig[index];
            }
        }

        void PODConversionTraits<double>::convert(float * copy, const double * orig, std::size_t count)
        {
            //__m128d m1, m2, m3, m4;
            //__m128 m5, m6, m7, m8;

            unsigned long orig_address = (unsigned long)orig;
            unsigned long orig_offset = orig_address % 16;

            unsigned long x_offset(orig_offset / 4);
            x_offset = (4 - x_offset) % 4;

            unsigned long quad_start = x_offset;
            unsigned long quad_end(count - 8 - ((count - quad_start) % 8));
            if (count < 24)
            {
                quad_end = 0;
                quad_start = 0;
            }

            for (unsigned long index(quad_start) ; index < quad_end ; index += 8)
            {
                /*m1 = _mm_load_pd(orig + index);
                m2 = _mm_load_pd(orig + index + 2);
                m3 = _mm_load_pd(orig + index + 4);
                m4 = _mm_load_pd(orig + index + 6);

                m5 = _mm_cvtpd_ps(m1);
                m6 = _mm_cvtpd_ps(m2);
                m7 = _mm_cvtpd_ps(m3);
                m8 = _mm_cvtpd_ps(m4);

                _mm_store_ps(copy + index, m5);
                _mm_store_ps(copy + index + 2, m6);
                _mm_store_ps(copy + index + 4, m7);
                _mm_store_ps(copy + index + 6, m8);*/

                copy[index] = orig[index];
                copy[index + 1] = orig[index + 1];
                copy[index + 2] = orig[index + 2];
                copy[index + 3] = orig[index + 3];
                copy[index + 4] = orig[index + 4];
                copy[index + 5] = orig[index + 5];
                copy[index + 6] = orig[index + 6];
                copy[index + 7] = orig[index + 7];
            }

            for (unsigned long index (0) ; index < quad_start ; index++)
            {
                copy[index] = orig[index];
            }
            for (unsigned long index(quad_end) ; index < count ; index++)
            {
                copy[index] = orig[index];
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

        void PODConversionTraits<float>::convert(double * copy, const float * orig, std::size_t count)
        {
            const float * s(orig);
            for (double * d(copy), * d_end(copy + count) ; d != d_end ; ++d, ++s)
            {
                *d = *s;
            }
        }

        void PODConversionTraits<double>::convert(float * copy, const double * orig, std::size_t count)
        {
            const double * s(orig);
            for (float * d(copy), * d_end(copy + count) ; d != d_end ; ++d, ++s)
            {
                *d = *s;
            }
        }
#endif

    }
}
