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

#if defined(__ALTIVEC__) || defined(HONEI_CELL)
#  include <altivec.h>
#  include <cstring>
#else /* Unoptimised */
#  include <cstring>
#endif

namespace honei
{
    namespace intern
    {
#if defined(__ALTIVEC__) || defined(HONEI_CELL)

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
        PODTraits<double>::fill(float * dest, std::size_t count, const double & proto)
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
