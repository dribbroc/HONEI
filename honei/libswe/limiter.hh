/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
 * Copyright (c) 2007 Joachim Messer <joachim.messer@t-online.de>
 * Copyright (c) 2007 Volker Jung <volker.m.jung@t-online.de>
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * This file is part of the SWE C++ library. LibSWE is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibSWE is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef LIBSWE_GUARD_LIMITER_HH
#define LIBSWE_GUARD_LIMITER_HH 1

//#include <libutil/tags.hh>

#include <algorithm>
#ifdef HONEI_SSE
#include <xmmintrin.h>
#include <emmintrin.h>
#endif

namespace honei
{
    using std::max;
    using std::min;

    /**
     * \brief MinModLimiter limits a value to the range [0, 1].
     *
     * \ingroup grplimiter
     **/
    class MinModLimiter
    {
        public:
            /// Evaluates with respect to the MinMod-Limiter specification.
            /// \{

            template <typename DataType_> DataType_ & operator() (DataType_ & scal) const
            {
                return scal =  max(DataType_(0), min(DataType_(1), scal));
            }

            template <typename DataType_> DataType_ operator() (const DataType_ & scal) const
            {
                return max(DataType_(0), min(DataType_(1), scal));
            }

            /// \}
    };
    static const MinModLimiter min_mod_limiter = MinModLimiter();

    /**
     * \brief SuperBeeLimiter limits a value to the range [0, 2].
     *
     * \ingroup grplibswe
     **/
    class SuperBeeLimiter
    {
        public:
            /// Evaluates with respect to the SuperBee-Limiter specification.
            /// \{

            template <typename DataType_> DataType_ & operator() (DataType_ & scal) const
            {
                return scal = max(DataType_(0), max(min(2 * scal, DataType_(1)), min(scal, DataType_(2))));
            }

            template <typename DataType_> DataType_ operator() (const DataType_ & scal) const
            {
                return max(DataType_(0), max(min(2 * scal, DataType_(1)), min(scal, DataType_(2))));
            }

            /// \}
    };
    static const SuperBeeLimiter super_bee_limiter = SuperBeeLimiter();

    /**
     * \brief MonotonizedCentralLimiter limits a value to the range [0, 2].
     *
     * \ingroup grplibswe
     **/
    class MonotonizedCentralLimiter
    {
        public:
            /// Evaluates with respect to the Monotonized-Central-Limiter specification.
            /// \{

            template <typename DataType_> DataType_ & operator() (DataType_ & scal) const
            {
                return scal = max(DataType_(0), min(min(2 * scal, DataType_(2)), (1 + scal ) / 2));
            }

            template <typename DataType_> DataType_ operator() (const DataType_ & scal) const
            {
                return max(DataType_(0), min(min(2 * scal, DataType_(2)), (1 + scal ) / 2));
            }

            /// \}
    };
    static const MonotonizedCentralLimiter monotonized_central_limiter = MonotonizedCentralLimiter();

    /**
     * \brief VanLeerLimiter limits a value to the range [0, 2].
     *
     * \ingroup grplibswe
     **/
    class VanLeerLimiter
    {
        public:
            /// Evaluates due to the VanLeer-Limiter specification.
            /// \{
            template <typename DataType_> DataType_ & operator() (DataType_ & scal) const
            {
                return scal = (scal <= DataType_(0) ? DataType_(0) : (2 * scal / (1 + scal)));
            }

            template <typename DataType_> DataType_ operator() (const DataType_ & scal) const
            {
                return (scal <= DataType_(0) ? DataType_(0) : (2 * scal / (1 + scal)));
            }

            /// \}
    };
    static const VanLeerLimiter van_leer_limiter = VanLeerLimiter();

#ifdef HONEI_SSE

    namespace intern
    {
        namespace sse
        {
            __m128  _mm_lim_ps(__m128 vec)
            {
                float __attribute__((aligned(16))) zero(0.f);
                float __attribute__((aligned(16))) one(1.f);

                __m128 mm0 = _mm_set_ps1(zero);
                __m128 mm1 = _mm_set_ps1(one);
                __m128 result = _mm_max_ps( mm0, _mm_min_ps(mm1, vec));

                return result;
            }

            __m128d _mm_lim_pd(__m128d &  vec)
            {
                double __attribute__((aligned(16))) zero(0.);
                double __attribute__((aligned(16))) one(1.);

                __m128d mm0 = _mm_set_pd1(zero);
                __m128d mm1 = _mm_set_pd1(one);
                __m128d result = _mm_max_pd( mm0, _mm_min_pd(mm1, vec));

                return result;
            }
        }
    }
#endif
}

#endif
