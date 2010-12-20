/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the LBM C++ library. LBM is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LBM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#pragma once
#ifndef LBM_GUARD_LBM_LIMITER_HH
#define LBM_GUARD_LBM_LIMITER_HH 1

#include <algorithm>
#include <limits>

namespace honei
{
    namespace lbm
    {
        /**
         * \brief MinModLimiter limits a value to the range [0, 1].
         *
         * \ingroup grplbmoperations
         **/

        template<typename Tag_>
            class MinModLimiter
            {
                public:
                    /// Evaluates with respect to the MinMod-Limiter specification.
                    /// \{

                    template <typename DataType_>
                        static DataType_ value(const DataType_ & scal)
                        {
                            return std::max(DataType_(0), std::min(DataType_(0.1), scal));
                        }

                    /// \}
            };
    }
}
#endif
