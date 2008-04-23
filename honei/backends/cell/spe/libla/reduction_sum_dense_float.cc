/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (C) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * This file is part of the LA C++ library. LibLa is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibLa is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/backends/cell/cell.hh>
#include <honei/backends/cell/spe/libla/operations.hh>

#include <spu_intrinsics.h>

namespace honei
{
    namespace cell
    {
        namespace implementation
        {
            vector float reduction_sum_dense_float(const vector float & accumulator, vector float * elements,
                    const unsigned size, const float optional_scalar)
            {
                vector float result(accumulator);

                for (unsigned i(0) ; i < size ; ++i)
                {
                    result = spu_add(elements[i], result);
                }

                return result;
            }
        }

        namespace operations
        {
            Operation<1, float, rtm_mail> reduction_sum_dense_float = {
                &zero_float,
                &implementation::reduction_sum_dense_float,
                &sum_float
            };
        }
    }
}

