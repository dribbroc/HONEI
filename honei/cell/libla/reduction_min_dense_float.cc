/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#include <honei/cell/cell.hh>
#include <honei/cell/libla/operations.hh>

#include <spu_intrinsics.h>

namespace honei
{
    namespace cell
    {
        namespace implementation
        {
            vector float reduction_min_dense_float(const vector float & accumulator, vector float * elements,
                    const unsigned size)
            {
                vector float result(accumulator);

                for (unsigned i(0) ; i < size ; ++i)
                {
                    vector unsigned is_smaller_than(spu_cmpgt(elements[i], result));
                    result = spu_sel(elements[i], result, is_smaller_than);
                }

                return result;
            }
        }

        namespace operations
        {
            Operation<1, float, rtm_mail> reduction_min_dense_float = {
                &biggest_float,
                &implementation::reduction_min_dense_float,
                &min_float
            };
        }
    }
}
