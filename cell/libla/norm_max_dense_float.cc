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

#include <cell/cell.hh>
#include <cell/libla/operations.hh>

#include <spu_intrinsics.h>

namespace honei
{
    namespace cell
    {
        namespace implementation
        {
            vector float norm_max_dense_float(const vector float & accumulator, vector float * elements, const unsigned size)
            {
                vector float result(accumulator);

                for (unsigned i(0) ; i < size ; ++i)
                {
                    vector unsigned int is_greater_than(spu_cmpabsgt(result, elements[i]));
                    result = spu_sel(elements[i], result, is_greater_than);
                }

                return result;
            }
        }

        namespace operations
        {
            Operation<1, float, rtm_mail> norm_max_dense_float = {
                &zero_float,
                &implementation::norm_max_dense_float,
                &max_float
            };
        }
    }
}

