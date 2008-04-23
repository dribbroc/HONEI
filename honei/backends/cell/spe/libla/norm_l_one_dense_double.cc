/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (C) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
            vector double norm_l_one_dense_double(const vector double & accumulator, vector double * elements, const unsigned size, const double optional_scalar)
            {
                const vector unsigned int mask = { 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF };

                vector double result(accumulator);

                for (unsigned i(0) ; i < size ; ++i)
                {
                    vector double absolute_values(
                            reinterpret_cast<vector double>(spu_and(reinterpret_cast<vector unsigned int>(elements[i]),
                                    mask)));

                    result = spu_add(absolute_values, result);
                }

                return result;
            }
        }

        namespace operations
        {
            Operation<1, double, rtm_mail> norm_l_one_dense_double = {
                &zero_double,
                &implementation::norm_l_one_dense_double,
                &sum_double
            };
        }
    }
}

