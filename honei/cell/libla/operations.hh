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

#ifndef CELL_GUARD_LIBLA_OPERATIONS_HH
#define CELL_GUARD_LIBLA_OPERATIONS_HH 1

#include <honei/cell/libutil/operations.hh>

namespace honei
{
    namespace cell
    {
        namespace operations
        {
            extern Operation<1, float, rtm_dma> element_inverse_float;

            extern Operation<1, float, rtm_mail> norm_max_dense_float;

            extern Operation<1, float, rtm_mail> norm_l_one_dense_float;

            extern Operation<1, float, rtm_mail> norm_l_two_dense_float;

            extern Operation<1, float, rtm_mail> reduction_max_dense_float;

            extern Operation<1, float, rtm_mail> reduction_min_dense_float;

            extern Operation<1, double, rtm_mail> reduction_sum_dense_double;

            extern Operation<1, float, rtm_mail> reduction_sum_dense_float;

            extern Operation<1, float, rtm_dma> sum_dense_scalar_float;

            extern Operation<1, float, rtm_dma> fill_float;

            extern Operation<1, double, rtm_dma> fill_double;

            extern Operation<2, float, rtm_dma> sum_dense_dense_float;

            extern Operation<2, float, rtm_dma> difference_dense_dense_float;

            extern Operation<2, float, rtm_dma> element_product_dense_dense_float;

        }
    }
}

#endif
