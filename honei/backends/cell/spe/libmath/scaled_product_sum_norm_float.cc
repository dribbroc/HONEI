/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2009 Danny van Dyk <dyk@honei.org>
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
#include <honei/backends/cell/spe/libutil/allocator.hh>
#include <honei/backends/cell/spe/libutil/transfer.hh>

#include <spu_intrinsics.h>

namespace honei
{
    namespace cell
    {
        namespace implementation
        {
            vector float scaled_product_sum_norm_float(const vector float & accumulator, vector float & b_carry, const vector float * a,
                    const vector float * b, unsigned size, unsigned b_offset, float alpha, float beta)
            {
                vector float acc = accumulator;
                const vector float alphav = spu_splats(alpha);
                const vector float betav = spu_splats(beta);

                for (unsigned i(0) ; i < size ; ++i)
                {
                    vector float b_with_offset = b_carry;
                    extract(b_with_offset, b[i], b_offset);

                    vector float sum = spu_madd(alphav, a[i], spu_mul(betav, b_with_offset));
                    acc = spu_madd(sum, sum, acc);

                    b_carry = b[i];
                }

                return acc;
            }
        }

        namespace operations
        {
            Operation<2, float, rtm_mail> scaled_product_sum_norm_float = {
                &zero_float,
                &implementation::scaled_product_sum_norm_float,
                &sum_float
            };
        }
    }
}
