/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
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

#include <honei/cell/cell.hh>
#include <honei/cell/libla/operations.hh>
#include <honei/cell/libutil/allocator.hh>
#include <honei/cell/libutil/transfer.hh>

#include <spu_intrinsics.h>

namespace honei
{
    namespace cell
    {
        namespace implementation
        {
            void product_banded_matrix_dense_vector_double(vector double * a, const vector double * b, const vector double * c,
                    const unsigned size, vector double & b_carry, const unsigned b_offset, vector double & c_carry,
                    const unsigned c_offset, const double optional_scalar)
            {
                vector double scalar_vector = spu_splats(optional_scalar);

                for (unsigned i(0) ; i < size ; ++i)
                {
                    extract(b_carry, b[i], b_offset);
                    extract(c_carry, c[i], c_offset);
                    a[i] = spu_madd(b_carry, c_carry, a[i]);
                    b_carry = b[i];
                    c_carry = c[i];
                }
            }
        }

        namespace operations
        {
            Operation<3, double, rtm_dma> product_banded_matrix_dense_vector_double = {
                &implementation::product_banded_matrix_dense_vector_double
            };
        }
    }
}
