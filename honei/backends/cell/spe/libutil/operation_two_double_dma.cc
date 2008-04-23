/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Sven Mallach <sven.mallach@honei.org>
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@web.de>
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
#include <honei/backends/cell/spe/libutil/allocator.hh>
#include <honei/backends/cell/spe/libutil/transfer.hh>
#include <honei/backends/cell/spe/libutil/operation_framework.hh>

#include <spu_intrinsics.h>
#include <spu_mfcio.h>

/*
 * operation_dense_dense_double
 *
 * Calculate an operation with two dense entities.
 *
 * \size Default transfer buffer size in bytes.
 * \operand a Base address of first entity.
 * \operand b Base address of second entity.
 * \operand c Number of transfers needed.
 * \operand d Last transfer buffer size in bytes.
 * \operand e Offset of entity b.
 * \operand f First element of carry_b (not transfered via dma).
 * \operand g Second element of carry_b (not transfered via dma).
 * \operand h Scalar value.
*/

namespace honei
{
    namespace cell
    {
        template <>
        void operation<Operation<2, double, rtm_dma> >(const Operation<2, double, rtm_dma> & operation,
                const Instruction & instruction)
        {
            EffectiveAddress ea_a(instruction.a.ea), ea_b(instruction.b.ea), ea_result(instruction.a.ea);

            Allocation * block_a[2] = { acquire_block(), acquire_block() };
            Allocation * block_b[2] = { acquire_block(), acquire_block() };

            Pointer<double> a[2] = { { block_a[0]->address }, { block_a[1]->address } };
            Pointer<double> b[2] = { { block_b[0]->address }, { block_b[1]->address } };

            unsigned counter(instruction.c.u);
            unsigned size(counter > 1 ? instruction.size : instruction.d.u);
            unsigned nextsize;
            unsigned current(0), next(1);

            double scalar = instruction.h.d;

            unsigned b_offset(instruction.e.u);
            vector double b_carry = { instruction.f.d, instruction.g.d };

            debug_get(ea_a, a[current].untyped, size);
            mfc_get(a[current].untyped, ea_a, size, current, 0, 0);
            debug_get(ea_b, b[current].untyped, size);
            mfc_get(b[current].untyped, ea_b, size, current, 0, 0);
            ea_a += size;
            ea_b += size;

            while (counter > 1)
            {
                nextsize = (counter == 2 ? instruction.d.u : instruction.size);

                debug_get(ea_a, a[next].untyped, nextsize);
                mfc_get(a[next].untyped, ea_a, nextsize, next, 0, 0);
                debug_get(ea_b, b[next].untyped, nextsize);
                mfc_get(b[next].untyped, ea_b, nextsize, next, 0, 0);
                ea_a += nextsize;
                ea_b += nextsize;

                mfc_write_tag_mask(1 << current);
                mfc_read_tag_status_all();

                operation.calculate(a[current].vectorised, b[current].vectorised, size / sizeof(vector double), b_carry, b_offset, scalar);

                debug_put(ea_result, a[current].untyped, size);
                mfc_putb(a[current].untyped, ea_result, size, current, 0, 0);
                ea_result += size;

                --counter;

                unsigned temp(next);
                next = current;
                current = temp;

                size = nextsize;
            }

            mfc_write_tag_mask(1 << current);
            mfc_read_tag_status_all();

            operation.calculate(a[current].vectorised, b[current].vectorised, size / sizeof(vector double), b_carry, b_offset, scalar);

            debug_put(ea_result, a[current].untyped, size);
            mfc_put(a[current].untyped, ea_result, size, current, 0, 0);
            mfc_write_tag_mask(1 << current);
            mfc_read_tag_status_all();

            release_all_blocks();
        }
    }
}

