/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Sven Mallach <sven.mallach@honei.org>
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
 * \operand c Base address of third entity.
 * \operand d Number of transfers needed.
 * \operand e Last transfer buffer size in bytes.
 * \operand f alignment offset of entity b.
 * \operand g alignment offset of entity c.
 * \operand h first element of carry_b (not transfered via dma)
 * \operand i second element of carry_b (not transfered via dma)
 * \operand j first element of carry_c (not transfered via dma)
 * \operand k second element of carry_c (not transfered via dma)
 * \operand l scalar
 * \operand m scalar
 * \operand n scalar
*/

namespace honei
{
    namespace cell
    {
        template <>
        void operation<Operation<3, double, rtm_dma> >(const Operation<3, double, rtm_dma> & operation,
                const Instruction & instruction)
        {
            EffectiveAddress ea_a(instruction.a.ea), ea_b(instruction.b.ea), ea_c(instruction.c.ea),
                ea_result(instruction.a.ea);

            Allocation * block_a[2] = { acquire_block(), acquire_block() };
            Allocation * block_b[2] = { acquire_block(), acquire_block() };
            Allocation * block_c[2] = { acquire_block(), acquire_block() };

            Pointer<double> a[2] = { { block_a[0]->address }, { block_a[1]->address } };
            Pointer<double> b[2] = { { block_b[0]->address }, { block_b[1]->address } };
            Pointer<double> c[2] = { { block_c[0]->address }, { block_c[1]->address } };

            unsigned counter(instruction.d.u);
            unsigned size(counter > 1 ? instruction.size : instruction.e.u);
            unsigned nextsize;
            unsigned current(0), next(1);

            double a_scalar = instruction.l.d; // optional scalar value to be computed.
            double b_scalar = instruction.m.d; // optional scalar value to be computed.
            double c_scalar = instruction.n.d; // optional scalar value to be computed.

            unsigned b_offset(instruction.f.u);
            unsigned c_offset(instruction.g.u);
            vector double b_carry = { instruction.h.d, instruction.i.d};
            vector double c_carry = { instruction.j.d, instruction.k.d};

            debug_get(ea_a, a[current].untyped, size);
            mfc_get(a[current].untyped, ea_a, size, current, 0, 0);
            debug_get(ea_b, b[current].untyped, size);
            mfc_get(b[current].untyped, ea_b, size, current, 0, 0);
            debug_get(ea_c, c[current].untyped, size);
            mfc_get(c[current].untyped, ea_c, size, current, 0, 0);
            ea_a += size;
            ea_b += size;
            ea_c += size;

            while (counter > 1)
            {
                nextsize = (counter == 2 ? instruction.e.u : instruction.size);

                debug_get(ea_a, a[next].untyped, nextsize);
                mfc_get(a[next].untyped, ea_a, nextsize, next, 0, 0);
                debug_get(ea_b, b[next].untyped, nextsize);
                mfc_get(b[next].untyped, ea_b, nextsize, next, 0, 0);
                debug_get(ea_c, c[next].untyped, nextsize);
                mfc_get(c[next].untyped, ea_c, nextsize, next, 0, 0);
                ea_a += nextsize;
                ea_b += nextsize;
                ea_c += nextsize;

                mfc_write_tag_mask(1 << current);
                mfc_read_tag_status_all();

                operation.calculate(a[current].vectorised, b[current].vectorised, c[current].vectorised,
                        size / sizeof(vector double), b_carry, b_offset, c_carry, c_offset, a_scalar, b_scalar, c_scalar);

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

            operation.calculate(a[current].vectorised, b[current].vectorised, c[current].vectorised,
                    size / sizeof(vector double), b_carry, b_offset, c_carry, c_offset, a_scalar, b_scalar, c_scalar);

            debug_put(ea_result, a[current].untyped, size);
            mfc_put(a[current].untyped, ea_result, size, current, 0, 0);
            mfc_write_tag_mask(1 << current);
            mfc_read_tag_status_all();

            release_all_blocks();
        }
    }
}

