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
#include <honei/backends/cell/spe/libutil/allocator.hh>
#include <honei/backends/cell/spe/libutil/transfer.hh>
#include <honei/backends/cell/spe/libutil/operation_framework.hh>

#include <spu_intrinsics.h>
#include <spu_mfcio.h>

/*
 * operation_dense_dense_float
 *
 * Calculate an operation with two dense entities.
 *
 * \size Default transfer buffer size in bytes.
 * \operand a Base address of first entity.
 * \operand b Base address of second entity.
 * \operand c Base address of third entity.
 * \operand d Base address of fourth entity.
 * \operand e Number of transfers needed.
 * \operand f Last transfer buffer size in bytes.
 * \operand g scalar
 * \operand h carry lower half
 * \operand i carry upper half
*/

namespace honei
{
    namespace cell
    {
        template <>
        void operation<Operation<4, float, rtm_dma> >(const Operation<4, float, rtm_dma> & operation,
                const Instruction & instruction)
        {
            EffectiveAddress ea_a(instruction.a.ea), ea_b(instruction.b.ea), ea_c(instruction.c.ea),
                ea_d(instruction.d.ea), ea_result(instruction.a.ea);

            Allocation * block_a[2] = { acquire_block(), acquire_block() };
            Allocation * block_b[2] = { acquire_block(), acquire_block() };
            Allocation * block_c[2] = { acquire_block(), acquire_block() };
            Allocation * block_d[2] = { acquire_block(), acquire_block() };

            Pointer<float> a[2] = { { block_a[0]->address }, { block_a[1]->address } };
            Pointer<float> b[2] = { { block_b[0]->address }, { block_b[1]->address } };
            Pointer<float> c[2] = { { block_c[0]->address }, { block_c[1]->address } };
            Pointer<float> d[2] = { { block_d[0]->address }, { block_d[1]->address } };

            unsigned counter(instruction.e.u);
            unsigned size(counter > 1 ? instruction.size : instruction.f.u);
            printf("size=%u, counter=%u\n", instruction.size, counter);
            unsigned nextsize;
            unsigned current(0), next(1);

            float scalar = instruction.g.f; // optional scalar value to be computed.
            vector float carry = { instruction.h.fa[0], instruction.h.fa[1], instruction.i.fa[0], instruction.i.fa[1] };

            debug_get(ea_a, a[current].untyped, size);
            mfc_get(a[current].untyped, ea_a, size, current, 0, 0);
            debug_get(ea_b, b[current].untyped, size);
            mfc_get(b[current].untyped, ea_b, size, current, 0, 0);
            debug_get(ea_c, c[current].untyped, size);
            mfc_get(c[current].untyped, ea_c, size, current, 0, 0);
            debug_get(ea_d, d[current].untyped, size);
            mfc_get(d[current].untyped, ea_d, size, current, 0, 0);

            ea_a += size;
            ea_b += size;
            ea_c += size;
            ea_d += size;

            operation.init(scalar, carry);

            while (counter > 1)
            {
                printf("loop started\n");
                nextsize = (counter == 2 ? instruction.f.u : instruction.size);

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

                printf("call calculate\n");

                operation.calculate(a[current].vectorised, b[current].vectorised, c[current].vectorised,
                        d[current].vectorised, size / sizeof(vector float), scalar);

                printf("calculate returned\n");

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
                    d[current].vectorised, size / sizeof(vector float), scalar);

            debug_put(ea_result, a[current].untyped, size);
            mfc_put(a[current].untyped, ea_result, size, current, 0, 0);
            mfc_write_tag_mask(1 << current);
            mfc_read_tag_status_all();

            release_all_blocks();
            printf("offdma: done(%u)\n", instruction.opcode);
        }
    }
}

