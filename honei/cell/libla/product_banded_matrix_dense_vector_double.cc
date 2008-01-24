/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
#include <honei/cell/libutil/allocator.hh>
#include <honei/cell/libutil/transfer.hh>

#include <spu_intrinsics.h>
#include <spu_mfcio.h>

using namespace honei::cell;

/*
 * product_banded_matrix_dense_vector_double
 *
 * \size Default transfer buffer size in bytes.
 * \operand a Base address of first entity.
 * \operand b Base address of second entity.
 * \operand c Base address of result entity.
 * \operand d Number of transfers needed.
 * \operand e Last transfer buffer size in bytes.
 * \operand f Band offset
 * \operand g SIMD x_offset
 */

void product_banded_matrix_dense_vector_double(const Instruction & inst)
{
    EffectiveAddress ea_a(inst.a.ea), ea_b(inst.b.ea), ea_r(inst.c.ea), ea_result(inst.c.ea);

    Allocation * block_a[2] = { acquire_block(), acquire_block() };
    Allocation * block_b[2] = { acquire_block(), acquire_block() };
    Allocation * block_r[2] = { acquire_block(), acquire_block() };

    Pointer<double> a[2] = { { block_a[0]->address }, { block_a[1]->address } };
    Pointer<double> b[2] = { { block_b[0]->address }, { block_b[1]->address } };
    Pointer<double> r[2] = { { block_r[0]->address }, { block_r[1]->address } };

    unsigned x_offset(inst.g.u);
    unsigned y_offset((2 - x_offset) % 2);

    unsigned counter(inst.d.u);
    unsigned size(counter > 1 ? inst.size : inst.e.u);
    unsigned nextsize;
    unsigned current(1), next(2);

    debug_get(ea_a, a[current -1].untyped, size * sizeof(double));
    mfc_get(a[current - 1].untyped, ea_a, size * sizeof(double), current, 0, 0);
    debug_get(ea_b, b[current -1].untyped, (size + x_offset + y_offset) * sizeof(double));
    mfc_get(b[current - 1].untyped, ea_b, (size + x_offset + y_offset) * sizeof(double), current, 0, 0);
    debug_get(ea_r, r[current -1].untyped, size * sizeof(double));
    mfc_get(r[current - 1].untyped, ea_r, size * sizeof(double), current, 0, 0);
    ea_a += size * sizeof(double);
    ea_b += size * sizeof(double);
    ea_r += size * sizeof(double);

    while (counter > 1)
    {
        nextsize = (counter == 2 ? inst.e.u : inst.size);

        debug_get(ea_a, a[next -1].untyped, nextsize * sizeof(double));
        mfc_get(a[next - 1].untyped, ea_a, nextsize * sizeof(double), next, 0, 0);
        debug_get(ea_b, b[next -1].untyped, (nextsize + x_offset + y_offset) * sizeof(double));
        mfc_get(b[next - 1].untyped, ea_b, (nextsize + x_offset + y_offset) * sizeof(double), next, 0, 0);
        debug_get(ea_r, r[next -1].untyped, nextsize * sizeof(double));
        mfc_get(r[next - 1].untyped, ea_r, nextsize * sizeof(double), next, 0, 0);
        ea_a += nextsize * sizeof(double);
        ea_b += nextsize * sizeof(double);
        ea_r += nextsize * sizeof(double);

        mfc_write_tag_mask(1 << current);
        mfc_read_tag_status_all();

        for (unsigned i(0) ; i < size / 2 ; i++)
        {
            extract(b[current - 1].vectorised[i], b[current - 1].vectorised[i + 1], x_offset);
            r[current - 1].vectorised[i] = spu_madd(a[current - 1].vectorised[i], b[current - 1].vectorised[i], r[current - 1].vectorised[i]);
        }
        debug_put(ea_result, r[current -1].untyped, size * sizeof(double));
        mfc_putb(r[current - 1].untyped, ea_result, size * sizeof(double), current, 0, 0);
        ea_result += size * sizeof(double);

        --counter;

        unsigned temp(next);
        next = current;
        current = temp;

        size = nextsize;
    }
    mfc_write_tag_mask(1 << current);
    mfc_read_tag_status_all();

    for (unsigned i(0) ; i < size / 2 ; i++)
    {
        extract(b[current - 1].vectorised[i], b[current - 1].vectorised[i + 1], x_offset);
        r[current - 1].vectorised[i] = spu_madd(a[current - 1].vectorised[i], b[current - 1].vectorised[i], r[current - 1].vectorised[i]);
    }

    debug_put(ea_result, r[current -1].untyped, size * sizeof(double));
    mfc_put(r[current - 1].untyped, ea_result, size * sizeof(double), current, 0, 0);

    mfc_write_tag_mask(1 << current);
    mfc_read_tag_status_all();

    release_block(*block_a[0]);
    release_block(*block_a[1]);
    release_block(*block_b[0]);
    release_block(*block_b[1]);
    release_block(*block_r[0]);
    release_block(*block_r[1]);
}
