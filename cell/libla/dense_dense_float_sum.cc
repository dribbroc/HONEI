/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Sven Mallach <sven.mallach@uni-dortmund.de>
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
#include <cell/libutil/allocator.hh>
#include <cell/libutil/debug.hh>

#include <spu_intrinsics.h>
#include <spu_mfcio.h>
#include <stdio.h>

using namespace honei;

/*
 * dense_dense_float_sum
 *
 * Calculate the sum of two dense entities.
 *
 * \size Default transfer buffer size in bytes.
 * \operand a Base address of first entity.
 * \operand b Base address of second entity.
 * \operand c Number of transfers needed.
 * \operand d Last transfer buffer size in bytes.
 */
int dense_dense_float_sum(const Instruction & inst)
{
    debug_dump();
    EffectiveAddress ea_a(inst.a.ea), ea_b(inst.b.ea), ea_result(inst.a.ea);
    unsigned counter(inst.c.u);

    allocator::Allocation * block_a[2] = { allocator::acquire_block(), allocator::acquire_block() };
    allocator::Allocation * block_b[2] = { allocator::acquire_block(), allocator::acquire_block() };

    Pointer<float> a[2] = { block_a[0]->address, block_a[1]->address };
    Pointer<float> b[2] = { block_b[0]->address, block_b[1]->address };

    unsigned current(1), next(2);
    debug_get(ea_a, a[current - 1].untyped, inst.size);
    mfc_get(a[current - 1].untyped, ea_a, inst.size, current, 0, 0);
    debug_get(ea_b, b[current - 1].untyped, inst.size);
    mfc_get(b[current - 1].untyped, ea_b, inst.size, current, 0, 0);

    while (counter > 1)
    {
        ea_a += inst.size;
        ea_b += inst.size;

        debug_get(ea_a, a[next - 1].untyped, inst.size);
        mfc_get(a[next - 1].untyped, ea_a, inst.size, next, 0, 0);
        debug_get(ea_b, b[next - 1].untyped, inst.size);
        mfc_get(b[next - 1].untyped, ea_b, inst.size, next, 0, 0);

        mfc_write_tag_mask(1 << current);
        mfc_read_tag_status_all();
        debug_dump();

        for (unsigned i(0) ; i < inst.size / sizeof(vector float) ; ++i)
        {
            a[current - 1].vectorised[i] = spu_add(a[current - 1].vectorised[i], b[current - 1].vectorised[i]);
        }

        debug_put(ea_result, a[current - 1].untyped, inst.size);
        mfc_putb(a[current - 1].untyped, ea_result, inst.size, next, 0, 0);

        --counter;
        ea_result += inst.size;

        unsigned temp(next);
        next = current;
        current = temp;
    }

    ea_a += inst.size;
    ea_b += inst.size;

    if (inst.d.u > 0)
    {
        debug_get(ea_a, a[next - 1].untyped, inst.d.u);
        mfc_get(a[next - 1].untyped, ea_a, inst.d.u, next, 0, 0);
        debug_get(ea_b, b[next - 1].untyped, inst.d.u);
        mfc_get(b[next - 1].untyped, ea_b, inst.d.u, next, 0, 0);
    }

    mfc_write_tag_mask(1 << current);
    mfc_read_tag_status_all();

    for (unsigned i(0) ; i < inst.size / sizeof(vector float) ; ++i)
    {
        a[current - 1].vectorised[i] = spu_add(a[current - 1].vectorised[i], b[current - 1].vectorised[i]);
    }

    debug_put(ea_result, a[current - 1].untyped, inst.size);
    mfc_putb(a[current - 1].untyped, ea_result, inst.size, next, 0, 0);

    --counter;
    ea_result += inst.size;

    unsigned temp(next);
    next = current;
    current = temp;

    if (inst.d.u > 0)
    {
        mfc_write_tag_mask(1 << current);
        mfc_read_tag_status_all();

        for (unsigned i(0) ; i < inst.d.u / sizeof(vector float) ; ++i)
        {
            a[current - 1].vectorised[i] = spu_add(a[current - 1].vectorised[i], b[current - 1].vectorised[i]);
        }

        debug_put(ea_result, a[current - 1].untyped, inst.d.u);
        mfc_putb(a[current - 1].untyped, ea_result, inst.d.u, next, 0, 0);

        mfc_write_tag_mask(1 << next);
        mfc_read_tag_status_all();
    }

    allocator::release_block(*block_a[0]);
    allocator::release_block(*block_a[1]);
    allocator::release_block(*block_b[0]);
    allocator::release_block(*block_b[1]);
}
