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
 * dense_dense_float_dot_product
 *
 * Calculate the dot product of two dense vectors.
 *
 * \size Default transfer buffer size in bytes.
 * \operand a Address of return value.
 * \operand b Base address of first vector.
 * \operand c Base address of second vector.
 * \operand d Number of transfers needed.
 * \operand e Last transfer buffer size in bytes.
 */
unsigned dense_dense_float_dot_product(const Instruction & inst)
{
    EffectiveAddress ea_a(inst.b.ea), ea_b(inst.c.ea);
    unsigned counter(inst.d.u);

    allocator::Allocation * block_a[2] = { allocator::acquire_block(), allocator::acquire_block() };
    allocator::Allocation * block_b[2] = { allocator::acquire_block(), allocator::acquire_block() };

    Pointer<float> a[2] = { block_a[0]->address, block_a[1]->address };
    Pointer<float> b[2] = { block_b[0]->address, block_b[1]->address };

    MailableResult<float> result;
    Subscriptable<float> intermediate = { spu_splats(0.0f) };

    unsigned current(1), next(2);
    //debug_get(ea_a, a[current - 1].untyped, inst.size);
    mfc_get(a[current - 1].untyped, ea_a, inst.size, current, 0, 0);
    //debug_get(ea_b, b[current - 1].untyped, inst.size);
    mfc_get(b[current - 1].untyped, ea_b, inst.size, current, 0, 0);

    while (counter > 1)
    {
        ea_a += inst.size;
        ea_b += inst.size;

        //debug_get(ea_a, a[next - 1].untyped, inst.size);
        mfc_get(a[next - 1].untyped, ea_a, inst.size, next, 0, 0);
        //debug_get(ea_b, b[next - 1].untyped, inst.size);
        mfc_get(b[next - 1].untyped, ea_b, inst.size, next, 0, 0);

        mfc_write_tag_mask(1 << current);
        mfc_read_tag_status_all();

        for (unsigned i(0) ; i < inst.size / sizeof(vector float) ; ++i)
        {
            intermediate.value = spu_madd(a[current - 1].vectorised[i], b[current - 1].vectorised[i], intermediate.value);
        }


        --counter;

        unsigned temp(next);
        next = current;
        current = temp;
    }

    ea_a += inst.size;
    ea_b += inst.size;

    if (inst.e.u > 0)
    {
        //debug_get(ea_a, a[next - 1].untyped, inst.d.u);
        mfc_get(a[next - 1].untyped, ea_a, inst.d.u, next, 0, 0);
        //debug_get(ea_b, b[next - 1].untyped, inst.d.u);
        mfc_get(b[next - 1].untyped, ea_b, inst.d.u, next, 0, 0);
    }

    mfc_write_tag_mask(1 << current);
    mfc_read_tag_status_all();

    for (unsigned i(0) ; i < inst.size / sizeof(vector float) ; ++i)
    {
        intermediate.value = spu_madd(a[current - 1].vectorised[i], b[current - 1].vectorised[i], intermediate.value);
    }

    unsigned temp(next);
    next = current;
    current = temp;

    if (inst.e.u > 0)
    {
        mfc_write_tag_mask(1 << current);
        mfc_read_tag_status_all();

        for (unsigned i(0) ; i < inst.e.u / sizeof(vector float) ; ++i)
        {
            intermediate.value = spu_madd(a[current - 1].vectorised[i], b[current - 1].vectorised[i], intermediate.value);
        }
    }

    result.value = intermediate.array[0] + intermediate.array[1] + intermediate.array[2] + intermediate.array[3];

    allocator::release_block(*block_a[0]);
    allocator::release_block(*block_a[1]);
    allocator::release_block(*block_b[0]);
    allocator::release_block(*block_b[1]);

    return result.mail;
}
