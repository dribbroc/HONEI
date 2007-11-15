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

    allocator::Allocation * block_a[2] = { allocator::acquire_block(), allocator::acquire_block() };
    allocator::Allocation * block_b[2] = { allocator::acquire_block(), allocator::acquire_block() };

    Pointer<float> a[2] = { block_a[0]->address, block_a[1]->address };
    Pointer<float> b[2] = { block_b[0]->address, block_b[1]->address };

    unsigned counter(inst.d.u);
    unsigned size(counter > 1 ? inst.size : inst.e.u);
    unsigned nextsize;
    unsigned current(0), next(1);

    mfc_get(a[current].untyped, ea_a, size, current, 0, 0);
    mfc_get(b[current].untyped, ea_b, size, current, 0, 0);
    ea_a += size;
    ea_b += size;

    Subscriptable<float> acc = { spu_splats(0.0f) };

    while (counter > 1)
    {
        nextsize = (counter == 2 ? inst.e.u : inst.size);

        mfc_get(a[next].untyped, ea_a, nextsize, next, 0, 0);
        mfc_get(b[next].untyped, ea_b, nextsize, next, 0, 0);
        ea_a += nextsize;
        ea_b += nextsize;

        mfc_write_tag_mask(1 << current);
        mfc_read_tag_status_all();

        for (unsigned i(0) ; i < size / sizeof(vector float) ; ++i)
        {
            acc.value = spu_madd(a[current].vectorised[i], b[current].vectorised[i], acc.value);
        }

        --counter;

        unsigned temp(next);
        next = current;
        current = temp;

        size = nextsize;
    }

    mfc_write_tag_mask(1 << current);
    mfc_read_tag_status_all();

    for (unsigned i(0) ; i < size / sizeof(vector float) ; ++i)
    {
        acc.value = spu_madd(a[current].vectorised[i], b[current].vectorised[i], acc.value);
    }

    allocator::release_block(*block_a[0]);
    allocator::release_block(*block_a[1]);
    allocator::release_block(*block_b[0]);
    allocator::release_block(*block_b[1]);

    MailableResult<float> result = { acc.array[0] + acc.array[1] + acc.array[2] + acc.array[3] };

    return result.mail;
}
