/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Sven Mallach <sven.mallach@honei.org>
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
 * node_distance_float
 *
 * Calculate the distances of the nodes of a position matrix.
 *
 * \size Default transfer buffer size in bytes.
 * \operand a Base address of position matrix.
 * \operand b Number of transfers needed.
 * \operand c Last transfer buffer size in bytes.
 * \operand d value x.
 * \operand e value y.
 * \operand f Base address of result row to be computed (not double buffered)
 * \operand g size of result row.
 */

void node_distance_float(const Instruction & inst)
{
    EffectiveAddress ea_a(inst.a.ea), ea_b(inst.f.u);

    Allocation * block_a[2] = { acquire_block(), acquire_block() };
    Allocation * block_b = { acquire_block() };

    Pointer<float> a[2] = { { block_a[0]->address }, { block_a[1]->address } };
    Pointer<float> b = { block_b->address };

    unsigned counter(inst.b.u);
    unsigned size(counter > 1 ? inst.size : multiple_of_sixteen(inst.c.u));
    unsigned nextsize;
    unsigned current(0), next(1);

    vector float xy = spu_splats(inst.d.f);
    xy = spu_insert(inst.e.f, xy, 1);
    xy = spu_insert(inst.d.f, xy, 2);
    xy = spu_insert(inst.e.f, xy, 3);

    mfc_get(a[current].untyped, ea_a, size, current, 0, 0);
    mfc_get(b.untyped, ea_b, inst.g.u, current, 0, 0);
    ea_a += size;

    unsigned ctr(0);

    while (counter > 1)
    {
        nextsize = (counter == 2 ? multiple_of_sixteen(inst.c.u) : inst.size);

        mfc_get(a[next].untyped, ea_a, nextsize, next, 0, 0);
        ea_a += nextsize;

        mfc_write_tag_mask(1 << current);
        mfc_read_tag_status_all();

        for (unsigned i(0) ; i < size / sizeof(vector float) ; ++i)
        {
            Subscriptable <float> temp = { spu_sub(xy, a[current].vectorised[i]) };
            temp.value = spu_mul(temp.value, temp.value);
            b.typed[ctr] = temp.array[0] + temp.array[1];
            b.typed[ctr + 1] = temp.array[2] + temp.array[3];
            ctr += 2;
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
        Subscriptable <float> temp = { spu_sub(xy, a[current].vectorised[i]) };
        temp.value = spu_mul(temp.value, temp.value);
        b.typed[ctr] = temp.array[0] + temp.array[1];
        b.typed[ctr + 1] = temp.array[2] + temp.array[3];
        ctr += 2;
    }

    mfc_put(b.untyped, ea_b, inst.g.u, 0, 0, 0);

    release_all_blocks();

    mfc_write_tag_mask(1 << 0);
    mfc_read_tag_status_all();

}
