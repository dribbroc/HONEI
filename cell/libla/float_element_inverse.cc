/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Sven Mallach <sven.mallach@honei.org>
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
#include <cell/libutil/transfer.hh>

#include <spu_intrinsics.h>
#include <spu_mfcio.h>
#include <stdio.h>

using namespace honei::cell;

void float_element_inverse(const Instruction & inst)
{
    EffectiveAddress ea_a(inst.a.ea), ea_result(inst.a.ea);

    Allocation * block_a[2] = { acquire_block(), acquire_block() };

    Pointer<float> a[2] = { { block_a[0]->address} , { block_a[1]->address } };

    unsigned counter(inst.b.u);
    unsigned size(counter > 1 ? inst.size : inst.c.u);
    unsigned nextsize;
    unsigned current(1), next(2);

    debug_get(ea_a, a[current - 1].untyped, size);
    mfc_get(a[current - 1].untyped, ea_a, size, current, 0, 0);
    ea_a += size;

    while (counter > 1)
    {
        nextsize = (counter == 2 ? inst.c.u : inst.size);

        debug_get(ea_a, a[next - 1].untyped, nextsize);
        mfc_get(a[next - 1].untyped, ea_a, nextsize, next, 0, 0);
        ea_a += nextsize;

        mfc_write_tag_mask(1 << current);
        mfc_read_tag_status_all();

        for (unsigned i(0) ; i < size / sizeof(vector float) ; i++)
        {
            a[current - 1].vectorised[i] = spu_re(a[current - 1].vectorised[i]);
        }

        mfc_putb(a[current - 1].untyped, ea_result, size, current, 0, 0);
        ea_result += size;

        --counter;

        unsigned temp(next);
        next = current;
        current = temp;

        size = nextsize;
    }

    mfc_write_tag_mask(1 << current);
    mfc_read_tag_status_all();

    for (unsigned i(0) ; i < size / sizeof(vector float) ; i++)
    {
        a[current - 1].vectorised[i] = spu_re(a[current - 1].vectorised[i]);
    }

    mfc_putb(a[current - 1].untyped, ea_result, size, current, 0, 0);
    mfc_write_tag_mask(1 << current);
    mfc_read_tag_status_all();

    release_block(*block_a[0]);
    release_block(*block_a[1]);
}
