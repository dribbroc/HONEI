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

using namespace honei::cell;

/*
 * dense_float_scale
 *
 * Scale a dense entity.
 *
 * \size Default transfer buffer size in bytes.
 * \operand a Base address of the entity.
 * \operand b Number of transfers needed.
 * \operand c Last transfer buffer size in bytes.
 * \operand d The scalar to use.
 */
void dense_float_scale(const Instruction & inst)
{
    EffectiveAddress ea_m(inst.a.ea), ea_r(inst.a.ea);

    Allocation * block_m[2] = { acquire_block(), acquire_block() };
    Pointer<float> m[2] = { { block_m[0]->address} , { block_m[1]->address } };

    unsigned counter(inst.b.u);
    unsigned size(counter > 1 ? inst.size : inst.c.u);
    unsigned nextsize;
    unsigned current(1), next(2);

    mfc_get(m[current - 1].untyped, ea_m, size, current, 0, 0);
    ea_m += size;

    Subscriptable<float> scale_vector = { spu_splats(static_cast<float>(inst.d.f)) };

    while (counter > 1)
    {
        nextsize = (counter == 2 ? inst.c.u : inst.size);

        mfc_get(m[next - 1].untyped, ea_m, nextsize, next, 0, 0);
        ea_m += nextsize;

        mfc_write_tag_mask(1 << current);
        mfc_read_tag_status_all();

        for (unsigned i(0); i < size / sizeof(vector float) ; ++i)
        {
            m[current - 1].vectorised[i] = spu_mul(m[current - 1].vectorised[i], scale_vector.value);
        }

        mfc_putb(m[current - 1].untyped, ea_r, size, current, 0, 0);
        ea_r += size;
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
        m[current - 1].vectorised[i] = spu_mul(m[current - 1].vectorised[i], scale_vector.value);
    }

    mfc_put(m[current - 1].untyped, ea_r, size, current, 0, 0);

    mfc_write_tag_mask(1 << current);
    mfc_read_tag_status_all();

    release_block(*block_m[0]);
    release_block(*block_m[1]);
}
