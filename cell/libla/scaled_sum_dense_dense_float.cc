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
#include <cell/libutil/debug.hh>
#include <cell/libutil/transfer.hh>

#include <spu_intrinsics.h>
#include <spu_mfcio.h>

using namespace honei::cell;

/*
 * scaled_sum_dense_dense_float
 *
 * Calculate the scaled sum of two dense entities.
 *
 * \size Default transfer buffer size in bytes.
 * \operand a Base address of first entity.
 * \operand b Base address of second entity.
 * \operand c Base address of the scalar.
 * \operand d Number of transfers needed.
 * \operand e Last transfer buffer size in bytes.
 * \operand f The alignment offset of second entity
 */
void scaled_sum_dense_dense_float(const Instruction & inst)
{
    EffectiveAddress ea_a(inst.a.ea), ea_b(inst.b.ea), ea_result(inst.a.ea);

    Allocation * block_a[2] = { acquire_block(), acquire_block() };
    Allocation * block_b[2] = { acquire_block(), acquire_block() };

    Pointer<float> a[2] = { { block_a[0]->address }, { block_a[1]->address } };
    Pointer<float> b[2] = { { block_b[0]->address }, { block_b[1]->address } };

    unsigned counter(inst.d.u);
    unsigned size(counter > 1 ? inst.size : inst.e.u);
    unsigned nextsize;
    unsigned current(1), next(2);

    unsigned b_offset(inst.f.u);

    vector float scalar_vector(spu_splats(inst.c.f));

    mfc_get(a[current - 1].untyped, ea_a, size, current, 0, 0);
    mfc_get(b[current - 1].untyped, ea_b, size, current, 0, 0);
    ea_a += size;
    ea_b += size;

    while (counter > 1)
    {
        nextsize = (counter == 2 ? inst.e.u : inst.size);

        ea_a -= 16;
        ea_b -= 16;
        mfc_get(a[next - 1].untyped, ea_a, nextsize, next, 0, 0);
        mfc_get(b[next - 1].untyped, ea_b, nextsize, next, 0, 0);
        ea_a += nextsize;
        ea_b += nextsize;

        mfc_write_tag_mask(1 << current);
        mfc_read_tag_status_all();

        for (unsigned i(0) ; i < (size - sizeof(vector float)) / sizeof(vector float) ; ++i)
        {
            extract(b[current - 1].vectorised[i], b[current - 1].vectorised[i + 1], b_offset);
            a[current - 1].vectorised[i] = spu_madd(b[current - 1].vectorised[i], scalar_vector, a[current - 1].vectorised[i]);
        }

        if (counter != inst.d.u)
            ea_result -= 16;

        mfc_putb(a[current - 1].untyped, ea_result, size, current, 0, 0);

        ea_result += size;

        --counter;

        unsigned temp(next);
        next = current;
        current = temp;

        size = nextsize;
    }

    mfc_write_tag_mask(1 << next); // Make sure earlier GETS/PUTS are done
    mfc_read_tag_status_all();

    // For calculation of "the end" (see below)
    mfc_get(a[next - 1].untyped, ea_a - 16, inst.d.u * 16, next, 0, 0);
    mfc_get(b[next - 1].untyped, ea_b - 16, (inst.d.u * 16) + 16, next, 0, 0);

    mfc_write_tag_mask(1 << current); // Make sure earlier GETS/PUTS are done
    mfc_read_tag_status_all();

    for (unsigned i(0) ; i < (size - sizeof(vector float)) / sizeof(vector float) ; ++i)
    {
        extract(b[current - 1].vectorised[i], b[current - 1].vectorised[i + 1], b_offset);
        a[current - 1].vectorised[i] = spu_madd(b[current - 1].vectorised[i], scalar_vector, a[current - 1].vectorised[i]);
    }

    if (counter != inst.d.u)
        ea_result -= 16;

    mfc_putb(a[current - 1].untyped, ea_result, size, current, 0, 0);

    // Now calculate the last vector of a and the shuffled last two vectors of b.
    mfc_write_tag_mask(1 << next); // Assure that GET(next) is done
    mfc_read_tag_status_all();

    for (unsigned i(0) ; i < inst.d.u ; i++)
    {
        extract(b[next - 1].vectorised[i], b[next - 1].vectorised[i + 1], b_offset);
        a[next - 1].vectorised[i] = spu_madd(b[next - 1].vectorised[i], scalar_vector, a[next - 1].vectorised[i]);
    }

    mfc_write_tag_mask(1 << current); // Assure that PUT(current) is done
    mfc_read_tag_status_all();

    mfc_put(a[next - 1].untyped, ea_result + size - 16, inst.d.u * 16, next, 0, 0);

    mfc_write_tag_mask(1 << next); // Assure that PUT(next) is done
    mfc_read_tag_status_all();

    release_all_blocks();
}
