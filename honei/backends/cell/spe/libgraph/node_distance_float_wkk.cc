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

#include <honei/backends/cell/cell.hh>
#include <honei/backends/cell/spe/libutil/allocator.hh>
#include <honei/backends/cell/spe/libutil/transfer.hh>

#include <spu_intrinsics.h>
#include <spu_mfcio.h>

using namespace honei::cell;

/*
 * node_distance_float for Weighted Kamada-Kawai
 *
 * Calculate the distances of the nodes of a position matrix.
 *
 * \size Default transfer buffer size in bytes.
 * \operand a Base address of the block of the position matrix.
 * \operand b Base address of position matrix.
 * \operand c Transfer size in bytes for the block
 * \operand d Number of transfers needed (for "B")
 * \operand e Last transfer buffer size in bytes.
 * \operand f Base address of result part.
 * \operand g Number of rows of position matrix.
 *
 */

void node_distance_float_wkk(const Instruction & inst)
{
    EffectiveAddress ea_a(inst.a.ea), ea_b(inst.b.ea), ea_r(inst.f.u);

    Allocation * block_a = acquire_block();
    Allocation * block_b[2] = { acquire_block(), acquire_block() };

    Pointer<float> a = { block_a->address };
    Pointer<float> b[2] = { { block_b[0]->address }, { block_b[1]->address } };

    unsigned result_part_size(inst.g.u * ((inst.c.u / 4) / 2));

    Allocation * block_d = acquire_block();
    Pointer<float> d = { block_d->address };

    unsigned b_counter(inst.d.u);
    unsigned a_size(inst.c.u);
    unsigned b_size(b_counter > 1 ? inst.size : multiple_of_sixteen(inst.e.u));
    unsigned b_nextsize(0), b_next_calc_size(0);
    unsigned b_current(0), b_next(1);
    unsigned b_calc_size = b_counter > 1 ? inst.size : inst.e.u;

    mfc_get(a.untyped, ea_a, a_size, 0, 0, 0);
    mfc_get(b[b_current].untyped, ea_b, b_size, b_current, 0, 0);

    ea_b += b_size;

    unsigned ctr(0);

    while(b_counter > 1)
    {
        b_nextsize = (b_counter == 2 ? multiple_of_sixteen(inst.e.u) : inst.size);
        b_next_calc_size = (b_counter == 2 ? inst.e.u : inst.size);

        mfc_get(b[b_next].untyped, ea_b, b_nextsize, b_next, 0, 0);
        ea_b += b_nextsize;

        mfc_write_tag_mask(1 << b_current);
        mfc_read_tag_status_all();

        unsigned rest(b_calc_size % sizeof(vector float));

        for (unsigned i(0) ; i < a_size / sizeof(float) ; i += 2)
        {
            vector float xy = spu_splats(a.typed[i]);
            xy = spu_insert(a.typed[i + 1], xy, 1);
            xy = spu_insert(a.typed[i + 1], xy, 3);

            unsigned j(0);
            for ( ; j < b_calc_size / sizeof(vector float) ; j++)
            {
                Subscriptable<float> temp = { spu_sub(xy, b[b_current].vectorised[j]) };

                temp.value = spu_mul(temp.value, temp.value);

                d.typed[ctr] = temp.array[0] + temp.array[1];
                d.typed[ctr + 1] = temp.array[2] + temp.array[3];

                ctr++; ctr++;
            }
            for(unsigned k(0) ; k < rest / sizeof(float) ; k += 2)
            {
                float temp_x = a.typed[i] - b[b_current].typed[(j * 4) + k];
                float temp_y = a.typed[i + 1] - b[b_current].typed[(j * 4) + k + 1];

                temp_x *= temp_x;
                temp_y *= temp_y;

                d.typed[ctr] = temp_x + temp_y;
                ctr++;
            }
        }
        --b_counter;

        unsigned b_temp(b_next);
        b_next = b_current;
        b_current = b_temp;

        b_size = b_nextsize;
        b_calc_size = b_next_calc_size;

    }

    mfc_write_tag_mask(1 << b_current);
    mfc_read_tag_status_all();

    unsigned rest = b_calc_size % sizeof(vector float);

    for (unsigned i(0) ; i < a_size / sizeof(float) ; i += 2)
    {
        vector float xy = spu_splats(a.typed[i]);
        xy = spu_insert(a.typed[i + 1], xy, 1);
        xy = spu_insert(a.typed[i + 1], xy, 3);

        unsigned j(0);
        for ( ; j < b_calc_size / sizeof(vector float) ; j++)
        {
            Subscriptable<float> temp = { spu_sub(xy, b[b_current].vectorised[j]) };

            temp.value = spu_mul(temp.value, temp.value);

            d.typed[ctr] = temp.array[0] + temp.array[1];
            d.typed[ctr + 1] = temp.array[2] + temp.array[3];

            ctr++; ctr++;
        }
        for(unsigned k(0) ; k < rest / sizeof(float) ; k += 2)
        {
            float temp_x = a.typed[i] - b[b_current].typed[(j * 4) + k];
            float temp_y = a.typed[i + 1] - b[b_current].typed[(j * 4) + k + 1];

            temp_x *= temp_x;
            temp_y *= temp_y;

            d.typed[ctr] = temp_x + temp_y;
            ctr++;
        }
    }

    mfc_put(d.untyped, ea_r, (4 * result_part_size), 0, 0, 0);
    mfc_write_tag_mask(1 << 0);
    mfc_read_tag_status_all();

    release_all_blocks();
}
