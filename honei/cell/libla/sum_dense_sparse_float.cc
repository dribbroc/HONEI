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

#include <honei/cell/cell.hh>
#include <honei/cell/libutil/allocator.hh>
#include <honei/cell/libutil/transfer.hh>

#include <spu_intrinsics.h>
#include <spu_mfcio.h>

using namespace honei::cell;

void sum_dense_sparse_float(const Instruction & inst)
{
    EffectiveAddress ea_a(inst.a.ea), ea_b(inst.b.ea), ea_c(inst.c.ea), ea_result(inst.a.ea);

    Allocation * block_a[2] = { acquire_block(), acquire_block() }; // a's elements
    Allocation * block_b[2] = { acquire_block(), acquire_block() }; // b's elements
    Allocation * block_c[2] = { acquire_block(), acquire_block() }; // b's indices

    Pointer<float> a[2] = { { block_a[0]->address }, { block_a[1]->address } };
    Pointer<float> b[2] = { { block_b[0]->address }, { block_b[1]->address } };
    Pointer<unsigned long long> c[2] = { { block_c[0]->address }, { block_c[1]->address } };

    unsigned actual_last_index_in_dv(inst.i.u);
    unsigned dv_already_finished_elements(0); // for getting the right index in current part of dv
    unsigned a_counter(inst.e.u);
    unsigned b_counter(inst.g.u);

    unsigned a_size(a_counter > 1 ? inst.size : multiple_of_sixteen(inst.f.u));
    unsigned b_size(b_counter > 1 ? inst.size/2 : multiple_of_sixteen(inst.h.u));
    unsigned a_nextsize(0), b_nextsize(0);
    unsigned a_current(1), b_current(1), a_next(2), b_next(2);

    mfc_get(a[a_current - 1].untyped, ea_a, a_size, a_current, 0, 0);
    debug_get(ea_a, a[a_current -1].untyped, a_size);
    mfc_get(b[b_current - 1].untyped, ea_b, b_size, b_current, 0, 0);
    debug_get(ea_b, b[b_current - 1].untyped, b_size);
    mfc_get(c[b_current - 1].untyped, ea_c, b_size*2, 3, 0, 0);
    debug_get(ea_c, c[b_current - 1].untyped, b_size*2);

    ea_a += a_size;
    ea_b += b_size;
    ea_c += b_size*2;

    if (a_counter > 1)
    {
        a_nextsize = (a_counter == 2 ? multiple_of_sixteen(inst.f.u) : inst.size);
        mfc_get(a[a_next - 1].untyped, ea_a, a_nextsize, a_next, 0, 0);
        debug_get(ea_a, a[a_next - 1].untyped, a_nextsize);
        ea_a += a_nextsize;
        a_counter--;
    }

    while (b_counter > 1)
    {
        b_nextsize = (b_counter == 2 ? multiple_of_sixteen(inst.h.u) : inst.size/2);

        mfc_get(b[b_next - 1].untyped, ea_b, b_nextsize, b_next, 0, 0);
        debug_get(ea_b, b[b_next -1].untyped, b_nextsize);
        ea_b += b_nextsize;
        mfc_get(c[b_next - 1].untyped, ea_c, b_nextsize*2, b_next, 0, 0);
        debug_get(ea_c, c[b_next -1].untyped, b_nextsize*2);
        ea_c += b_nextsize*2;

        mfc_write_tag_mask(1 << b_current);
        mfc_read_tag_status_all();

        for (unsigned long long i(0); i < b_size / 4 ; i++)
        {
            if (c[b_current - 1].typed[i] > actual_last_index_in_dv)
            {
                if (c[b_current - 1].typed[i] > inst.j.u-1)
                    break;

                debug_put(ea_result, a[a_current -1].untyped, a_size);
                mfc_putb(a[a_current - 1].untyped, ea_result, a_size, a_current, 0, 0);
                ea_result += a_size;

                unsigned temp(a_next);
                a_next = a_current;
                a_current = temp;
                actual_last_index_in_dv += a_nextsize / 4;
                dv_already_finished_elements += a_size / 4;
                a_nextsize = (a_counter == 1 ? multiple_of_sixteen(inst.f.u) : inst.size);

                if (a_counter > 1)
                {
                    debug_get(ea_a, a[a_next -1].untyped, a_nextsize);
                    mfc_get(a[a_next - 1].untyped, ea_a, a_nextsize, a_next, 0, 0);
                    ea_a += a_nextsize;

                    a_counter--;
                }

                a_size = a_nextsize;

                mfc_write_tag_mask (1 << a_current);
                mfc_read_tag_status_all();
            }
            a[a_current - 1].typed[c[b_current - 1].typed[i] - dv_already_finished_elements] += b[b_current - 1].typed[i];
        }

        --b_counter;

        unsigned temp(b_next);
        b_next = b_current;
        b_current = temp;

        b_size = b_nextsize;
    }

    mfc_write_tag_mask(1 << b_current);
    mfc_read_tag_status_all();

    mfc_write_tag_mask(1 << 3);
    mfc_read_tag_status_any();

        for (unsigned long long i(0); i < b_size / 4 ; i++)
        {
            if (c[b_current - 1].typed[i] > actual_last_index_in_dv)
            {
                if (c[b_current - 1].typed[i] > inst.j.u-1)
                    break;

                debug_put(ea_result, a[a_current -1].untyped, a_size);
                mfc_putb(a[a_current - 1].untyped, ea_result, a_size, a_current, 0, 0);
                ea_result += a_size;

                unsigned temp(a_next);
                a_next = a_current;
                a_current = temp;
                actual_last_index_in_dv += a_nextsize / 4;
                dv_already_finished_elements += a_size / 4;
                a_nextsize = (a_counter == 1 ? multiple_of_sixteen(inst.f.u) : inst.size);

                if (a_counter > 1)
                {
                    debug_get(ea_a, a[a_next -1].untyped, a_nextsize);
                    mfc_get(a[a_next - 1].untyped, ea_a, a_nextsize, a_next, 0, 0);
                    ea_a += a_nextsize;

                    a_counter--;
                }

                a_size = a_nextsize;

                mfc_write_tag_mask (1 << a_current);
                mfc_read_tag_status_all();
            }
            a[a_current - 1].typed[c[b_current - 1].typed[i] - dv_already_finished_elements] += b[b_current - 1].typed[i];
        }

    mfc_put(a[a_current - 1].untyped, ea_result, a_size, a_current, 0, 0);
    debug_put(ea_result, a[a_current -1].untyped, a_size);

    mfc_write_tag_mask(1 << a_current);
    mfc_read_tag_status_all();

    release_all_blocks();
}
