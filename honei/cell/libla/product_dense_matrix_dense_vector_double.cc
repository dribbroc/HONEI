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
 * dense_dense_double_matrix_vector_product
 *
 * Calculate the product of a dense matrix and
 * a dense vector.
 *
 * \size The default transfer size for the matrix (representing full rows),
 * \operand a Base address of the matrix.
 * \operand b Base address of vector to multiply with.
 * \operand c Base address of result vector.
 * \operand d The number of transfers for the matrix.
 * \operand e The last transfer size for the matrix.
 * \operand f The number of transfers for the vector.
 * \operand g The last transfer size for the vector.
 * \operand h The number of columns of the matrix.
 * \operand i matrix.columns() modulo 4
 */

void product_dense_matrix_dense_vector_double(const Instruction & inst)
{
    EffectiveAddress ea_a(inst.a.ea), ea_x(inst.b.ea), ea_r(inst.c.ea);

    Allocation * block_a[2] = { acquire_block(), acquire_block() };
    Allocation * block_x[2] = { acquire_block(), acquire_block() };
    Allocation * block_r = acquire_block();

    Pointer<double> a[2] = { { block_a[0]->address }, { block_a[1]->address } };
    Pointer<double> x[2] = { { block_x[0]->address }, { block_x[1]->address } };
    Pointer<double> r = { block_r->address};

    unsigned a_cols(inst.h.u);
    unsigned a_counter(inst.d.u);
    unsigned x_counter(inst.f.u);
    unsigned a_size(a_counter > 1 ? inst.size : multiple_of_sixteen(inst.e.u));
    unsigned x_size(x_counter > 1 ? 16384 : multiple_of_sixteen(inst.g.u));
    unsigned r_size(((a_counter - 1) * (inst.size / a_cols)) + inst.e.u / a_cols);
    unsigned a_nextsize, x_nextsize;
    unsigned a_current(0), a_next(1);
    unsigned x_current(0), x_next(1);
    unsigned a_offset(inst.i.u);
    unsigned act_offset(0);

    debug_get(inst.a.ea, a[a_current].untyped, a_size);
    mfc_get(a[a_current].untyped, inst.a.ea, a_size, 1, 0, 0);
    debug_get(inst.b.ea, x[x_current].untyped, x_size);
    mfc_get(x[x_current].untyped, inst.b.ea, x_size, 1, 0, 0);
    mfc_write_tag_mask(1 << 1);
    mfc_read_tag_status_all();

    ea_a += a_size;
    ea_x += x_size;

    unsigned a_vecs(a_cols / 2);
    unsigned a_idx(0);
    unsigned r_idx(0);

    while (a_counter > 1)
    {
        a_nextsize = (a_counter == 2 ? multiple_of_sixteen(inst.e.u) : a_size);

        debug_get(ea_a, a[a_next].untyped, a_nextsize);
        mfc_get(a[a_next].untyped, ea_a, a_nextsize, a_next, 0, 0);
        ea_a += a_nextsize;

        x_counter = inst.f.u;
        a_idx = 0;

        mfc_write_tag_mask(1 << a_current);
        mfc_read_tag_status_all();

        while (x_counter >= 1)
        {
            if (x_counter == 2)
            {
                x_nextsize = multiple_of_sixteen(inst.g.u);
            }
            else
            {
                x_nextsize = x_size;
            }

            if (x_counter == 1)
            {
                ea_x = inst.b.ea;
                x_nextsize = inst.f.u > 1 ? 16384 : multiple_of_sixteen(inst.g.u);
            }

            debug_get(ea_x, x[x_next].untyped, x_nextsize);
            mfc_get(x[x_next].untyped, ea_x, x_nextsize, x_next, 0, 0);
            ea_x += x_nextsize;

            mfc_write_tag_mask(1 << x_current);
            mfc_read_tag_status_all();

            for (unsigned i(0) ; i < (a_size / a_cols / 8) ; i++) // for rows
            {
                r.typed[r_idx] = double(0);
                Subscriptable<double> temp = { spu_splats(r.typed[r_idx]) };

                for (unsigned j(0) ; j < a_vecs ; j++, a_idx++) // for vecs in the row
                {
                    vector double a_cur = a[a_current].vectorised[a_idx];
                    extract(a_cur, a[a_current].vectorised[a_idx + 1], act_offset);
                    temp.value = spu_madd(a_cur, x[x_current].vectorised[j], temp.value);
                }

                r.typed[r_idx] = temp.array[0] + temp.array[1];

                act_offset = (act_offset + a_offset) % 2;

                for (unsigned j(0) ; j < a_offset ; j++)
                {
                    r.typed[r_idx] += a[a_current].typed[(2 * a_idx) + j] * x[x_current].typed[(2 * a_vecs) + j];
                }

                if ( i != 0 && i % 2 == 0 && act_offset > 0)
                {
                    a_idx++;
                }

                r_idx++;
            }

            unsigned x_temp(x_current);
            x_current = x_next;
            x_next = x_temp;

            x_counter--;
        }

        unsigned a_temp(a_current);
        a_current = a_next;
        a_next = a_temp;
        a_counter--;
    }

    mfc_write_tag_mask(1 << a_current);
    mfc_read_tag_status_all();
    x_counter = inst.f.u;
    a_idx = 0;

    while (x_counter >= 1)
    {
        if (x_counter == 2)
        {
            x_nextsize = multiple_of_sixteen(inst.g.u);
        }
        else
        {
            x_nextsize = x_size;
        }

        debug_get(ea_x, x[x_next].untyped, x_nextsize);
        mfc_get(x[x_next].untyped, ea_x, x_nextsize, x_next, 0, 0);
        ea_x += x_nextsize;

        mfc_write_tag_mask(1 << x_current);
        mfc_read_tag_status_all();

        for (unsigned i(0) ; i < (a_size / a_cols / 8) ; i++) // for rows
        {
            r.typed[r_idx] = double(0);
            Subscriptable<double> temp = { spu_splats(r.typed[r_idx]) };

            for (unsigned j(0) ; j < a_vecs ; j++, a_idx++) // for vecs in the row
            {
                vector double a_cur = a[a_current].vectorised[a_idx];
                extract(a_cur, a[a_current].vectorised[a_idx + 1], act_offset);
                temp.value = spu_madd(a_cur, x[x_current].vectorised[j], temp.value);
            }

            r.typed[r_idx] = temp.array[0] + temp.array[1];
            act_offset = (act_offset + a_offset) % 2;

            for (unsigned j(0) ; j < a_offset ; j++)
            {
                r.typed[r_idx] += a[a_current].typed[(2 * a_idx) + j] * x[x_current].typed[(2 * a_vecs) + j];
            }

            if ( i != 0 && i % 2 == 0 && act_offset > 0)
            {
                a_idx++;
            }

            r_idx++;
        }

         unsigned x_temp(x_current);
         x_current = x_next;
         x_next = x_temp;
         x_counter--;
    }

    debug_put(ea_r, r.untyped, r_size);
    mfc_put(r.untyped, ea_r, r_size, 2, 0, 0);
    mfc_write_tag_mask(1 << 2);
    mfc_read_tag_status_all();

    release_all_blocks();
}
