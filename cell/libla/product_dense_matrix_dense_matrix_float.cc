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
 * dense_dense_float_matrix_product
 *
 * Calculate the product of two dense matrices.
 *
 * \size default transfer size for a
 * \operand a Base address of first entity.
 * \operand b Base address of second entity.
 * \operand c Base address of result entity.
 * \operand d a.columns()
 * \operand e b.columns()
 * \operand f Number of transfers needed for a and r
 * \operand g Number of transfers needed for b
 * \operand h Size of the last transfer for a
 * \operand i Size of the last transfer for b
 * \operand j Size of the last transfer for r
 * \operand k default transfer size for r
 */
void product_dense_matrix_dense_matrix_float(const Instruction & inst)
{
    EffectiveAddress ea_a(inst.a.ea), ea_b(inst.b.ea), ea_r(inst.c.ea);
    EffectiveAddress ea_r_start(inst.c.ea);

    Allocation * block_a[2] = { acquire_block(), acquire_block() };
    Allocation * block_r[2] = { acquire_block(), acquire_block() };

    Pointer<float> a[2] = { { block_a[0]->address }, { block_a[1]->address } };
    Pointer<float> r[2] = { { block_r[0]->address }, { block_r[1]->address } };

    Allocation * block_b[inst.g.u];
    Pointer<float> b[inst.g.u];
    // full transfers
    for (unsigned i(0) ; i < inst.g.u - 1 ; i++)
    {
        block_b[i] = acquire_block();
        b[i].untyped = block_b[i]->address;

        mfc_get(b[i].untyped, ea_b, 16384, 1, 0, 0);
        debug_get(ea_b, b[i].untyped, 16384);
        ea_b += 16384;
    }
    // last transfer
    block_b[inst.g.u - 1] = acquire_block();
    b[inst.g.u - 1].untyped = block_b[inst.g.u - 1]->address;
    mfc_get(b[inst.g.u - 1].untyped, ea_b, multiple_of_sixteen(inst.i.u), 1, 0, 0);
    debug_get(ea_b, b[inst.g.u - 1].untyped, multiple_of_sixteen(inst.i.u));

    Pointer<float> b_comp = b[0];

    unsigned last_a_size(multiple_of_sixteen(inst.h.u));
    unsigned last_r_size(multiple_of_sixteen(inst.j.u));
    unsigned counter(inst.f.u);
    unsigned size(counter > 1 ? inst.size : last_a_size);
    unsigned r_size(counter > 1 ? inst.k.u : last_r_size);
    unsigned nextsize, r_nextsize;
    unsigned current(1), next(2);

    mfc_get(a[current - 1].untyped, ea_a, size, current, 0, 0);
    debug_get(ea_a, a[current - 1].untyped, size);
    mfc_get(r[current - 1].untyped, ea_r, r_size, current, 0, 0);
    debug_get(ea_r, r[current - 1].untyped, r_size);

    ea_a += size;
    ea_r += r_size;

    unsigned long b_offset(inst.e.u % 4); // The number of elements BEHIND the last complete vector in every column.

    union {
        unsigned long v;
        float f;
    } u = {~0}; // in two's complement the negaton of 0 is -1 = 111...111 (way through ~0 necessary cause we're unsigned)

    const vector float mask_vector[4] = { // b_offset is the index into this vector
        { u.f, u.f, u.f, u.f },
        { u.f, 0x0000, 0x0000, 0x0000 },
        { u.f, u.f, 0x0000, 0x0000 },
        { u.f, u.f, u.f, 0x0000 }
    };

    // shuffle-offsets depend on b.columns() % 4 = b_offset and actual row in b % 4
    const unsigned extract_offsets[4][4] = {  // [b_offset][b_row % 4]
        { 0, 0, 0, 0 }, // b_offset = 0, b_row % 4 = 0 - 3
        { 0, 1, 2, 3 }, // b_offset = 1, ...
        { 0, 2, 0, 2 },
        { 0, 3, 2, 1 }
    };

    // the indices-offset of the FIRST vector that contains elements of a row in b assuming size < 4
    const unsigned vector_indices[4][4] = { // [b_offset][b_row % 4]
        { 0, 0, 0, 0 },
        { 0, 0, 0, 0 },
        { 0, 0, 1, 1 },
        { 0, 0, 1, 2 }
    };

    // the number of vectors in a row of b for b.columns() % 4 assuming size < 4
    const unsigned vectors_per_row[4] = { 0, 1, 1, 1}; // Use b_offset as index into this array.

    unsigned long a_elem(0); // The actual considered element of matrix a
    vector float a_elem_v;
    unsigned act_a_row, act_b_row, vecs_in_act_b_row, b_vec_idx, r_vec_idx;
    unsigned b_column_vecs(inst.e.u / 4);

    while (counter > 1)
    {
        a_elem = 0;
        nextsize = (counter == 2 ? last_a_size : inst.size);
        r_nextsize = (counter == 2 ? last_r_size : inst.k.u);

        debug_get(ea_a, a[next - 1].untyped, nextsize);
        mfc_get(a[next - 1].untyped, ea_a, nextsize, next, 0, 0);
        debug_get(ea_r, r[next - 1].untyped, r_nextsize);
        mfc_get(r[next - 1].untyped, ea_r, r_nextsize, next, 0, 0);
        ea_a += nextsize;
        ea_r += r_nextsize;

        mfc_write_tag_mask(1 << current);
        mfc_read_tag_status_all();

        for( ; a_elem < inst.size / sizeof(float) ; a_elem++)
        {
            a_elem_v = spu_splats(a[current - 1].typed[a_elem]);
            act_a_row = a_elem / inst.d.u; // a_elem is in row a_elem / inst.size = a.columns()
            act_b_row = a_elem % inst.d.u; // The row to multiply with is the index of a_elem modulo the number of rows of b

            vecs_in_act_b_row = vectors_per_row[b_offset] + b_column_vecs; // Number of vectors in actual row of b
            b_vec_idx = (act_b_row * b_column_vecs) + ((act_b_row / 4) * b_offset) + vector_indices[b_offset][act_b_row % 4];

            r_vec_idx = (act_a_row * b_column_vecs) + ((act_a_row / 4) * b_offset) + vector_indices[b_offset][act_a_row % 4];

            for(unsigned i(0) ; i < vecs_in_act_b_row - 1 ; i++) // Make all computations except the last
            {
                vector float b_temp(b_comp.vectorised[b_vec_idx]); // temp version needed, cause original matrix must not be changed!
                vector float r_temp(r[current - 1].vectorised[r_vec_idx + i]); // result matrix must either never be changed!
                extract(b_temp, b_comp.vectorised[b_vec_idx + 1], extract_offsets[b_offset][act_b_row % 4]);
                extract(r_temp, r[current - 1].vectorised[r_vec_idx + i + 1], extract_offsets[b_offset][act_a_row % 4]);

                r_temp = spu_madd(a_elem_v, b_temp, r_temp);

                insert(r[current - 1].vectorised[r_vec_idx + i], r[current - 1].vectorised[r_vec_idx + i + 1],
                        r_temp, extract_offsets[b_offset][act_a_row % 4]);

                b_vec_idx++;
            }
                // Now we handle the last vector in row where we perhaps must cut some elements that don't belong to the current row

                r_vec_idx += (vecs_in_act_b_row - 1); // Take vector for the last computation of the possibly incomplete vector

                vector float b_temp(b_comp.vectorised[b_vec_idx]);
                vector float r_temp(r[current - 1].vectorised[r_vec_idx]);
                extract(b_temp, b_comp.vectorised[b_vec_idx + 1], extract_offsets[b_offset][act_b_row % 4]);
                extract(r_temp, r[current - 1].vectorised[r_vec_idx + 1], extract_offsets[b_offset][act_a_row % 4]);

                vector float v(spu_and(b_temp, mask_vector[b_offset]));
                r_temp = spu_madd(a_elem_v, v, r_temp);

                insert(r[current -1].vectorised[r_vec_idx], r[current - 1].vectorised[r_vec_idx + 1], r_temp, extract_offsets[b_offset][act_a_row % 4]);
            }

        debug_put(ea_r_start, r[current - 1].untyped, r_size);
        mfc_putb(r[current - 1].untyped, ea_r_start, r_size, current, 0, 0);
        ea_r_start += r_size;
        --counter;

        unsigned temp(next);
        next = current;
        current = temp;
        size = nextsize;
        r_size = r_nextsize;
    }

    mfc_write_tag_mask(1 << current);
    mfc_read_tag_status_all();

    a_elem = 0; // The actual considered element of matrix a

    for( ; a_elem < inst.size / sizeof(float) ; a_elem++)
    {
        a_elem_v = spu_splats(a[current -1].typed[a_elem]);
        act_a_row = a_elem / inst.d.u; // a_elem is in row a_elem / inst.size = a.columns()
        act_b_row = a_elem % inst.d.u; // The row to multiply with is the index of a_elem modulo the number of rows of b

        vecs_in_act_b_row = vectors_per_row[b_offset] + b_column_vecs; // Number of vectors in actual row of b
        b_vec_idx = (act_b_row * b_column_vecs) + ((act_b_row / 4) * b_offset) + vector_indices[b_offset][act_b_row % 4];
        r_vec_idx = (act_a_row * b_column_vecs) + ((act_a_row / 4) * b_offset) + vector_indices[b_offset][act_a_row % 4];

        for(unsigned i(0) ; i < vecs_in_act_b_row - 1 ; i++) // Make all computations except the last
        {
            vector float b_temp(b_comp.vectorised[b_vec_idx]); // temp version needed, cause original matrix must not be changed!
            vector float r_temp(r[current - 1].vectorised[r_vec_idx + i]); // result matrix must either never be changed!
            extract(b_temp, b_comp.vectorised[b_vec_idx + 1], extract_offsets[b_offset][act_b_row % 4]);
            extract(r_temp, r[current - 1].vectorised[r_vec_idx + i + 1], extract_offsets[b_offset][act_a_row % 4]);

            r_temp = spu_madd(a_elem_v, b_temp, r_temp);

            insert(r[current - 1].vectorised[r_vec_idx + i], r[current - 1].vectorised[r_vec_idx + i + 1],
                    r_temp, extract_offsets[b_offset][act_a_row % 4]);

            b_vec_idx++;
        }
        // Now we handle the last vector in row where we perhaps must cut some elements that don't belong to the current row

        r_vec_idx += (vecs_in_act_b_row - 1); // Take vector for the last computation of the possibly incomplete vector

        vector float b_temp(b_comp.vectorised[b_vec_idx]);
        vector float r_temp(r[current - 1].vectorised[r_vec_idx]);
        extract(b_temp, b_comp.vectorised[b_vec_idx + 1], extract_offsets[b_offset][act_b_row % 4]);
        extract(r_temp, r[current - 1].vectorised[r_vec_idx + 1], extract_offsets[b_offset][act_a_row % 4]);

        vector float v(spu_and(b_temp, mask_vector[b_offset]));
        r_temp = spu_madd(a_elem_v, v, r_temp);

        insert(r[current - 1].vectorised[r_vec_idx], r[current - 1].vectorised[r_vec_idx + 1], r_temp, extract_offsets[b_offset][act_a_row % 4]);
    }

    debug_put(ea_r_start, r[current - 1].untyped, r_size);
    mfc_put(r[current - 1].untyped, ea_r_start, r_size, current, 0, 0);
    mfc_write_tag_mask(1 << current);
    mfc_read_tag_status_all();

    release_all_blocks();
}
