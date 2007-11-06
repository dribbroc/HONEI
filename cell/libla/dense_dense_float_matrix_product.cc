/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
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
#include <cell/libutil/transfer.hh>

#include <spu_intrinsics.h>
#include <spu_mfcio.h>
#include <stdio.h>

using namespace honei;

void dense_dense_float_matrix_product(const Instruction & inst)
{
    //printf("dense_dense_float_matrix_product:\n");
    // inst.size = a.columns() = b.rows()
    // inst.d.u = a.rows()
    // inst.e.u = b.columns()

    allocator::Allocation * block_a(allocator::acquire_block());
    allocator::Allocation * block_b(allocator::acquire_block());
    allocator::Allocation * block_r(allocator::acquire_block());

    Pointer<float> a = { block_a->address };
    Pointer<float> b = { block_b->address };
    Pointer<float> r = { block_r->address };

    mfc_get(a.untyped, inst.a.ea, multiple_of_sixteen(inst.size * inst.d.u * sizeof(float)), 1, 0, 0);
    mfc_write_tag_mask(1 << 1);
    mfc_read_tag_status_any();
    mfc_get(b.untyped, inst.b.ea, multiple_of_sixteen(inst.size * inst.e.u * sizeof(float)), 2, 0, 0);
    mfc_write_tag_mask(1 << 2);
    mfc_read_tag_status_any();
    mfc_get(r.untyped, inst.c.ea, multiple_of_sixteen(inst.d.u * inst.e.u * sizeof(float)), 3, 0, 0);
    mfc_write_tag_mask(1 << 3);
    mfc_read_tag_status_any();

    unsigned long b_offset(inst.e.u % 4); // The number of elements BEHIND the last complete vector in every column.

    union {
        unsigned long v;
        float f;
    } u = {~0}; // im 2er-Kompl ist Negation von 0 = -1 = 111...111 (umweg ueber ~0 noetig da unsigned)

    const vector float mask_vector[4] = { // b_offset is the index into this vector
        { u.f, u.f, u.f, u.f },
        { u.f, 0x0000, 0x0000, 0x0000 },
        { u.f, u.f, 0x0000, 0x0000 },
        { u.f, u.f, u.f, 0x0000 }
    };

    // shuffle-offsets depend on b.columns() % 4 = b_offset and actual row in b % 4
    const unsigned int extract_offsets[4][4] = {  // [b_offset][b_row % 4]
        { 0, 0, 0, 0 }, // b_offset = 0, b_row % 4 = 0 - 3
        { 0, 1, 2, 3 }, // b_offset = 1, ...
        { 0, 2, 0, 2 },
        { 0, 3, 2, 1 }
    };

    // the indices-offset of the FIRST vector that contains elements of a row in b assuming size < 4
    const unsigned int vector_indices[4][4] = { // [b_offset][b_row % 4]
        { 0, 0, 0, 0 },
        { 0, 0, 0, 0 },
        { 0, 0, 1, 1 },
        { 0, 0, 1, 2 }
    };

    // the number of vectors in a row of b for b.columns() % 4 assuming size < 4
    const unsigned int vectors_per_row[4][4] = {  // [b_offset][b_row % 4]
        { 0, 0, 0, 0},
        { 1, 1, 1, 1},
        { 1, 1, 1, 1},
        { 1, 1, 1, 1}
    };

    unsigned long a_elem(0); // The actual considered element of matrix a

    for( ; a_elem < inst.d.u * inst.size ; a_elem++)
    {
        unsigned long act_a_row = a_elem / inst.size; // a_elem is in row a_elem / inst.size = a.columns()
        unsigned long act_b_row = a_elem % inst.size; // The row to multiply with is the index of a_elem modulo the number of rows of b

        unsigned long vecs_in_act_b_row = vectors_per_row[b_offset][act_b_row % 4] + (inst.e.u / 4); // Number of vectors in actual row of b
        unsigned long b_vec_idx = (act_b_row * (inst.e.u / 4)) + ((act_b_row / 4) * b_offset) + vector_indices[b_offset][act_b_row % 4];

        unsigned long r_vec_idx = (act_a_row * (inst.e.u / 4)) + ((act_a_row / 4) * b_offset) + vector_indices[b_offset][act_a_row % 4];

        for(unsigned i(0) ; i < vecs_in_act_b_row - 1 ; i++) // Make all computations except the last
        {
            //r_vec_idx += i;

            vector float temp = b.vectorised[b_vec_idx]; // temp version needed, cause original matrix must not be changed!
            Subscriptable<float> r_temp = { r.vectorised[r_vec_idx + i] }; // result matrix must either never be changed!
            extract(temp, b.vectorised[b_vec_idx+1], extract_offsets[b_offset][act_b_row % 4]);
            extract(r_temp.value, r.vectorised[r_vec_idx + i + 1], extract_offsets[b_offset][act_a_row % 4]);

            r_temp.value = spu_madd(spu_splats(a.typed[a_elem]), temp, r_temp.value);

            r.typed[(r_vec_idx + i) * 4 + extract_offsets[b_offset][act_a_row % 4]] = r_temp.array[0];
            r.typed[(r_vec_idx + i) * 4 + extract_offsets[b_offset][act_a_row % 4] + 1] = r_temp.array[1];
            r.typed[(r_vec_idx + i) * 4 + extract_offsets[b_offset][act_a_row % 4] + 2] = r_temp.array[2];
            r.typed[(r_vec_idx + i) * 4 + extract_offsets[b_offset][act_a_row % 4] + 3] = r_temp.array[3];

            b_vec_idx++;
        }
        // Now we handle the last vector in row where we perhaps must cut some elements that don't belong to the current row

        r_vec_idx += (vecs_in_act_b_row - 1); // Take vector for the last computation of the possibly incomplete vector

        vector float temp = b.vectorised[b_vec_idx];
        Subscriptable<float> r_temp = { r.vectorised[r_vec_idx] };
        extract(temp, b.vectorised[b_vec_idx+1], extract_offsets[b_offset][act_b_row % 4]);
        extract(r_temp.value, r.vectorised[r_vec_idx+1], extract_offsets[b_offset][act_a_row % 4]);

        Subscriptable<float> v = { spu_and(temp, mask_vector[b_offset]) };

        r_temp.value = spu_madd(spu_splats(a.typed[a_elem]), v.value, r_temp.value);

        r.typed[r_vec_idx * 4 + extract_offsets[b_offset][act_a_row % 4]] = r_temp.array[0];
        r.typed[r_vec_idx * 4 + extract_offsets[b_offset][act_a_row % 4] + 1] = r_temp.array[1];
        r.typed[r_vec_idx * 4 + extract_offsets[b_offset][act_a_row % 4] + 2] = r_temp.array[2];
        r.typed[r_vec_idx * 4 + extract_offsets[b_offset][act_a_row % 4] + 3] = r_temp.array[3];
    }

    mfc_put(r.untyped, inst.c.ea, multiple_of_sixteen(inst.d.u * inst.e.u * sizeof(float)), 4, 0, 0);
    mfc_write_tag_mask(1 << 4);
    mfc_read_tag_status_any();

    allocator::release_block(*block_a);
    allocator::release_block(*block_b);
    allocator::release_block(*block_r);

}
