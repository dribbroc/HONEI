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

void dense_dense_float_matrix_vector_product(const Instruction & inst)
{
    //printf("dense_dense_float_matrix_vector_product:\n");

    // inst.size = vector.size = matrix.columns
    // inst.d.u = matrix.rows

    allocator::Allocation * block_a(allocator::acquire_block());
    allocator::Allocation * block_x(allocator::acquire_block());
    allocator::Allocation * block_r(allocator::acquire_block());

    Pointer<float> a = { block_a->address };
    Pointer<float> x = { block_x->address };
    Pointer<float> r = { block_r->address };

    mfc_get(a.untyped, inst.a.ea, multiple_of_sixteen(inst.size * inst.d.u * sizeof(float)), 1, 0, 0);
    mfc_get(x.untyped, inst.b.ea, multiple_of_sixteen(inst.size * sizeof(float)), 2, 0, 0);
    mfc_write_tag_mask(1 << 2 | 1 << 1);
    mfc_read_tag_status_all();

    unsigned offset = inst.size % 4;
    union {
        unsigned long v;
        float f;
    } u = {~0}; // im 2er-Kompl ist Negation von 0 = -1 = 111...111 (umweg ueber ~0 noetig da unsigned)

    const vector float mask_vector[4] = { // offset is the index into this vector
        { u.f, u.f, u.f, u.f },
        { u.f, 0x0000, 0x0000, 0x0000 },
        { u.f, u.f, 0x0000, 0x0000 },
        { u.f, u.f, u.f, 0x0000 }
    };

    // shuffle-offsets depend on a.columns() % 4 = offset and actual row in a % 4
    const unsigned int extract_offsets[4][4] = {  // [offset][a_row % 4]
        { 0, 0, 0, 0 }, // offset = 0, a_row % 4 = 0 - 3
        { 0, 1, 2, 3 }, // offset = 1, ...
        { 0, 2, 0, 2 },
        { 0, 3, 2, 1 }
    };

    // the indices-offset of the FIRST vector that contains elements of a row in a assuming size < 4
    const unsigned int vector_indices[4][4] = { // [offset][a_row % 4]
        { 0, 0, 0, 0 },
        { 0, 0, 0, 0 },
        { 0, 0, 1, 1 },
        { 0, 0, 1, 2 }
    };

    // the number of vectors in a row of a for a.columns() % 4 assuming size < 4
    const unsigned int vectors_per_row[4][4] = {  // [offset][a_row % 4]
        { 0, 0, 0, 0},
        { 1, 1, 1, 1},
        { 1, 1, 1, 1},
        { 1, 1, 1, 1}
    };

    unsigned long r_elem(0); // The actual considered element of the result vector
    unsigned long a_row(0);
    for( ; r_elem < inst.d.u ; r_elem++)
    {
        unsigned long vecs_in_a_row = vectors_per_row[offset][a_row % 4] + (inst.size / 4); // Number of vectors in actual row of a
        unsigned long a_vec_idx = (a_row * (inst.size / 4)) + ((a_row / 4) * offset) + vector_indices[offset][a_row % 4];
        Subscriptable<float> res = { spu_splats(0.0f) };
        for(unsigned i(0) ; i < vecs_in_a_row - 1 ; i++) // Make all computations except the last
        {
            vector float temp = a.vectorised[a_vec_idx + i]; // temp version needed, cause original matrix must not be changed!
            extract(temp, a.vectorised[a_vec_idx + i + 1], extract_offsets[offset][a_row % 4]);

            res.value = spu_madd(temp, x.vectorised[i], res.value);
            Subscriptable<float> ausgabe = { temp };
            Subscriptable<float> ausgabe2 = { x.vectorised[i] };
       }
        // Now we handle the last vector in row where we perhaps must cut some elements that don't belong to the current row
        unsigned i = (vecs_in_a_row - 1); // Take vector for the last computation of the possibly incomplete vector

        vector float temp = a.vectorised[a_vec_idx + i];
        extract(temp, a.vectorised[a_vec_idx + i + 1], extract_offsets[offset][a_row % 4]);

        Subscriptable<float> v = { spu_and(x.vectorised[i], mask_vector[offset]) };

        res.value = spu_madd(temp, v.value, res.value);
        Subscriptable<float> ausgabe = { temp };

        // Finished computation for one row of a -> one element of r: Cleaning up
        r.typed[r_elem] = res.array[0] + res.array[1] + res.array[2] + res.array[3];
        a_row++;
    }

    mfc_put(r.untyped, inst.c.ea, multiple_of_sixteen(inst.d.u * sizeof(float)), 3, 0, 0);
    mfc_write_tag_mask(1 << 3);
    mfc_read_tag_status_all();
}
