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

#include <cell/libutil/allocator.hh>
#include <cell/cell.hh>

#include <spu_intrinsics.h>
#include <spu_mfcio.h>
#include <stdio.h>

using namespace honei;

void dense_dense_float_matrix_product(const Instruction & inst)
{
    printf("dense_dense_float_matrix_product:\n");
    // inst.size = a.rows() = b.columns()
    // inst.d.u = a.columns()
    // inst.e.u = b.rows()

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
    mfc_get(r.untyped, inst.c.ea, multiple_of_sixteen(inst.size * inst.size * sizeof(float)), 3, 0, 0);
    mfc_write_tag_mask(1 << 3);
    mfc_read_tag_status_any();

    printf("a.rows: %lu, a.columns: %llu, b.rows: %llu, b.columns: %u \n ", inst.size, inst.d.u, inst.e.u, inst.size);

    unsigned long a_elem(0); // The actual considered element of matrix a
    for (unsigned a_row(0) ; a_row < inst.size ; ++a_row) // Go through all rows of a
    {
        for(unsigned b_row(0) ; b_row < inst.e.u ; ++b_row) // Go through all rows of b
        {
            unsigned b_vecs(0);
            for( ; b_vecs < inst.size / 4 ; ++b_vecs) // Take all vectors of row b_row
            {
                r.vectorised[b_vecs + (a_row * inst.size / 4)] = spu_madd(spu_splats(a.typed[a_elem]), b.vectorised[b_vecs + (b_row * inst.size / 4)],
                        r.vectorised[b_vecs + (a_row * inst.size / 4)]); // Multiply-add a_elem with all elements of the current vector of matrix b.

                /* printf("Writing: a_row: %u, b_row: %u, b_vecs: %u, a_elem: %lu, vector in R: %lu, vector in B: %lu \n ", 
                        a_row, b_row, b_vecs, a_elem, a_row*inst.size/4 + b_vecs, b_row*inst.size/4 + b_vecs); */
            }

            for(unsigned b_vecs2(0) ; b_vecs2 < inst.size % 4 ; ++b_vecs2) // Take elements behind last complete vector of row b_row
            {
               unsigned index_a( (inst.size * a_row) + (b_vecs * 4) + b_vecs2 );
               unsigned index_b( (inst.size * b_row) + (b_vecs * 4) + b_vecs2 );
                //printf("Out of vector: a_elem: %lu, b_vecs: %u, b_vecs2: %u, index_r: %u, index_b: %u \n", a_elem, b_vecs, b_vecs2, index_a, index_b);
                r.typed[index_a] += a.typed[a_elem] *  b.typed[index_b]; // Multiply-add a_elem with all elements of the current vector of matrix b.
            }
            a_elem++; // After completing one row in b we can go on for the next element of matrix a.
        }
    }

    mfc_put(r.untyped, inst.c.ea, multiple_of_sixteen(inst.size * inst.size * sizeof(float)), 4, 0, 0);
    mfc_write_tag_mask(1 << 4);
    mfc_read_tag_status_any();
}
