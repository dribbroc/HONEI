/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
 * product_banded_matrix_dense_vector_float
 *
 * \size Default transfer buffer size in bytes.
 * \operand a Base address of first entity.
 * \operand b Base address of second entity.
 * \operand c Base address of result entity.
 * \operand d Number of transfers needed.
 * \operand e Last transfer buffer size in bytes.
 * \operand f Band offset
 * \operand g SIMD x_offset
 */

void banded_dense_float_matrix_vector_product(const Instruction & inst)
{
    Allocation * block_a(acquire_block());
    Allocation * block_b(acquire_block());
    Allocation * block_r(acquire_block());

    Pointer<float> a = { block_a->address };
    Pointer<float> b = { block_b->address };
    Pointer<float> r = { block_r->address };

    unsigned x_offset(inst.g.u);
    unsigned y_offset((4 - x_offset) % 4);

    mfc_get(a.untyped, inst.a.ea, inst.size * sizeof(float), 1, 0, 0);
    mfc_get(b.untyped, inst.b.ea, (inst.size + x_offset + y_offset) * sizeof(float), 2, 0, 0);
    mfc_get(r.untyped, inst.c.ea, inst.size * sizeof(float), 3, 0, 0);
    mfc_write_tag_mask(1 << 3 | 1 << 2 | 1 << 1);
    mfc_read_tag_status_all();

    for (unsigned i(0) ; i < inst.size / 4 ; i++)
    {
        vector float temp = b.vectorised[i]; // temp version needed?
        extract(temp, b.vectorised[i + 1], x_offset);
        r.vectorised[i] = spu_madd(a.vectorised[i], temp, r.vectorised[i]);
    }

    mfc_put(r.untyped, inst.c.ea, inst.size * sizeof(float), 3, 0, 0);
    mfc_write_tag_mask(1 << 3);
    mfc_read_tag_status_all();

    release_block(*block_a);
    release_block(*block_b);
    release_block(*block_r);
}
