/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Dirk Ribbrock <d_ribbrock@web.de>
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
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

#include <spu_intrinsics.h>
#include <spu_mfcio.h>
#include <stdio.h>

using namespace honei;

void dense_dense_float_matrix_vector_product(const Instruction & inst)
{
    printf("dense_dense_float_matrix_vector_product:\n");

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

    for (unsigned row(0) ; row < inst.d.u ; ++row)
    {
        printf("row: %u\n", row);
        Subscriptable<float> t = { spu_splats(0.0f) };
        for (unsigned i(0) ; i < inst.size / 4 ; ++i)
        {
            t.value = spu_madd(a.vectorised[i + row * inst.size / 4], x.vectorised[i], t.value);
        }

        t.array[0] += t.array[1] + t.array[2] + t.array[3];
        printf("value = %f\n", t.array[0]);
        r.typed[row] = t.array[0];
    }

    mfc_put(r.untyped, inst.c.ea, multiple_of_sixteen(inst.d.u * sizeof(float)), 3, 0, 0);
    mfc_write_tag_mask(1 << 3);
    mfc_read_tag_status_all();
}
