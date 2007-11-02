/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

unsigned dense_dense_float_dot_product(const Instruction & inst)
{
    printf("dense_dense_float_dot_product:\n");

    allocator::Allocation * block_a(allocator::acquire_block());
    allocator::Allocation * block_b(allocator::acquire_block());

    Pointer<float> a = { block_a->address };
    Pointer<float> b = { block_b->address };

    mfc_get(a.untyped, inst.a.ea, multiple_of_sixteen(inst.size * sizeof(float)), 1, 0, 0);
    mfc_get(b.untyped, inst.b.ea, multiple_of_sixteen(inst.size * sizeof(float)), 2, 0, 0);
    mfc_write_tag_mask(1 << 2 | 1 << 1);
    mfc_read_tag_status_all();

    Subscriptable<float> t = { spu_splats(0.0f) };

    unsigned i(0);
    for ( ; i < inst.size / 4 ; ++i)
    {
        t.value = spu_madd(a.vectorised[i], b.vectorised[i], t.value);
    }

    for (unsigned j(0) ; j < inst.size % 4 ; j += 1)
    {
        t.array[0] += a.typed[i * 4 + j] * b.typed[i * 4 + j];
    }

    t.array[0] += t.array[1] + t.array[2] + t.array[3];

    MailableResult<float> result;
    result.value = t.array[0];

    return result.mail;
}
