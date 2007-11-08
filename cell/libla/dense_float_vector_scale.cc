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

void dense_float_vector_scale(const Instruction & inst)
{
    // inst.size = the size of the vector
    // inst.b.f = the scalar

    allocator::Allocation * block_m(allocator::acquire_block());

    Pointer<float> m = { block_m->address };

    mfc_get(m.untyped, inst.a.ea, multiple_of_sixteen(inst.size * sizeof(float)), 1, 0, 0);
    mfc_write_tag_mask(1 << 1);
    mfc_read_tag_status_all();

    unsigned offset = inst.size % 4;
    vector float scale_vector(spu_splats(static_cast<float>(inst.b.f)));

    unsigned i(0);
    for ( ; i < inst.size / 4 ; i++)
    {
        m.vectorised[i] = spu_mul(m.vectorised[i], scale_vector);
    }

    for (unsigned j(0); j < offset ; j++)
    {
        m.typed[i * 4 + j] = m.typed[i * 4 + j] * inst.b.f;
    }

    mfc_put(m.untyped, inst.a.ea, multiple_of_sixteen(inst.size * sizeof(float)), 2, 0, 0);
    mfc_write_tag_mask(1 << 2);
    mfc_read_tag_status_all();

    allocator::release_block(*block_m);
}
