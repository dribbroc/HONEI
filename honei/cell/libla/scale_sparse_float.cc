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

#include <honei/attributes.hh>
#include <honei/cell/cell.hh>
#include <honei/cell/libutil/allocator.hh>
#include <honei/cell/libutil/transfer.hh>

#include <spu_intrinsics.h>
#include <spu_mfcio.h>

using namespace honei::cell;

/*
 * scale_sparse_float
 *
 * Scale a sparse entity.
 *
 * \size Default transfer buffer size in bytes for mfc_getl
 * \operand a Base address of the transfer list.
 * \operand b Effective address (high) of the transfers.
 * \operand c The scalar to use.
 * \operand d The accumulated size of the transfered elements.
 * \operand e The transfer buffer size in bytes for mfc_get
 */

void scale_sparse_float(const Instruction & inst)
{
    EffectiveAddress ea_list(inst.a.ea);
    vector float scalar = spu_splats(inst.c.f);

    ListElement list[inst.size] HONEI_ATTRIBUTE(aligned(8));

    debug_get(ea_list, list, sizeof(ListElement) * inst.e.u);
    mfc_get(list, ea_list, sizeof(ListElement) * inst.e.u, 0, 0, 0);
    mfc_write_tag_mask(1 << 0);
    mfc_read_tag_status_all();

    Allocation * block_data = acquire_block();
    Pointer<float> p = { block_data->address };

    debug_getl(inst.b.ea, block_data->address, inst.size * sizeof(ListElement));
    mfc_getl(block_data->address, inst.b.ea, &list[0], inst.size * sizeof(ListElement), 1, 0, 0);

    mfc_write_tag_mask(1 << 1);
    mfc_read_tag_status_all();

    unsigned i(0);
    for ( ; i < inst.d.u / sizeof(vector float) ; ++i)
    {
        p.vectorised[i] = spu_mul(p.vectorised[i], scalar);
    }

    debug_putl(inst.b.ea, block_data->address, inst.size * sizeof(ListElement));
    mfc_putl(block_data->address, inst.b.ea, &list[0], inst.size * sizeof(ListElement), 2, 0, 0);
    mfc_write_tag_mask(1 << 2);
    mfc_read_tag_status_all();

    release_block(*block_data);
}
