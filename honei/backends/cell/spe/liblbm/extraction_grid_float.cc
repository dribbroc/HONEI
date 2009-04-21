/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2009 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/util/attributes.hh>
#include <honei/backends/cell/cell.hh>
#include <honei/backends/cell/spe/libutil/allocator.hh>
#include <honei/backends/cell/spe/libutil/transfer.hh>

#include <spu_intrinsics.h>
#include <spu_mfcio.h>

#include <stdio.h>

/*
 * extraction_grid_float
 *
 * Calculate lbm grid extraction
 *
 * \size Default transfer buffer size in bytes.
 * \operand a Base address of real data.
*/

using namespace honei::cell;

void extraction_grid_float(const Instruction & instruction)
{
    Operand operands[32] __attribute__((aligned(16)));
    mfc_get(operands, instruction.a.ea, sizeof(Operand) * 32, 0, 0, 0);
    mfc_write_tag_mask(1 << 0);
    mfc_read_tag_status_any();

    EffectiveAddress ea_h(operands[0].ea), ea_u(operands[1].ea), ea_v(operands[2].ea),
                     ea_f_0(operands[3].ea), ea_f_1(operands[4].ea),
                     ea_f_2(operands[5].ea), ea_f_3(operands[6].ea),
                     ea_f_4(operands[7].ea), ea_f_5(operands[8].ea),
                     ea_f_6(operands[9].ea), ea_f_7(operands[10].ea),
                     ea_f_8(operands[11].ea);

    vector float dx_0(spu_splats(operands[12].fa[0]));
    vector float dx_1(spu_splats(operands[13].fa[0]));
    vector float dx_2(spu_splats(operands[14].fa[0]));
    vector float dx_3(spu_splats(operands[15].fa[0]));
    vector float dx_4(spu_splats(operands[16].fa[0]));
    vector float dx_5(spu_splats(operands[17].fa[0]));
    vector float dx_6(spu_splats(operands[18].fa[0]));
    vector float dx_7(spu_splats(operands[19].fa[0]));
    vector float dx_8(spu_splats(operands[20].fa[0]));
    vector float dy_0(spu_splats(operands[21].fa[1]));
    vector float dy_1(spu_splats(operands[22].fa[1]));
    vector float dy_2(spu_splats(operands[23].fa[1]));
    vector float dy_3(spu_splats(operands[24].fa[1]));
    vector float dy_4(spu_splats(operands[25].fa[1]));
    vector float dy_5(spu_splats(operands[26].fa[1]));
    vector float dy_6(spu_splats(operands[27].fa[1]));
    vector float dy_7(spu_splats(operands[28].fa[1]));
    vector float dy_8(spu_splats(operands[29].fa[1]));


    Allocation * block_h = { acquire_block() };
    Allocation * block_u = { acquire_block() };
    Allocation * block_v = { acquire_block() };
    Allocation * block_fa[2] = { acquire_block(), acquire_block() };
    Allocation * block_fb[2] = { acquire_block(), acquire_block() };
    Allocation * block_fc[2] = { acquire_block(), acquire_block() };

    Pointer<float> hp = { block_h->address };
    Pointer<float> up = { block_u->address };
    Pointer<float> vp = { block_v->address};
    Pointer<float> fa[2] = { { block_fa[0]->address} , { block_fa[1]->address } };
    Pointer<float> fb[2] = { { block_fb[0]->address} , { block_fb[1]->address } };
    Pointer<float> fc[2] = { { block_fc[0]->address} , { block_fc[1]->address } };

    unsigned counter(instruction.b.u);
    unsigned size(counter > 1 ? instruction.size : instruction.c.u);
    unsigned nextsize;
    unsigned current(0), next(1);

    //float scalar = instruction.d.f; // optional scalar value to be computed.

    // f1 f2 f3
    debug_get(ea_f_0, fa[current].untyped, size);
    mfc_get(fa[current].untyped, ea_f_0, size, current, 0, 0);
    ea_f_0 += size;

    debug_get(ea_f_1, fb[current].untyped, size);
    mfc_get(fb[current].untyped, ea_f_1, size, current, 0, 0);
    ea_f_1 += size;

    debug_get(ea_f_2, fc[current].untyped, size);
    mfc_get(fc[current].untyped, ea_f_2, size, current, 0, 0);
    ea_f_2 += size;

    //printf("%f %f \n", h[current].typed[5], h[current].typed[257]);
    while (counter > 1)
    {
        nextsize = (counter == 2 ? instruction.c.u : instruction.size);

        // f4 f5 f6
        debug_get(ea_f_3, fa[next].untyped, nextsize);
        mfc_get(fa[next].untyped, ea_f_3, nextsize, next, 0, 0);
        ea_f_3 += nextsize;

        debug_get(ea_f_4, fb[next].untyped, nextsize);
        mfc_get(fb[next].untyped, ea_f_4, nextsize, next, 0, 0);
        ea_f_4 += nextsize;

        debug_get(ea_f_5, fc[next].untyped, nextsize);
        mfc_get(fc[next].untyped, ea_f_5, nextsize, next, 0, 0);
        ea_f_5 += nextsize;

        mfc_write_tag_mask(1 << current);
        mfc_read_tag_status_all();

        vector float * h(hp.vectorised);
        vector float * u(up.vectorised);
        vector float * v(vp.vectorised);
        vector float * f_0(fa[current].vectorised);
        vector float * f_1(fb[current].vectorised);
        vector float * f_2(fc[current].vectorised);

        for (unsigned i(0) ; i < size / sizeof(vector float) ; ++i)
        {
            h[i] = spu_add(f_0[i], f_1[i]);
            h[i] = spu_add(h[i], f_2[i]);
            u[i] = spu_mul(f_0[i], dx_0);
            u[i] = spu_madd(f_1[i], dx_1, u[i]);
            u[i] = spu_madd(f_2[i], dx_2, u[i]);
            v[i] = spu_mul(f_0[i], dy_0);
            v[i] = spu_madd(f_1[i], dy_1, v[i]);
            v[i] = spu_madd(f_2[i], dy_2, v[i]);
        }

        unsigned temp(next);
        next = current;
        current = temp;

        // f6 f7 f8
        debug_get(ea_f_6, fa[next].untyped, nextsize);
        mfc_get(fa[next].untyped, ea_f_6, nextsize, next, 0, 0);
        ea_f_6 += nextsize;

        debug_get(ea_f_7, fb[next].untyped, nextsize);
        mfc_get(fb[next].untyped, ea_f_7, nextsize, next, 0, 0);
        ea_f_7 += nextsize;

        debug_get(ea_f_8, fc[next].untyped, nextsize);
        mfc_get(fc[next].untyped, ea_f_8, nextsize, next, 0, 0);
        ea_f_8 += nextsize;

        mfc_write_tag_mask(1 << current);
        mfc_read_tag_status_all();

        vector float * f_3(fa[current].vectorised);
        vector float * f_4(fb[current].vectorised);
        vector float * f_5(fc[current].vectorised);

        for (unsigned i(0) ; i < size / sizeof(vector float) ; ++i)
        {
            h[i] = spu_add(h[i], f_3[i]);
            h[i] = spu_add(h[i], f_4[i]);
            h[i] = spu_add(h[i], f_5[i]);
            u[i] = spu_madd(f_3[i], dx_3, u[i]);
            u[i] = spu_madd(f_4[i], dx_4, u[i]);
            u[i] = spu_madd(f_5[i], dx_5, u[i]);
            v[i] = spu_madd(f_3[i], dy_3, v[i]);
            v[i] = spu_madd(f_4[i], dy_4, v[i]);
            v[i] = spu_madd(f_5[i], dy_5, v[i]);
        }

        temp = next;
        next = current;
        current = temp;


        // prefetch data for next size iteration
        debug_get(ea_f_0, fa[next].untyped, nextsize);
        mfc_get(fa[next].untyped, ea_f_0, nextsize, next, 0, 0);
        ea_f_0 += nextsize;

        debug_get(ea_f_1, fb[next].untyped, nextsize);
        mfc_get(fb[next].untyped, ea_f_1, nextsize, next, 0, 0);
        ea_f_1 += nextsize;

        debug_get(ea_f_2, fc[next].untyped, nextsize);
        mfc_get(fc[next].untyped, ea_f_2, nextsize, next, 0, 0);
        ea_f_2 += nextsize;

        mfc_write_tag_mask(1 << current);
        mfc_read_tag_status_all();

        vector float * f_6(fa[current].vectorised);
        vector float * f_7(fb[current].vectorised);
        vector float * f_8(fc[current].vectorised);
        vector float vd, id, td;
        const vector float one = spu_splats(1.f);

        for (unsigned i(0) ; i < size / sizeof(vector float) ; ++i)
        {
            h[i] = spu_add(h[i], f_6[i]);
            h[i] = spu_add(h[i], f_7[i]);
            h[i] = spu_add(h[i], f_8[i]);
            u[i] = spu_madd(f_6[i], dx_6, u[i]);
            u[i] = spu_madd(f_7[i], dx_7, u[i]);
            u[i] = spu_madd(f_8[i], dx_8, u[i]);
            v[i] = spu_madd(f_6[i], dy_6, v[i]);
            v[i] = spu_madd(f_7[i], dy_7, v[i]);
            v[i] = spu_madd(f_8[i], dy_8, v[i]);

            vd = h[i];
            id = spu_re(vd);
            td = spu_nmsub(vd, id, one);
            vd = spu_madd(td, id, id);
            u[i] = spu_mul(u[i], vd);
            v[i] = spu_mul(v[i], vd);
        }


        mfc_putb(hp.untyped, ea_h, size, current, 0, 0);
        ea_h += size;

        mfc_putb(up.untyped, ea_u, size, current, 0, 0);
        ea_u += size;

        mfc_putb(vp.untyped, ea_v, size, current, 0, 0);
        ea_v += size;

        --counter;

        temp = next;
        next = current;
        current = temp;

        size = nextsize;
    }
    nextsize = size;
    mfc_write_tag_mask(1 << current);
    mfc_read_tag_status_all();

    // f4 f5 f6
    debug_get(ea_f_3, fa[next].untyped, nextsize);
    mfc_get(fa[next].untyped, ea_f_3, nextsize, next, 0, 0);
    ea_f_3 += nextsize;

    debug_get(ea_f_4, fb[next].untyped, nextsize);
    mfc_get(fb[next].untyped, ea_f_4, nextsize, next, 0, 0);
    ea_f_4 += nextsize;

    debug_get(ea_f_5, fc[next].untyped, nextsize);
    mfc_get(fc[next].untyped, ea_f_5, nextsize, next, 0, 0);
    ea_f_5 += nextsize;

    mfc_write_tag_mask(1 << current);
    mfc_read_tag_status_all();

    vector float * h(hp.vectorised);
    vector float * u(up.vectorised);
    vector float * v(vp.vectorised);
    vector float * f_0(fa[current].vectorised);
    vector float * f_1(fb[current].vectorised);
    vector float * f_2(fc[current].vectorised);

    for (unsigned i(0) ; i < size / sizeof(vector float) ; ++i)
    {
        h[i] = spu_add(f_0[i], f_1[i]);
        h[i] = spu_add(h[i], f_2[i]);
        u[i] = spu_mul(f_0[i], dx_0);
        u[i] = spu_madd(f_1[i], dx_1, u[i]);
        u[i] = spu_madd(f_2[i], dx_2, u[i]);
        v[i] = spu_mul(f_0[i], dy_0);
        v[i] = spu_madd(f_1[i], dy_1, v[i]);
        v[i] = spu_madd(f_2[i], dy_2, v[i]);
    }

    unsigned temp(next);
    next = current;
    current = temp;

    // f6 f7 f8
    debug_get(ea_f_6, fa[next].untyped, nextsize);
    mfc_get(fa[next].untyped, ea_f_6, nextsize, next, 0, 0);
    ea_f_6 += nextsize;

    debug_get(ea_f_7, fb[next].untyped, nextsize);
    mfc_get(fb[next].untyped, ea_f_7, nextsize, next, 0, 0);
    ea_f_7 += nextsize;

    debug_get(ea_f_8, fc[next].untyped, nextsize);
    mfc_get(fc[next].untyped, ea_f_8, nextsize, next, 0, 0);
    ea_f_8 += nextsize;

    mfc_write_tag_mask(1 << current);
    mfc_read_tag_status_all();

    vector float * f_3(fa[current].vectorised);
    vector float * f_4(fb[current].vectorised);
    vector float * f_5(fc[current].vectorised);

    for (unsigned i(0) ; i < size / sizeof(vector float) ; ++i)
    {
        h[i] = spu_add(h[i], f_3[i]);
        h[i] = spu_add(h[i], f_4[i]);
        h[i] = spu_add(h[i], f_5[i]);
        u[i] = spu_madd(f_3[i], dx_3, u[i]);
        u[i] = spu_madd(f_4[i], dx_4, u[i]);
        u[i] = spu_madd(f_5[i], dx_5, u[i]);
        v[i] = spu_madd(f_3[i], dy_3, v[i]);
        v[i] = spu_madd(f_4[i], dy_4, v[i]);
        v[i] = spu_madd(f_5[i], dy_5, v[i]);
    }

    temp = next;
    next = current;
    current = temp;

    mfc_write_tag_mask(1 << current);
    mfc_read_tag_status_all();

    vector float * f_6(fa[current].vectorised);
    vector float * f_7(fb[current].vectorised);
    vector float * f_8(fc[current].vectorised);
    vector float vd, id, td;
    const vector float one = spu_splats(1.f);

    for (unsigned i(0) ; i < size / sizeof(vector float) ; ++i)
    {
        h[i] = spu_add(h[i], f_6[i]);
        h[i] = spu_add(h[i], f_7[i]);
        h[i] = spu_add(h[i], f_8[i]);
        u[i] = spu_madd(f_6[i], dx_6, u[i]);
        u[i] = spu_madd(f_7[i], dx_7, u[i]);
        u[i] = spu_madd(f_8[i], dx_8, u[i]);
        v[i] = spu_madd(f_6[i], dy_6, v[i]);
        v[i] = spu_madd(f_7[i], dy_7, v[i]);
        v[i] = spu_madd(f_8[i], dy_8, v[i]);

        vd = h[i];
        id = spu_re(vd);
        td = spu_nmsub(vd, id, one);
        vd = spu_madd(td, id, id);
        u[i] = spu_mul(u[i], vd);
        v[i] = spu_mul(v[i], vd);
    }


    mfc_putb(hp.untyped, ea_h, size, current, 0, 0);
    ea_h += size;

    mfc_putb(up.untyped, ea_u, size, current, 0, 0);
    ea_u += size;

    mfc_putb(vp.untyped, ea_v, size, current, 0, 0);
    ea_v += size;

    mfc_write_tag_mask(1 << current);
    mfc_read_tag_status_all();

    release_all_blocks();
}
