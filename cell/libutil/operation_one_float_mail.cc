/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Sven Mallach <sven.mallach@honei.org>
 *
 * This file is part of the Utility C++ library. LibUtil is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibUtil is distributed in the hope that it will be useful, but WITHOUT ANY
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
#include <cell/libutil/operations.hh>

#include <spu_intrinsics.h>
#include <spu_mfcio.h>

namespace honei
{
    namespace cell
    {
        template <>
        unsigned operation<Operation<1, float, rtm_mail> >(const Operation<1, float, rtm_mail> & operation,
                const Instruction & instruction)
        {
            EffectiveAddress ea_a(instruction.b.ea);

            Allocation * block_a[2] = { acquire_block(), acquire_block() };

            Pointer<float> a[2] = { { block_a[0]->address} , { block_a[1]->address } };

            unsigned counter(instruction.c.u);
            unsigned size(counter > 1 ? instruction.size : instruction.d.u);
            unsigned nextsize;
            unsigned current(0), next(1);

            debug_get(ea_a, a[current].untyped, size);
            mfc_get(a[current].untyped, ea_a, size, current, 0, 0);
            ea_a += size;

            Subscriptable<float> acc = { operation.init() };

            while (counter > 1)
            {
                nextsize = (counter == 2 ? instruction.d.u : instruction.size);

                debug_get(ea_a, a[next].untyped, nextsize);
                mfc_get(a[next].untyped, ea_a, nextsize, next, 0, 0);
                ea_a += nextsize;

                mfc_write_tag_mask(1 << current);
                mfc_read_tag_status_all();

                acc.value = operation.calculate(acc.value, a[current].vectorised, size / sizeof(vector float));

                --counter;

                unsigned temp(next);
                next = current;
                current = temp;

                size = nextsize;
            }

            mfc_write_tag_mask(1 << current);
            mfc_read_tag_status_all();

            acc.value = operation.calculate(acc.value, a[current].vectorised, size / sizeof(vector float));

            release_block(*block_a[0]);
            release_block(*block_a[1]);

            MailableResult<float> result = { operation.finish(acc.value) };

            return result.mail;
        }
    }
}
