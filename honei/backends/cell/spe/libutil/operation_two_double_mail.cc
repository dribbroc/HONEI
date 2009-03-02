/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2008 Sven Mallach <sven.mallach@honei.org>
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/backends/cell/cell.hh>
#include <honei/backends/cell/spe/libutil/allocator.hh>
#include <honei/backends/cell/spe/libutil/transfer.hh>
#include <honei/backends/cell/spe/libutil/operation_framework.hh>

#include <spu_intrinsics.h>
#include <spu_mfcio.h>

namespace honei
{
    namespace cell
    {
        template <>
        unsigned operation<Operation<2, double, rtm_mail> >(const Operation<2, double, rtm_mail> & operation,
                const Instruction & instruction)
        {
            EffectiveAddress ea_a(instruction.b.ea), ea_b(instruction.c.ea);

            Allocation * block_a[2] = { acquire_block(), acquire_block() };
            Allocation * block_b[2] = { acquire_block(), acquire_block() };

            Pointer<double> a[2] = { { block_a[0]->address }, { block_a[1]->address } };
            Pointer<double> b[2] = { { block_b[0]->address }, { block_b[1]->address } };

            unsigned counter(instruction.d.u);
            unsigned size(counter > 1 ? instruction.size : instruction.e.u);
            unsigned nextsize;
            unsigned current(0), next(1);

            double scalar1 = instruction.i.d;
            double scalar2 = instruction.j.d;

            unsigned offset(instruction.f.u);

            mfc_get(a[current].untyped, ea_a, size, current, 0, 0);
            mfc_get(b[current].untyped, ea_b, size, current, 0, 0);
            ea_a += size;
            ea_b += size;

            vector double accumulator(operation.init());
            vector double carry = { instruction.g.d, instruction.h.d };

            while (counter > 1)
            {
                nextsize = (counter == 2 ? instruction.e.u : instruction.size);

                mfc_get(a[next].untyped, ea_a, nextsize, next, 0, 0);
                mfc_get(b[next].untyped, ea_b, nextsize, next, 0, 0);
                ea_a += nextsize;
                ea_b += nextsize;

                mfc_write_tag_mask(1 << current);

                mfc_read_tag_status_all();

                accumulator = operation.calculate(accumulator, carry, a[current].vectorised, b[current].vectorised,
                        size / sizeof(vector double), offset, scalar1, scalar2);

                --counter;

                unsigned temp(next);
                next = current;
                current = temp;

                size = nextsize;
            }

            mfc_write_tag_mask(1 << current);
            mfc_read_tag_status_all();

            accumulator = operation.calculate(accumulator, carry, a[current].vectorised, b[current].vectorised,
                    size / sizeof(vector double), offset, scalar1, scalar2);
            release_all_blocks();

            MailableResult<double> result = { operation.finish(accumulator) };

            return result.mail;
        }
    }
}

