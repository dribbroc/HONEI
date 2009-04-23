/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2009 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
#include <honei/backends/cell/spe/liblbm/operations.hh>
#include <honei/backends/cell/spe/libutil/allocator.hh>
#include <honei/backends/cell/spe/libutil/transfer.hh>

#include <spu_intrinsics.h>

namespace honei
{
    namespace cell
    {
        namespace implementation
        {
            void collide_stream_grid_float(vector float * f_temp, const vector float * f, const vector float * f_eq,
                    const unsigned size, vector float & f_carry, const unsigned f_offset, vector float & f_eq_carry,
                    const unsigned f_eq_offset, const float tau, const float offsetf, const float timesf)
            {
                vector float temp, vd, id, td, tauv, one;
                tauv = spu_splats(tau);
                one = spu_splats(float(1));
                for (unsigned i(0) ; i < size ; ++i)
                {
                    extract(f_carry, f[i], f_offset);
                    extract(f_eq_carry, f_eq[i], f_eq_offset);

                    temp = spu_sub(f_carry, f_eq_carry);
                    vd = tauv;
                    id = spu_re(vd);
                    td = spu_nmsub(vd, id, one);
                    vd = spu_madd(td, id, id);
                    temp = spu_mul(temp, vd);
                    f_temp[i] = spu_sub(f_carry, temp);

                    f_carry = f[i];
                    f_eq_carry = f_eq[i];
                }
            }
            void collide_stream_grid_times_float(vector float * f_temp, const vector float * f, const vector float * f_eq,
                    const unsigned size, vector float & f_carry, const unsigned f_offset, vector float & f_eq_carry,
                    const unsigned f_eq_offset, const float tau, const float offsetf, const float timesf)
            {
                vector float temp, talt, vd, id, td, tauv, one;
                vector unsigned mod, zero, four;
                tauv = spu_splats(tau);
                one = spu_splats(float(1));
                four = spu_splats(unsigned(4));
                zero = spu_splats(unsigned(0));
                vector unsigned times = spu_splats((unsigned)timesf);
                vector unsigned offset = spu_splats((unsigned)offsetf);
                vector unsigned real_index;
                vector unsigned zero123;
                zero123 = spu_insert((unsigned)0, zero123, 0);
                zero123 = spu_insert((unsigned)1, zero123, 1);
                zero123 = spu_insert((unsigned)2, zero123, 2);
                zero123 = spu_insert((unsigned)3, zero123, 3);
                vector unsigned iv = zero123;
                vector unsigned diff;

                for (unsigned i(0) ; i < size ; ++i)
                {
                    real_index = spu_add(offset, iv);
                    // real_index > times
                    mod = spu_cmpgt(real_index, times);
                    diff = spu_sub(real_index, times);
                    real_index = spu_sel(real_index, diff, mod);

                    // times % real_index
                    mod = spu_sub(times, real_index);

                    //mod = 1 - mod
                    mod = spu_cmpeq(mod, zero);

                    talt = f_temp[i];
                    extract(f_carry, f[i], f_offset);
                    extract(f_eq_carry, f_eq[i], f_eq_offset);

                    temp = spu_sub(f_carry, f_eq_carry);
                    vd = tauv;
                    id = spu_re(vd);
                    td = spu_nmsub(vd, id, one);
                    vd = spu_madd(td, id, id);
                    temp = spu_mul(temp, vd);
                    temp = spu_sub(f_carry, temp);

                    f_carry = f[i];
                    f_eq_carry = f_eq[i];

                    // mod:1 restore | 0 do nothing
                    f_temp[i] = spu_sel(temp, talt, mod);

                    iv = spu_add(iv, four);
                }
            }
        }

        namespace operations
        {
            Operation<3, float, rtm_dma> collide_stream_grid_float = {
                &implementation::collide_stream_grid_float
            };
            Operation<3, float, rtm_dma> collide_stream_grid_times_float = {
                &implementation::collide_stream_grid_times_float
            };
        }
    }
}
