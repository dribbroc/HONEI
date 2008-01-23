/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#include <honei/cell/cell.hh>
#include <honei/cell/libla/operations.hh>
#include <honei/cell/libutil/debug.hh>
#include <honei/cell/libutil/profiler.hh>

#include <spu_intrinsics.h>

namespace honei
{
    namespace cell
    {
        namespace implementation
        {
            void sqrt_dense_float(vector float * elements, const unsigned size, const float scalar)
            {
                const vector float c0(spu_splats(+2.44992202153062e-01f));
                const vector float c1(spu_splats(+1.24103772746144e+00f));
                const vector float c2(spu_splats(-8.32603540870303e-01f));
                const vector float c3(spu_splats(+4.18723092448099e-01f));
                const vector float c4(spu_splats(+9.90260860770548e-02f));
                const vector float c5(spu_splats(-4.12165065694119e-01f));
                const vector float c6(spu_splats(+4.15153097188655e-01f));
                const vector float c7(spu_splats(-2.57347272326779e-01f));
                const vector float c8(spu_splats(+1.11073254285445e-01f));
                const vector float c9(spu_splats(-3.45276978305398e-02f));
                const vector float c10(spu_splats(+7.74109744341578e-03f));
                const vector float c11(spu_splats(-1.22485258337576e-03f));
                const vector float c12(spu_splats(+1.29996363219277e-04f));
                const vector float c13(spu_splats(-8.31284792897932e-06f));
                const vector float c14(spu_splats(+2.42244634980064e-07f));
                const vector float one(spu_splats(1.0f));
                const vector float two(spu_splats(2.0f));
                const vector unsigned exponent_mask(spu_splats(0x7F800000U));
                const vector unsigned exponent_sign_mask(spu_splats(0x40000000U));
                const vector unsigned exponent_even_mask(spu_splats(0x00800000U));
                const vector unsigned sign_mask(spu_splats(0x80000000U));

                profiler_start();

                for (unsigned i(0) ; i < size ; ++i)
                {
                    vector unsigned exponent(spu_and(
                                exponent_mask,
                                reinterpret_cast<vector unsigned>(elements[i])));
                    vector unsigned exponent_is_negative(spu_cmpeq(
                                exponent_sign_mask,
                                spu_and(exponent, exponent_sign_mask)));
                    vector unsigned exponent_is_even(spu_cmpeq(
                                exponent_even_mask,
                                spu_and(exponent, exponent_even_mask)));

                    // Compute the input's exponent.
                    vector unsigned input_exponent(spu_sel(
                                reinterpret_cast<vector unsigned>(two),
                                reinterpret_cast<vector unsigned>(one),
                                exponent_is_even));

                    // Compute the result's exponent.
                    vector unsigned result_exponent(spu_rlmask(exponent, -23));
                    result_exponent = spu_add(result_exponent, -127);

                    vector unsigned result_exponent_is_negative(spu_cmpeq(
                                spu_and(result_exponent, sign_mask),
                                sign_mask));

                    vector unsigned result_exponent_carry(spu_add(result_exponent, 1U));
                    result_exponent = spu_sel(
                            result_exponent,
                            result_exponent_carry,
                            result_exponent_is_negative);

                    result_exponent = spu_rlmask(result_exponent, -1);

                    vector unsigned result_exponent_odd_negative = spu_add(result_exponent, 126);
                    vector unsigned result_exponent_odd_positive = spu_add(result_exponent, 127);
                    result_exponent = spu_add(result_exponent, 127);
                    vector unsigned result_exponent_odd(spu_sel(
                                result_exponent_odd_negative,
                                result_exponent_odd_positive,
                                exponent_is_negative));
                    result_exponent = spu_sel(
                            result_exponent,
                            result_exponent_odd,
                            spu_nor(
                                spu_splats(0U),
                                spu_xor(exponent_is_even, exponent_is_negative)));

                    result_exponent = reinterpret_cast<vector unsigned>(si_roti(reinterpret_cast<vector signed char>(result_exponent), 23));

                    // Interpolate using a polynomial of grade 14.
                    vector float x1(reinterpret_cast<vector float>(spu_or(input_exponent, spu_and(reinterpret_cast<vector unsigned>(elements[i]),
                                        spu_nor(spu_splats(0U), exponent_mask)))));
                    vector float x2(spu_mul(x1, x1));
                    vector float x4(spu_mul(x2, x2));
                    vector float x8(spu_mul(x4, x4));
                    vector float x12(spu_mul(x4, x8));

                    vector float r12(spu_madd(x1, c13, c12));
                    r12 = spu_madd(x2, c14, r12);
                    r12 = spu_mul(x12, r12);

                    vector float r8(spu_madd(x1, c11, c10));
                    r8 = spu_madd(x2, r8, c8);
                    r8 = spu_madd(x1, c9, r8);
                    r8 = spu_mul(x8, r8);

                    vector float r4(spu_madd(x1, c7, c6));
                    r4 = spu_madd(x2, r4, c4);
                    r4 = spu_madd(x1, c5, r4);
                    r4 = spu_mul(x4, r4);

                    vector float r2(spu_madd(x1, c3, c2));
                    r2 = spu_mul(x2, r2);

                    vector float r1(spu_madd(x1, c1, c0));

                    vector float y0(spu_add(spu_add(spu_add(spu_add(r12, r8), r4), r2), r1));

                    // y0 is our first guess. Enhance precision by a Newton step.
                    // We're searching for y that fullfils: 0 = f(y) = x - y^2,
                    // where x is the original value whose square-root we're computing.
                    // y1 = y0 + f(y0) / f'(y0) = y0 + (x - y0^2) / (2 * y0)
                    vector float inv_v(spu_mul(two, y0));
                    vector float inv_i(spu_re(inv_v));
                    vector float inv_t(spu_nmsub(inv_v, inv_i, one));
                    inv_v = spu_madd(inv_t, inv_i, inv_i);

                    vector float y1(spu_madd(inv_v, spu_sub(x1, spu_mul(y0, y0)), y0));

                    // Reinsert the exponent.
                    vector unsigned value(spu_and(
                                reinterpret_cast<vector unsigned>(y1),
                                spu_nor(spu_splats(0U), exponent_mask)));
                    y1 = reinterpret_cast<vector float>(spu_or(
                                value,
                                spu_and(exponent_mask, result_exponent)));

                    elements[i] = y1;
                }

                profiler_stop();
            }
        }

        namespace operations
        {
            Operation<1, float, rtm_dma> sqrt_dense_float = {
                &implementation::sqrt_dense_float,
            };
        }
    }
}

