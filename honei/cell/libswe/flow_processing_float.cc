/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Sven Mallach <sven.mallach@honei.org>
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
#include <honei/cell/libswe/operations.hh>
#include <honei/cell/libswe/triple.hh>

#include <spu_intrinsics.h>

namespace honei
{
    namespace cell
    {
        namespace implementation
        {
            void flow_processing_x_float(vector float * elements, const unsigned size, const float)
            {
                const vector float one(spu_splats(1.0f));
                const vector float gravity_by_two(spu_splats(9.81f / 2.0f));
                const vector float eps(spu_splats(1.19209e-07f));
                const vector float eps_inv(spu_splats(8.3886e+06f));

                for (unsigned i(0) ; i + 2 < size ; i += 3)
                {
                    vector float h(extract_h(elements[i], elements[i + 1], elements[i + 2]));
                    vector float q1(extract_q1(elements[i], elements[i + 1], elements[i + 2]));
                    vector float q2(extract_q2(elements[i], elements[i + 1], elements[i + 2]));

                    vector float h_inv(spu_re(h));
                    vector float h_inv_t(spu_nmsub(h, h_inv, one));
                    h_inv = spu_madd(h_inv_t, h_inv, h_inv);

                    vector unsigned h_smaller_than_eps(spu_cmpabsgt(eps, h));

                    vector float & result_h(q1);
                    vector float result_q1(spu_madd(spu_mul(q1, h_inv), q1, spu_mul(spu_mul(h, h), gravity_by_two)));
                    vector float result_q2(spu_mul(spu_mul(q1, h_inv), q2));

                    vector float result_q1_eps(spu_mul(spu_mul(q1, eps_inv), q1));
                    vector float result_q2_eps(spu_mul(spu_mul(q1, eps_inv), q2));

                    result_q1 = spu_sel(result_q1, result_q1_eps, h_smaller_than_eps);
                    result_q2 = spu_sel(result_q2, result_q2_eps, h_smaller_than_eps);

                    elements[i] = merge_result1(result_h, result_q1, result_q2);
                    elements[i + 1] = merge_result2(result_h, result_q1, result_q2);
                    elements[i + 2] = merge_result3(result_h, result_q1, result_q2);
                }
            }

            void flow_processing_y_float(vector float * elements, const unsigned size, const float scalar)
            {
                const vector float one(spu_splats(1.0f));
                const vector float gravity_by_two(spu_splats(9.81f / 2.0f));
                const vector float eps(spu_splats(1.19209e-07f));
                const vector float eps_inv(spu_splats(8.3886e+06f));

                for (unsigned i(0) ; i + 2 < size ; i += 3)
                {
                    vector float h(extract_h(elements[i], elements[i + 1], elements[i + 2]));
                    vector float q1(extract_q1(elements[i], elements[i + 1], elements[i + 2]));
                    vector float q2(extract_q2(elements[i], elements[i + 1], elements[i + 2]));

                    vector float h_inv(spu_re(h));
                    vector float h_inv_t(spu_nmsub(h, h_inv, one));
                    h_inv = spu_madd(h_inv_t, h_inv, h_inv);

                    vector unsigned h_smaller_than_eps(spu_cmpabsgt(eps, h));

                    vector float & result_h(q2);
                    vector float result_q1(spu_mul(spu_mul(q1, h_inv), q2));
                    vector float result_q2(spu_madd(spu_mul(q2, h_inv), q2, spu_mul(spu_mul(h, h), gravity_by_two)));

                    vector float result_q1_eps(spu_mul(spu_mul(q1, eps_inv), q2));
                    vector float result_q2_eps(spu_mul(spu_mul(q2, eps_inv), q2));

                    result_q1 = spu_sel(result_q1, result_q1_eps, h_smaller_than_eps);
                    result_q2 = spu_sel(result_q2, result_q2_eps, h_smaller_than_eps);

                    elements[i] = merge_result1(result_h, result_q1, result_q2);
                    elements[i + 1] = merge_result2(result_h, result_q1, result_q2);
                    elements[i + 2] = merge_result3(result_h, result_q1, result_q2);
                }
            }
        }

        namespace operations
        {
            Operation<1, float, rtm_dma> flow_processing_x_float = {
                &implementation::flow_processing_x_float,
            };

            Operation<1, float, rtm_dma> flow_processing_y_float = {
                &implementation::flow_processing_y_float,
            };
        }
    }
}

