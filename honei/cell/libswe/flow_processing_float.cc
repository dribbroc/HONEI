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
#include <honei/cell/libutil/transfer.hh>

#include <spu_intrinsics.h>

namespace honei
{
    namespace cell
    {
        namespace implementation
        {
            void flow_processing_x_float(vector float * elements, const unsigned size, const float scalar)
            {
                const vector float one(spu_splats(1.0f));
                const vector float gravity_by_two(spu_splats(9.81f / 2.0f));
                const vector unsigned char extract_patterns[2] = {
                    { 0x00, 0x01, 0x02, 0x03, 0x0C, 0x0D, 0x0E, 0x0F, 0x18, 0x19, 0x1A, 0x1B, 0x1F, 0x1F, 0x1F, 0x1F },
                    { 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A, 0x0B, 0x14, 0x15, 0x16, 0x17 }
                };
                const vector unsigned char insert_patterns[2] = {
                    { 0x00, 0x01, 0x02, 0x03, 0x10, 0x11, 0x12, 0x13, 0xFF, 0xFF, 0xFF, 0xFF, 0x04, 0x05, 0x06, 0x07 },
                    { 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x10, 0x11, 0x12, 0x13, 0x0C, 0x0D, 0x0E, 0x0F }
                };

                for (unsigned i(0) ; i + 2 < size ; i += 3)
                {
                    vector float h(spu_shuffle(elements[i], elements[i + 1], extract_patterns[0]));
                    h = spu_shuffle(h, elements[i + 2], extract_patterns[1]);

                    vector float q1(spu_shuffle(elements[i + 1], elements[i + 2], extract_patterns[0]));
                    q1 = spu_shuffle(q1, elements[i], extract_patterns[1]);

                    vector float q2(spu_shuffle(elements[i + 2], elements[i], extract_patterns[0]));
                    q2 = spu_shuffle(q2, elements[i + 1], extract_patterns[1]);

                    vector float h_inv(spu_re(h));
                    vector float h_inv_t(spu_nmsub(h, h_inv, one));
                    h_inv = spu_madd(h_inv_t, h_inv, h_inv);

                    vector float & result_h(q1);
                    vector float result_q1(spu_madd(spu_mul(q1, h_inv), q1, spu_mul(spu_mul(h, h), gravity_by_two)));
                    vector float result_q2(spu_mul(spu_mul(q1, h_inv), q2));

                    elements[i] = spu_shuffle(result_h, result_q1, insert_patterns[0]);
                    elements[i] = spu_shuffle(elements[i], result_q2, insert_patterns[1]);

                    elements[i + 1] = spu_shuffle(result_q1, result_q2, insert_patterns[0]);
                    elements[i + 1] = spu_shuffle(elements[i + 1], result_h, insert_patterns[1]);

                    elements[i + 2] = spu_shuffle(result_q2, result_h, insert_patterns[0]);
                    elements[i + 2] = spu_shuffle(elements[i + 2], result_q1, insert_patterns[1]);
                }
            }

            void flow_processing_y_float(vector float * elements, const unsigned size, const float scalar)
            {
                const vector float one(spu_splats(1.0f));
                const vector float gravity_by_two(spu_splats(9.81f / 2.0f));
                const vector unsigned char extract_patterns[2] = {
                    { 0x00, 0x01, 0x02, 0x03, 0x0C, 0x0D, 0x0E, 0x0F, 0x18, 0x19, 0x1A, 0x1B, 0x1F, 0x1F, 0x1F, 0x1F },
                    { 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A, 0x0B, 0x14, 0x15, 0x16, 0x17 }
                };
                const vector unsigned char insert_patterns[2] = {
                    { 0x00, 0x01, 0x02, 0x03, 0x10, 0x11, 0x12, 0x13, 0xFF, 0xFF, 0xFF, 0xFF, 0x04, 0x05, 0x06, 0x07 },
                    { 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x10, 0x11, 0x12, 0x13, 0x0C, 0x0D, 0x0E, 0x0F }
                };

                for (unsigned i(0) ; i + 2 < size ; i += 3)
                {
                    vector float h(spu_shuffle(elements[i], elements[i + 1], extract_patterns[0]));
                    h = spu_shuffle(h, elements[i + 2], extract_patterns[1]);

                    vector float q1(spu_shuffle(elements[i + 1], elements[i + 2], extract_patterns[0]));
                    q1 = spu_shuffle(q1, elements[i], extract_patterns[1]);

                    vector float q2(spu_shuffle(elements[i + 2], elements[i], extract_patterns[0]));
                    q2 = spu_shuffle(q2, elements[i + 1], extract_patterns[1]);

                    vector float h_inv(spu_re(h));
                    vector float h_inv_t(spu_nmsub(h, h_inv, one));
                    h_inv = spu_madd(h_inv_t, h_inv, h_inv);

                    vector float & result_h(q2);
                    vector float result_q1(spu_mul(spu_mul(q1, h_inv), q2));
                    vector float result_q2(spu_madd(spu_mul(q2, h_inv), q2, spu_mul(spu_mul(h, h), gravity_by_two)));

                    elements[i] = spu_shuffle(result_h, result_q1, insert_patterns[0]);
                    elements[i] = spu_shuffle(elements[i], result_q2, insert_patterns[1]);

                    elements[i + 1] = spu_shuffle(result_q1, result_q2, insert_patterns[0]);
                    elements[i + 1] = spu_shuffle(elements[i + 1], result_h, insert_patterns[1]);

                    elements[i + 2] = spu_shuffle(result_q2, result_h, insert_patterns[0]);
                    elements[i + 2] = spu_shuffle(elements[i + 2], result_q1, insert_patterns[1]);
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

