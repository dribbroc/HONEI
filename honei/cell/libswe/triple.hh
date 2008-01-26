/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
 *
 * This file is part of the SWE C++ library. LibSWE is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibSWE is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef CELL_GUARD_LIBSWE_TRIPLE_HH
#define CELL_GUARD_LIBSWE_TRIPLE_HH 1

#include <spu_intrinsics.h>

namespace honei
{
    namespace cell
    {
        namespace intern
        {
            const vector unsigned char extract_patterns_h[2] = {
                { 0x00, 0x01, 0x02, 0x03, 0x0C, 0x0D, 0x0E, 0x0F, 0x18, 0x19, 0x1A, 0x1B, 0xFF, 0xFF, 0xFF, 0xFF },
                { 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A, 0x0B, 0x14, 0x15, 0x16, 0x17 }
            };
            const vector unsigned char extract_patterns_q1[2] = {
                { 0xFF, 0xFF, 0xFF, 0xFF, 0x00, 0x01, 0x02, 0x03, 0x0C, 0x0D, 0x0E, 0x0F, 0x18, 0x19, 0x1A, 0x1B },
                { 0x14, 0x15, 0x16, 0x17, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x0E, 0x0F }
            };
            const vector unsigned char extract_patterns_q2[2] = {
                { 0x18, 0x19, 0x1A, 0x1B, 0xFF, 0xFF, 0xFF, 0xFF, 0x00, 0x01, 0x02, 0x03, 0x0C, 0x0D, 0x0E, 0x0F },
                { 0x00, 0x01, 0x02, 0x03, 0x14, 0x15, 0x16, 0x17, 0x08, 0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x0E, 0x0F }
            };

            const vector unsigned char merge_patterns_result1[2] = {
                { 0x00, 0x01, 0x02, 0x03, 0x10, 0x11, 0x12, 0x13, 0xFF, 0xFF, 0xFF, 0xFF, 0x04, 0x05, 0x06, 0x07 },
                { 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x10, 0x11, 0x12, 0x13, 0x0C, 0x0D, 0x0E, 0x0F }
            };
            const vector unsigned char merge_patterns_result2[2] = {
                { 0x04, 0x05, 0x06, 0x07, 0x14, 0x15, 0x16, 0x17, 0xFF, 0xFF, 0xFF, 0xFF, 0x08, 0x09, 0x0A, 0x0B },
                { 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x18, 0x19, 0x1A, 0x1B, 0x0C, 0x0D, 0x0E, 0x0F }
            };
            const vector unsigned char merge_patterns_result3[2] = {
                { 0x08, 0x09, 0x0A, 0x0B, 0x1C, 0x1D, 0x1E, 0x1F, 0xFF, 0xFF, 0xFF, 0xFF, 0x0C, 0x0D, 0x0E, 0x0F },
                { 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x1C, 0x1D, 0x1E, 0x1F, 0x0C, 0x0D, 0x0E, 0x0F }
            };
        }

        inline vector float extract_h(const vector float & a, const vector float & b, const vector float & c)
        {
            vector float result(spu_shuffle(a, b, intern::extract_patterns_h[0]));

            return spu_shuffle(result, c, intern::extract_patterns_h[1]);
        }

        inline vector float extract_q1(const vector float & a, const vector float & b, const vector float & c)
        {
            vector float result(spu_shuffle(b, c, intern::extract_patterns_q1[0]));

            return spu_shuffle(result, a, intern::extract_patterns_q1[1]);
        }

        inline vector float extract_q2(const vector float & a, const vector float & b, const vector float & c)
        {
            vector float result(spu_shuffle(c, a, intern::extract_patterns_q2[0]));

            return spu_shuffle(result, b, intern::extract_patterns_q2[1]);
        }

        inline vector float merge_result1(const vector float & h, const vector float & q1, const vector float & q2)
        {
            vector float result(spu_shuffle(h, q1, intern::merge_patterns_result1[0]));

            return spu_shuffle(result, q2, intern::merge_patterns_result1[1]);
        }

        inline vector float merge_result2(const vector float & h, const vector float & q1, const vector float & q2)
        {
            vector float result(spu_shuffle(q1, q2, intern::merge_patterns_result2[0]));

            return spu_shuffle(result, h, intern::merge_patterns_result2[1]);
        }

        inline vector float merge_result3(const vector float & h, const vector float & q1, const vector float & q2)
        {
            vector float result(spu_shuffle(q2, h, intern::merge_patterns_result3[0]));

            return spu_shuffle(result, q1, intern::merge_patterns_result3[1]);
        }
    }
}

#endif
