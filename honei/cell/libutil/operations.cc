/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#include <honei/cell/types.hh>
#include <honei/cell/libutil/operations.hh>

namespace honei
{
    namespace cell
    {
        vector float biggest_float()
        {
            static const vector unsigned result = { 0x7F7FFFFF, 0x7F7FFFFF, 0x7F7FFFFF, 0x7F7FFFFF };

            return reinterpret_cast<vector float>(result);
        }

        float max_float(const vector float & value)
        {
            Subscriptable<float> temp = { value };
            float result(temp.array[0]);

            result = (result > temp.array[1]) ? result : temp.array[1];
            result = (result > temp.array[2]) ? result : temp.array[2];
            result = (result > temp.array[3]) ? result : temp.array[3];

            return result;
        }

        float min_float(const vector float & value)
        {
            Subscriptable<float> temp = { value };
            float result(temp.array[0]);

            result = (result < temp.array[1]) ? result : temp.array[1];
            result = (result < temp.array[2]) ? result : temp.array[2];
            result = (result < temp.array[3]) ? result : temp.array[3];

            return result;
        }

        vector float smallest_float()
        {
            static const vector unsigned result = { 0x00800000, 0x00800000, 0x00800000, 0x00800000 };

            return reinterpret_cast<vector float>(result);
        }

        float sum_float(const vector float & value)
        {
            Subscriptable<float> temp = { value };

            return temp.array[0] + temp.array[1] + temp.array[2] + temp.array[3];
        }

        vector float zero_float()
        {
            static const vector float result = { 0.0f, 0.0f, 0.0f, 0.0f };

            return result;
        }

    }
}
