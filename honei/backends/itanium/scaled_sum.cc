/* vim: set sw=4 sts=4 et nofoldenable : */

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


namespace honei
{
    namespace itanium
    {
        template <typename DT1_, typename DT2_>
        void scaled_sum(DT1_ * x, const DT2_ * y, DT2_ b, unsigned long size)
        {
            for (unsigned long index(0) ; index < size ; ++index)
            {
                x[index] += y[index] * b;
            }

        }
    }
}

template void honei::itanium::scaled_sum<float, float>(float *, const float *, float, unsigned long);
template void honei::itanium::scaled_sum<double, double>(double *, const double *, double, unsigned long);
