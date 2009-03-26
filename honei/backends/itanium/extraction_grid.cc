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
        template <typename DT1_>
        void extraction_grid_wet(unsigned long begin, unsigned long end,
                DT1_ * distribution_x, DT1_ * distribution_y,
                DT1_ * h, DT1_ * u, DT1_ * v,
                DT1_ * f_0, DT1_ * f_1, DT1_ * f_2,
                DT1_ * f_3, DT1_ * f_4, DT1_ * f_5,
                DT1_ * f_6, DT1_ * f_7, DT1_ * f_8, DT1_ epsilon)
        {
            for (unsigned long index(begin) ; index < end ; ++index)
            {
                h[index] = f_0[index] + f_1[index] + f_2[index] + f_3[index] + f_4[index] +
                    f_5[index] + f_6[index] + f_7[index] + f_8[index];

                u[index] = (distribution_x[0] * f_0[index] +
                        distribution_x[1] * f_1[index] +
                        distribution_x[2] * f_2[index] +
                        distribution_x[3] * f_3[index] +
                        distribution_x[4] * f_4[index] +
                        distribution_x[5] * f_5[index] +
                        distribution_x[6] * f_6[index] +
                        distribution_x[7] * f_7[index] +
                        distribution_x[8] * f_8[index]) / h[index];

                v[index] = (distribution_y[0] * f_0[index] +
                        distribution_y[1] * f_1[index] +
                        distribution_y[2] * f_2[index] +
                        distribution_y[3] * f_3[index] +
                        distribution_y[4] * f_4[index] +
                        distribution_y[5] * f_5[index] +
                        distribution_y[6] * f_6[index] +
                        distribution_y[7] * f_7[index] +
                        distribution_y[8] * f_8[index]) / h[index];
            }
        }
    }
}

template void honei::itanium::extraction_grid_wet<float>(unsigned long, unsigned long,
        float *, float *,
        float *, float *, float *,
        float *, float *, float *,
        float *, float *, float *,
        float *, float *, float *, float);
template void honei::itanium::extraction_grid_wet<double>(unsigned long, unsigned long,
        double *, double *,
        double *, double *, double *,
        double *, double *, double *,
        double *, double *, double *,
        double *, double *, double *, double);
