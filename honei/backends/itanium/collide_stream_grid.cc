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
        void collide_stream_grid_dir_0(unsigned long begin, unsigned long end, DT1_ tau,
                DT1_ * f_temp_0, DT1_ * f_0, DT1_ * f_eq_0)
        {
            for (unsigned long index(begin) ; index < end ; ++index)
            {
                f_temp_0[index] = f_0[index] - (f_0[index] - f_eq_0[index])/tau;
            }
        }

        template <typename DT1_>
        void collide_stream_grid_dir_n(unsigned long end, DT1_ tau,
                unsigned long * dir, unsigned long * dir_index,
                DT1_ * f_temp, DT1_ * f, DT1_ * f_eq)
        {
            for (unsigned long begin(0), half(0) ; begin < end; begin+=2, ++half)
            {
                unsigned long start(dir_index[begin]);
                unsigned long offset(0);
                for (unsigned long index(start) ; index < dir_index[begin + 1] ; ++index, ++offset)
                {
                    f_temp[dir[half] + offset] = f[index] - (f[index] - f_eq[index])/tau;
                }
            }
        }
    }
}

template void honei::itanium::collide_stream_grid_dir_0<float>(unsigned long, unsigned long, float,
        float *, float *, float *);
template void honei::itanium::collide_stream_grid_dir_0<double>(unsigned long, unsigned long, double,
        double *, double *, double *);

template void honei::itanium::collide_stream_grid_dir_n(unsigned long, float,
        unsigned long *, unsigned long *,
        float*, float *, float *);
template void honei::itanium::collide_stream_grid_dir_n(unsigned long, double,
        unsigned long *, unsigned long *,
        double*, double *, double *);
