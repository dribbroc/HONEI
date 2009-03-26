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
        void eq_dist_grid_dir_0(unsigned long begin, unsigned long end,
                DT1_ g, DT1_ e,
                DT1_ * h, DT1_ * u, DT1_ * v,
                DT1_ * f_eq_0)
        {

            DT1_ e2(e);
            DT1_ e23(DT1_(3.) * e2);
            DT1_ e26(DT1_(6.) * e2);
            for (unsigned long index(begin) ; index < end ; ++index)
            {
                DT1_ u2(u[index] * u[index]);
                DT1_ v2(v[index] * v[index]);
                DT1_ h2(h[index] * h[index]);

                DT1_ t1, t2;

                t1 = (DT1_(5.) * g * h2) / e26;
                t2 = (DT1_(2.) * h[index]) / e23 * (u2 + v2);
                f_eq_0[index] = h[index] - t1 - t2;
            }
        }

        template <typename DT1_>
        void eq_dist_grid_dir_odd(unsigned long begin, unsigned long end,
                DT1_ g, DT1_ e,
                DT1_ * h, DT1_ * u, DT1_ * v,
                DT1_ * distribution_x, DT1_ * distribution_y,
                DT1_ * f_eq,
                unsigned long dir)
        {

            DT1_ e2(e);
            DT1_ e42(DT1_(2.) * e2 * e2);
            DT1_ e23(DT1_(3.) * e2);
            DT1_ e26(DT1_(6.) * e2);
            for (unsigned long index(begin) ; index < end ; ++index)
            {
                DT1_ u2(u[index] * u[index]);
                DT1_ v2(v[index] * v[index]);
                DT1_ h2(h[index] * h[index]);

                DT1_ dxu, dyv;
                DT1_ t1, t2, t3, t4;
                dxu = distribution_x[dir] * u[index];
                dyv = distribution_y[dir] * v[index];
                t1 = (g * h2) / e26;
                t2 = (h[index] / e23) * (dxu + dyv);
                t3 = (h[index] / e42) * (dxu * dxu + DT1_(2.) * dxu * dyv + dyv * dyv);
                t4 = (h[index] / e26) * (u2 + v2);
                f_eq[index] = t1 + t2 + t3 - t4;
            }
        }

        template <typename DT1_>
        void eq_dist_grid_dir_even(unsigned long begin, unsigned long end,
                DT1_ g, DT1_ e,
                DT1_ * h, DT1_ * u, DT1_ * v,
                DT1_ * distribution_x, DT1_ * distribution_y,
                DT1_ * f_eq,
                unsigned long dir)
        {
            DT1_ e2(e);
            DT1_ e48(DT1_(8.) * e2 * e2);
            DT1_ e212(DT1_(12.) * e2);
            DT1_ e224(DT1_(24.) * e2);
            for (unsigned long index(begin) ; index < end ; ++index)
            {
                DT1_ u2(u[index] * u[index]);
                DT1_ v2(v[index] * v[index]);
                DT1_ h2(h[index] * h[index]);

                DT1_ dxu, dyv;
                DT1_ t1, t2, t3, t4;
                dxu = distribution_x[dir] * u[index];
                dyv = distribution_y[dir] * v[index];
                t1 = (g * h2) / e224;
                t2 = (h[index] / e212) * (dxu + dyv);
                t3 = (h[index] / e48) * (dxu * dxu + DT1_(2.) * dxu * dyv + dyv * dyv);
                t4 = (h[index] / e224) * (u2 + v2);
                f_eq[index] =  t1 + t2 + t3 - t4;
            }
        }
    }
}
template void honei::itanium::eq_dist_grid_dir_0<float>(unsigned long, unsigned long,
        float, float,
        float*, float*, float *,
        float*);
template void honei::itanium::eq_dist_grid_dir_0<double>(unsigned long, unsigned long,
        double, double,
        double*, double*, double *,
        double*);

template void honei::itanium::eq_dist_grid_dir_odd<float>(unsigned long, unsigned long,
        float, float,
        float *, float *,float *,
        float *, float *,
        float *,
        unsigned long);
template void honei::itanium::eq_dist_grid_dir_odd<double>(unsigned long, unsigned long,
        double, double,
        double *, double *,double *,
        double *, double *,
        double *,
        unsigned long);

template void honei::itanium::eq_dist_grid_dir_even<float>(unsigned long, unsigned long,
        float, float,
        float *, float *,float *,
        float *, float *,
        float *,
        unsigned long);
template void honei::itanium::eq_dist_grid_dir_even<double>(unsigned long, unsigned long,
        double, double,
        double *, double *,double *,
        double *, double *,
        double *,
        unsigned long);
