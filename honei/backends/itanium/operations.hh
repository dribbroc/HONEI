/* vim: set sw=4 sts=4 et foldmethod=syntax : */

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

#ifndef ITANIUM_GUARD_OPERATIONS_HH
#define ITANIUM_GUARD_OPERATIONS_HH 1

namespace honei
{
    namespace itanium
    {
        ///////////// LA
        template <typename DT1_, typename DT2_>
        void scaled_sum(DT1_ * x, const DT2_ * y, DT2_ b, unsigned long size);

        ///////////// LBM
        template <typename DT1_>
        void eq_dist_grid_dir_0(unsigned long begin, unsigned long end,
                DT1_ g, DT1_ e,
                DT1_ * h, DT1_ * u, DT1_ * v,
                DT1_ * f_eq_0);

        template <typename DT1_>
        void eq_dist_grid_dir_odd(unsigned long begin, unsigned long end,
                DT1_ g, DT1_ e,
                DT1_ * h, DT1_ * u, DT1_ * v,
                DT1_ * distribution_x, DT1_ * distribution_y,
                DT1_ * f_eq,
                unsigned long dir);

        template <typename DT1_>
        void eq_dist_grid_dir_even(unsigned long begin, unsigned long end,
                DT1_ g, DT1_ e,
                DT1_ * h, DT1_ * u, DT1_ * v,
                DT1_ * distribution_x, DT1_ * distribution_y,
                DT1_ * f_eq,
                unsigned long dir);

        template <typename DT1_>
        void collide_stream_grid_dir_0(unsigned long begin, unsigned long end, DT1_ tau,
                DT1_ * f_temp_0, DT1_ * f_0, DT1_ * f_eq_0);


        template <typename DT1_>
        void collide_stream_grid_dir_n(unsigned long end, DT1_ tau,
                unsigned long * dir, unsigned long * dir_index,
                DT1_ * f_temp, DT1_ * f, DT1_ * f_eq);

        template <typename DT1_>
        void extraction_grid_wet(unsigned long begin, unsigned long end,
                DT1_ * distribution_x, DT1_ * distribution_y,
                DT1_ * h, DT1_ * u, DT1_ * v,
                DT1_ * f_0, DT1_ * f_1, DT1_ * f_2,
                DT1_ * f_3, DT1_ * f_4, DT1_ * f_5,
                DT1_ * f_6, DT1_ * f_7, DT1_ * f_8, DT1_ epsilon);
    }
}

#endif
