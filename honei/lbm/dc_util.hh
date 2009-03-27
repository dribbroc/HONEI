/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2009 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the LBM C++ library. LBM is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LBM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef LBM_GUARD_DC_UTIL_HH
#define LBM_GUARD_DC_UTIL_HH 1

namespace honei
{
    namespace lbm
    {
        class VELOC;

        class Reynolds
        {
            public:
                static double value(double veloc, double g_w, double d_x, double d_t, double tau)
                {
                    return veloc * (g_w * d_x) / (d_x * d_x * (2. * tau - 1.)/ (6. * d_t));
                }

                static double adjust_veloc(double reynolds, double veloc, double g_w, double d_x, double d_t, double tau)
                {
                    return (reynolds * (d_x * d_x / d_t) * (2. * tau - 1.))/(6. * g_w * d_x);
                }
        };
    }
}
#endif
