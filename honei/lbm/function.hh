/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2012 Markus Geveler <apryde@gmx.de>
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

#pragma once
#ifndef LBM_GUARD_FUNCTION_HH
#define LBM_GUARD_FUNCTION_HH 1

#include<cmath>

using namespace std;

namespace honei
{

    ///2D Gaussian
    template<int angle_ = 0>      ///rotation angle in grad
    struct Gaussian2D
    {
        template<typename DT_>
        static DT_ value(DT_ x, DT_ y, DT_ amplitude = 1., DT_ x_0 = 0., DT_ y_0 = 0., DT_ x_spread = 1., DT_ y_spread = 1.)
        {
            double pi = 4 * atan(1);
            DT_ angle((DT_)angle_ * DT_(2) * pi / 360);
            DT_ a(pow(cos(angle), DT_(2.))/(DT_(2) * x_spread * x_spread)  + pow(sin(angle), 2.)/(DT_(2) * y_spread * y_spread));
            DT_ b(-sin(DT_(2) * angle)/(DT_(4) * x_spread * x_spread)  + sin(DT_(2) * angle)/(DT_(4) * y_spread * y_spread));
            DT_ c(pow(sin(angle), DT_(2.))/(DT_(2) * x_spread * x_spread)  + pow(cos(angle), 2.)/(DT_(2) * y_spread * y_spread));

            return amplitude * exp( -( a * (x - x_0) * (x - x_0) + DT_(2) * b * (x - x_0) * (y - y_0) + c * (y - y_0) * (y - y_0) ) );
        }
    };
}

#endif
