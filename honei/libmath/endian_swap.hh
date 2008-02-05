/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the Math C++ library. LibMath is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibMath is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */


#ifndef LIBMATH_GUARD_ENDIAN_SWAP_HH
#define LIBMATH_GUARD_ENDIAN_SWAP_HH 1

#include <honei/libutil/stringify.hh>
#include <iostream>

using namespace honei;
using namespace tests;
using namespace std;
namespace honei
{
    float FloatSwap( float f )
    {
        union
        {
            float f;
            unsigned char b[4];
        } dat1, dat2;

        dat1.f = f;
        dat2.b[0] = dat1.b[3];
        dat2.b[1] = dat1.b[2];
        dat2.b[2] = dat1.b[1];
        dat2.b[3] = dat1.b[0];
        return dat2.f;
    }

     double DoubleSwap( double d )
    {
        union
        {
            double d;
            unsigned char b[8];
        } dat1, dat2;

        dat1.d = d;
        dat2.b[0] = dat1.b[7];
        dat2.b[1] = dat1.b[6];
        dat2.b[2] = dat1.b[5];
        dat2.b[3] = dat1.b[4];
        dat2.b[4] = dat1.b[3];
        dat2.b[5] = dat1.b[2];
        dat2.b[6] = dat1.b[1];
        dat2.b[7] = dat1.b[0];

        return dat2.d;
    } 
}

#endif
