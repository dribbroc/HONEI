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

#include <honei/util/attributes.hh>

#include <xmmintrin.h>
#include <emmintrin.h>

namespace honei
{
    namespace sse
    {
        void up_vel_dir_grid(unsigned long begin, unsigned long end, unsigned long * limits, unsigned long * types,
                float * f_temp_1, float * f_temp_2,  float * f_temp_3,
                float * f_temp_4, float * f_temp_5,  float * f_temp_6,
                float * f_temp_7, float * f_temp_8)
        {
            for (unsigned long index(begin) ; index < end ; ++index)
            {
                if((types[index] & 1<<0) == 1<<0)
                    for (unsigned long i(limits[index]) ; i != limits[index + 1] ; ++i)
                    {
                        f_temp_5[i] = f_temp_1[i];
                    }
                if((types[index] & 1<<1) == 1<<1)
                    for (unsigned long i(limits[index]) ; i != limits[index + 1] ; ++i)
                    {
                        f_temp_6[i] = f_temp_2[i];
                    }
                if((types[index] & 1<<2) == 1<<2)
                    for (unsigned long i(limits[index]) ; i != limits[index + 1] ; ++i)
                    {
                        f_temp_7[i] = f_temp_3[i];
                    }
                if((types[index] & 1<<3) == 1<<3)
                    for (unsigned long i(limits[index]) ; i != limits[index + 1] ; ++i)
                    {
                        f_temp_8[i] = f_temp_4[i];
                    }
                if((types[index] & 1<<4) == 1<<4)
                    for (unsigned long i(limits[index]) ; i != limits[index + 1] ; ++i)
                    {
                        f_temp_1[i] = f_temp_5[i];
                    }
                if((types[index] & 1<<5) == 1<<5)
                    for (unsigned long i(limits[index]) ; i != limits[index + 1] ; ++i)
                    {
                        f_temp_2[i] = f_temp_6[i];
                    }
                if((types[index] & 1<<6) == 1<<6)
                    for (unsigned long i(limits[index]) ; i != limits[index + 1] ; ++i)
                    {
                        f_temp_3[i] = f_temp_7[i];
                    }
                if((types[index] & 1<<7) == 1<<7)
                    for (unsigned long i(limits[index]) ; i != limits[index + 1] ; ++i)
                    {
                        f_temp_4[i] = f_temp_8[i];
                    }

                // Corners
                /*if((types[index] & 1<<2) == 1<<2 && (types[index] & 1<<4) == 1<<4)
                    for (unsigned long i(limits[index]) ; i != limits[index + 1] ; ++i)
                    {
                        f_temp_2[i] = f_temp_8[i];
                        f_temp_6[i] = f_temp_8[i];
                    }
                if((types[index] & 1<<4) == 1<<4 && (types[index] & 1<<6) == 1<<6)
                    for (unsigned long i(limits[index]) ; i != limits[index + 1] ; ++i)
                    {
                        f_temp_4[i] = f_temp_2[i];
                        f_temp_8[i] = f_temp_2[i];
                    }
                if((types[index] & 1<<0) == 1<<0 && (types[index] & 1<<6) == 1<<6)
                    for (unsigned long i(limits[index]) ; i != limits[index + 1] ; ++i)
                    {
                        f_temp_2[i] = f_temp_4[i];
                        f_temp_6[i] = f_temp_4[i];
                    }
                if((types[index] & 1<<0) == 1<<0 && (types[index] & 1<<2) == 1<<2)
                    for (unsigned long i(limits[index]) ; i != limits[index + 1] ; ++i)
                    {
                        f_temp_4[i] = f_temp_6[i];
                        f_temp_8[i] = f_temp_6[i];
                    }*/
                if(((types)[index] & 1<<2) == 1<<2 &&
                        ((types)[index] & 1<<4) == 1<<4 &&
                        ((types)[index] & 1<<1) == 1<<1 &&
                        ((types)[index] & 1<<5) == 1<<5)
                {
                    (f_temp_2)[index] = (f_temp_8)[index];
                    (f_temp_6)[index] = (f_temp_8)[index];
                }
                if(((types)[index] & 1<<4) == 1<<4 &&
                        ((types)[index] & 1<<6) == 1<<6 &&
                        ((types)[index] & 1<<7) == 1<<7 &&
                        ((types)[index] & 1<<3) == 1<<3)
                {
                    (f_temp_4)[index] = (f_temp_2)[index];
                    (f_temp_8)[index] = (f_temp_2)[index];
                }
                if(((types)[index] & 1<<0) == 1<<0 &&
                        ((types)[index] & 1<<6) == 1<<6 &&
                        ((types)[index] & 1<<1) == 1<<1 &&
                        ((types)[index] & 1<<5) == 1<<5)
                {
                    (f_temp_2)[index] = (f_temp_4)[index];
                    (f_temp_6)[index] = (f_temp_4)[index];
                }
                if(((types)[index] & 1<<0) == 1<<0 &&
                        ((types)[index] & 1<<2) == 1<<2 &&
                        ((types)[index] & 1<<7) == 1<<7 &&
                        ((types)[index] & 1<<3) == 1<<3)
                {
                    (f_temp_4)[index] = (f_temp_6)[index];
                    (f_temp_8)[index] = (f_temp_6)[index];
                }

            }
        }

        void up_vel_dir_grid(unsigned long begin, unsigned long end, unsigned long * limits, unsigned long * types,
                double * f_temp_1, double * f_temp_2,  double * f_temp_3,
                double * f_temp_4, double * f_temp_5,  double * f_temp_6,
                double * f_temp_7, double * f_temp_8)
        {
            for (unsigned long index(begin) ; index < end ; ++index)
            {
                if((types[index] & 1<<0) == 1<<0)
                    for (unsigned long i(limits[index]) ; i != limits[index + 1] ; ++i)
                    {
                        f_temp_5[i] = f_temp_1[i];
                    }
                if((types[index] & 1<<1) == 1<<1)
                    for (unsigned long i(limits[index]) ; i != limits[index + 1] ; ++i)
                    {
                        f_temp_6[i] = f_temp_2[i];
                    }
                if((types[index] & 1<<2) == 1<<2)
                    for (unsigned long i(limits[index]) ; i != limits[index + 1] ; ++i)
                    {
                        f_temp_7[i] = f_temp_3[i];
                    }
                if((types[index] & 1<<3) == 1<<3)
                    for (unsigned long i(limits[index]) ; i != limits[index + 1] ; ++i)
                    {
                        f_temp_8[i] = f_temp_4[i];
                    }
                if((types[index] & 1<<4) == 1<<4)
                    for (unsigned long i(limits[index]) ; i != limits[index + 1] ; ++i)
                    {
                        f_temp_1[i] = f_temp_5[i];
                    }
                if((types[index] & 1<<5) == 1<<5)
                    for (unsigned long i(limits[index]) ; i != limits[index + 1] ; ++i)
                    {
                        f_temp_2[i] = f_temp_6[i];
                    }
                if((types[index] & 1<<6) == 1<<6)
                    for (unsigned long i(limits[index]) ; i != limits[index + 1] ; ++i)
                    {
                        f_temp_3[i] = f_temp_7[i];
                    }
                if((types[index] & 1<<7) == 1<<7)
                    for (unsigned long i(limits[index]) ; i != limits[index + 1] ; ++i)
                    {
                        f_temp_4[i] = f_temp_8[i];
                    }

                // Corners
                if((types[index] & 1<<2) == 1<<2 && (types[index] & 1<<4) == 1<<4)
                    for (unsigned long i(limits[index]) ; i != limits[index + 1] ; ++i)
                    {
                        f_temp_2[i] = f_temp_8[i];
                        f_temp_6[i] = f_temp_8[i];
                    }
                if((types[index] & 1<<4) == 1<<4 && (types[index] & 1<<6) == 1<<6)
                    for (unsigned long i(limits[index]) ; i != limits[index + 1] ; ++i)
                    {
                        f_temp_4[i] = f_temp_2[i];
                        f_temp_8[i] = f_temp_2[i];
                    }
                if((types[index] & 1<<0) == 1<<0 && (types[index] & 1<<6) == 1<<6)
                    for (unsigned long i(limits[index]) ; i != limits[index + 1] ; ++i)
                    {
                        f_temp_2[i] = f_temp_4[i];
                        f_temp_6[i] = f_temp_4[i];
                    }
                if((types[index] & 1<<0) == 1<<0 && (types[index] & 1<<2) == 1<<2)
                    for (unsigned long i(limits[index]) ; i != limits[index + 1] ; ++i)
                    {
                        f_temp_4[i] = f_temp_6[i];
                        f_temp_8[i] = f_temp_6[i];
                    }
            }
        }
    }
}
