/* vim: set sw=4 sts=4 et nofoldenable filetype=cpp : */

/*
* Copyright (c) 2013 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
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

#ifndef __CPU__
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif

__kernel void update_velocity_directions_grid_f(const unsigned long start, const unsigned long end,
                                      __global unsigned long * types, __global float * f_temp_1, __global float * f_temp_2,
                                      __global float * f_temp_3, __global float * f_temp_4, __global float * f_temp_5,
                                      __global float * f_temp_6, __global float * f_temp_7, __global float * f_temp_8)
{
uint idx = get_global_id(0);
idx += start;
if (idx < end)
{
unsigned long i = idx;
if((types[i] & 1<<0) == 1<<0)
f_temp_5[i] = f_temp_1[i];
if((types[i] & 1<<1) == 1<<1)
f_temp_6[i] = f_temp_2[i];
if((types[i] & 1<<2) == 1<<2)
f_temp_7[i] = f_temp_3[i];
if((types[i] & 1<<3) == 1<<3)
f_temp_8[i] = f_temp_4[i];
if((types[i] & 1<<4) == 1<<4)
f_temp_1[i] = f_temp_5[i];
if((types[i] & 1<<5) == 1<<5)
f_temp_2[i] = f_temp_6[i];
if((types[i] & 1<<6) == 1<<6)
f_temp_3[i] = f_temp_7[i];
if((types[i] & 1<<7) == 1<<7)
f_temp_4[i] = f_temp_8[i];

// Corners
if((types[i] & 1<<2) == 1<<2 && (types[i] & 1<<4) == 1<<4)
{
f_temp_2[i] = f_temp_8[i];
f_temp_6[i] = f_temp_8[i];
}
if((types[i] & 1<<4) == 1<<4 && (types[i] & 1<<6) == 1<<6)
{
f_temp_4[i] = f_temp_2[i];
f_temp_8[i] = f_temp_2[i];
}
if((types[i] & 1<<0) == 1<<0 && (types[i] & 1<<6) == 1<<6)
{
f_temp_2[i] = f_temp_4[i];
f_temp_6[i] = f_temp_4[i];
}
if((types[i] & 1<<0) == 1<<0 && (types[i] & 1<<2) == 1<<2)
{
f_temp_4[i] = f_temp_6[i];
f_temp_8[i] = f_temp_6[i];
}
}
}

__kernel void update_velocity_directions_grid_d(const unsigned long start, const unsigned long end,
                                      __global unsigned long * types, __global double * f_temp_1, __global double * f_temp_2,
                                      __global double * f_temp_3, __global double * f_temp_4, __global double * f_temp_5,
                                      __global double * f_temp_6, __global double * f_temp_7, __global double * f_temp_8)
{
uint idx = get_global_id(0);
idx += start;
if (idx < end)
{
unsigned long i = idx;
if((types[i] & 1<<0) == 1<<0)
f_temp_5[i] = f_temp_1[i];
if((types[i] & 1<<1) == 1<<1)
f_temp_6[i] = f_temp_2[i];
if((types[i] & 1<<2) == 1<<2)
f_temp_7[i] = f_temp_3[i];
if((types[i] & 1<<3) == 1<<3)
f_temp_8[i] = f_temp_4[i];
if((types[i] & 1<<4) == 1<<4)
f_temp_1[i] = f_temp_5[i];
if((types[i] & 1<<5) == 1<<5)
f_temp_2[i] = f_temp_6[i];
if((types[i] & 1<<6) == 1<<6)
f_temp_3[i] = f_temp_7[i];
if((types[i] & 1<<7) == 1<<7)
f_temp_4[i] = f_temp_8[i];

// Corners
if((types[i] & 1<<2) == 1<<2 && (types[i] & 1<<4) == 1<<4)
{
f_temp_2[i] = f_temp_8[i];
f_temp_6[i] = f_temp_8[i];
}
if((types[i] & 1<<4) == 1<<4 && (types[i] & 1<<6) == 1<<6)
{
f_temp_4[i] = f_temp_2[i];
f_temp_8[i] = f_temp_2[i];
}
if((types[i] & 1<<0) == 1<<0 && (types[i] & 1<<6) == 1<<6)
{
f_temp_2[i] = f_temp_4[i];
f_temp_6[i] = f_temp_4[i];
}
if((types[i] & 1<<0) == 1<<0 && (types[i] & 1<<2) == 1<<2)
{
f_temp_4[i] = f_temp_6[i];
f_temp_8[i] = f_temp_6[i];
}
}
}

