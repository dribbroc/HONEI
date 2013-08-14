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

#ifdef HONEI_OPENCL_GPU
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif

__kernel void collide_stream_grid_0_f(const unsigned long start, const unsigned long end,
                                      __global float * f_eq_0,
                                      __global float * f_0,
                                      __global float * f_temp_0,
                                      const float tau)
{
uint idx = get_global_id(0);

idx += start;
if (idx < end)
{
unsigned long i = idx;

f_temp_0[i] = f_0[i] - (f_0[i] - f_eq_0[i]) / tau;
}
}

__kernel void collide_stream_grid_n_f(const unsigned long start, const unsigned long end,
                                      __global unsigned long * dir,
                                      __global float * f_eq,
                                      __global float * f,
                                      __global float * f_temp,
                                      const float tau, const unsigned long size)
{
uint idx = get_global_id(0);

idx += start;
if (idx < end)
{
unsigned long i = idx;

if (dir[i] < size)
f_temp[dir[i]] = f[i] - (f[i] - f_eq[i])/tau;
}
}

__kernel void collide_stream_grid_0_d(const unsigned long start, const unsigned long end,
                                      __global double * f_eq_0,
                                      __global double * f_0,
                                      __global double * f_temp_0,
                                      const double tau)
{
uint idx = get_global_id(0);

idx += start;
if (idx < end)
{
unsigned long i = idx;

f_temp_0[i] = f_0[i] - (f_0[i] - f_eq_0[i]) / tau;
}
}

__kernel void collide_stream_grid_n_d(const unsigned long start, const unsigned long end,
                                      __global unsigned long * dir,
                                      __global double * f_eq,
                                      __global double * f,
                                      __global double * f_temp,
                                      const double tau, const unsigned long size)
{
uint idx = get_global_id(0);

idx += start;
if (idx < end)
{
unsigned long i = idx;

if (dir[i] < size)
f_temp[dir[i]] = f[i] - (f[i] - f_eq[i])/tau;
}
}
