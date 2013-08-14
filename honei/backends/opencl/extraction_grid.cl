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

__kernel void extraction_grid_dry_f(const unsigned long start, const unsigned long end,
                                      __global float * f_0, __global float * f_1, __global float * f_2,
                                      __global float * f_3, __global float * f_4, __global float * f_5,
                                      __global float * f_6, __global float * f_7, __global float * f_8,
                                      __global float * h, __global float * u, __global float * v,
                                      __global float * distribution_x_data, __global float * distribution_y_data, const float epsilon,
                                      __local float * distribution_cache)
{
__local float* distribution_x = distribution_cache;
__local float* distribution_y = distribution_cache + 9;

uint idx = get_global_id(0);

uint midx = idx % get_local_size(0);
if (midx < 9)
{
distribution_x[midx] = distribution_x_data[midx];
distribution_y[midx] = distribution_y_data[midx];
}
barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

idx += start;
if (idx < end)
{
unsigned long i = idx;

h[i] = f_0[i] +
f_1[i] +
f_2[i] +
f_3[i] +
f_4[i] +
f_5[i] +
f_6[i] +
f_7[i] +
f_8[i];

float lax_upper = epsilon;
float lax_lower = -lax_upper;

if(h[i] < lax_lower || h[i] > lax_upper)
{
u[i] = (distribution_x[0] * f_0[i] +
                          distribution_x[1] * f_1[i] +
                          distribution_x[2] * f_2[i] +
                          distribution_x[3] * f_3[i] +
                          distribution_x[4] * f_4[i] +
                          distribution_x[5] * f_5[i] +
                          distribution_x[6] * f_6[i] +
                          distribution_x[7] * f_7[i] +
                          distribution_x[8] * f_8[i]) / h[i];

v[i] = (distribution_y[0] * f_0[i] +
                          distribution_y[1] * f_1[i] +
                          distribution_y[2] * f_2[i] +
                          distribution_y[3] * f_3[i] +
                          distribution_y[4] * f_4[i] +
                          distribution_y[5] * f_5[i] +
                          distribution_y[6] * f_6[i] +
                          distribution_y[7] * f_7[i] +
                          distribution_y[8] * f_8[i]) / h[i];
}
else
{
h[i] = 0;
u[i] = 0;
v[i] = 0;
}
h[i] = max(0.f, min(1.f, h[i]));
}
}

__kernel void extraction_grid_dry_d(const unsigned long start, const unsigned long end,
                                      __global double * f_0, __global double * f_1, __global double * f_2,
                                      __global double * f_3, __global double * f_4, __global double * f_5,
                                      __global double * f_6, __global double * f_7, __global double * f_8,
                                      __global double * h, __global double * u, __global double * v,
                                      __global double * distribution_x_data, __global double * distribution_y_data, const double epsilon,
                                      __local double * distribution_cache)
{
__local double* distribution_x = distribution_cache;
__local double* distribution_y = distribution_cache + 9;

uint idx = get_global_id(0);

uint midx = idx % get_local_size(0);
if (midx < 9)
{
distribution_x[midx] = distribution_x_data[midx];
distribution_y[midx] = distribution_y_data[midx];
}
barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

idx += start;
if (idx < end)
{
unsigned long i = idx;

h[i] = f_0[i] +
f_1[i] +
f_2[i] +
f_3[i] +
f_4[i] +
f_5[i] +
f_6[i] +
f_7[i] +
f_8[i];

double lax_upper = epsilon;
double lax_lower = -lax_upper;

if(h[i] < lax_lower || h[i] > lax_upper)
{
u[i] = (distribution_x[0] * f_0[i] +
                          distribution_x[1] * f_1[i] +
                          distribution_x[2] * f_2[i] +
                          distribution_x[3] * f_3[i] +
                          distribution_x[4] * f_4[i] +
                          distribution_x[5] * f_5[i] +
                          distribution_x[6] * f_6[i] +
                          distribution_x[7] * f_7[i] +
                          distribution_x[8] * f_8[i]) / h[i];

v[i] = (distribution_y[0] * f_0[i] +
                          distribution_y[1] * f_1[i] +
                          distribution_y[2] * f_2[i] +
                          distribution_y[3] * f_3[i] +
                          distribution_y[4] * f_4[i] +
                          distribution_y[5] * f_5[i] +
                          distribution_y[6] * f_6[i] +
                          distribution_y[7] * f_7[i] +
                          distribution_y[8] * f_8[i]) / h[i];
}
else
{
h[i] = 0;
u[i] = 0;
v[i] = 0;
}
h[i] = max(0., min(1., h[i]));
}
}


__kernel void extraction_grid_wet_f(const unsigned long start, const unsigned long end,
                                      __global float * f_0, __global float * f_1, __global float * f_2,
                                      __global float * f_3, __global float * f_4, __global float * f_5,
                                      __global float * f_6, __global float * f_7, __global float * f_8,
                                      __global float * h, __global float * u, __global float * v,
                                      __global float * distribution_x_data, __global float * distribution_y_data, const float epsilon,
                                      __local float * distribution_cache)
{
__local float* distribution_x = distribution_cache;
__local float* distribution_y = distribution_cache + 9;

uint idx = get_global_id(0);

uint midx = idx % get_local_size(0);
if (midx < 9)
{
distribution_x[midx] = distribution_x_data[midx];
distribution_y[midx] = distribution_y_data[midx];
}
barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

idx += start;
if (idx < end)
{
unsigned long i = idx;

h[i] = f_0[i] +
f_1[i] +
f_2[i] +
f_3[i] +
f_4[i] +
f_5[i] +
f_6[i] +
f_7[i] +
f_8[i];

u[i] = (distribution_x[0] * f_0[i] +
                          distribution_x[1] * f_1[i] +
                          distribution_x[2] * f_2[i] +
                          distribution_x[3] * f_3[i] +
                          distribution_x[4] * f_4[i] +
                          distribution_x[5] * f_5[i] +
                          distribution_x[6] * f_6[i] +
                          distribution_x[7] * f_7[i] +
                          distribution_x[8] * f_8[i]) / h[i];

v[i] = (distribution_y[0] * f_0[i] +
                          distribution_y[1] * f_1[i] +
                          distribution_y[2] * f_2[i] +
                          distribution_y[3] * f_3[i] +
                          distribution_y[4] * f_4[i] +
                          distribution_y[5] * f_5[i] +
                          distribution_y[6] * f_6[i] +
                          distribution_y[7] * f_7[i] +
                          distribution_y[8] * f_8[i]) / h[i];
}
}

__kernel void extraction_grid_wet_d(const unsigned long start, const unsigned long end,
                                      __global double * f_0, __global double * f_1, __global double * f_2,
                                      __global double * f_3, __global double * f_4, __global double * f_5,
                                      __global double * f_6, __global double * f_7, __global double * f_8,
                                      __global double * h, __global double * u, __global double * v,
                                      __global double * distribution_x_data, __global double * distribution_y_data, const double epsilon,
                                      __local double * distribution_cache)
{
__local double* distribution_x = distribution_cache;
__local double* distribution_y = distribution_cache + 9;

uint idx = get_global_id(0);

uint midx = idx % get_local_size(0);
if (midx < 9)
{
distribution_x[midx] = distribution_x_data[midx];
distribution_y[midx] = distribution_y_data[midx];
}
barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

idx += start;
if (idx < end)
{
unsigned long i = idx;

h[i] = f_0[i] +
f_1[i] +
f_2[i] +
f_3[i] +
f_4[i] +
f_5[i] +
f_6[i] +
f_7[i] +
f_8[i];

u[i] = (distribution_x[0] * f_0[i] +
                          distribution_x[1] * f_1[i] +
                          distribution_x[2] * f_2[i] +
                          distribution_x[3] * f_3[i] +
                          distribution_x[4] * f_4[i] +
                          distribution_x[5] * f_5[i] +
                          distribution_x[6] * f_6[i] +
                          distribution_x[7] * f_7[i] +
                          distribution_x[8] * f_8[i]) / h[i];

v[i] = (distribution_y[0] * f_0[i] +
                          distribution_y[1] * f_1[i] +
                          distribution_y[2] * f_2[i] +
                          distribution_y[3] * f_3[i] +
                          distribution_y[4] * f_4[i] +
                          distribution_y[5] * f_5[i] +
                          distribution_y[6] * f_6[i] +
                          distribution_y[7] * f_7[i] +
                          distribution_y[8] * f_8[i]) / h[i];
}
}
