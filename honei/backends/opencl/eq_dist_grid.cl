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

__kernel void eq_dist_grid_f(const unsigned long start, const unsigned long end,
                                      __global float * f_eq_0, __global float * f_eq_1, __global float * f_eq_2,
                                      __global float * f_eq_3, __global float * f_eq_4, __global float * f_eq_5,
                                      __global float * f_eq_6, __global float * f_eq_7, __global float * f_eq_8,
                                      __global float * h, __global float * u, __global float * v,
                                      __global float * distribution_x_data, __global float * distribution_y_data,
                                      const float g, const float e,
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
unsigned long index = idx;

float e2=(e);
float e23=(3. * e2);
float e26=(6. * e2);
float e42=(2. * e2 * e2);
float e48=(8. * e2 * e2);
float e212=(12. * e2);
float e224=(24. * e2);

float u2=(u[index] * u[index]);
float v2=(v[index] * v[index]);
float gh=(g * h[index]);

float dxu, dyv;
float t1, t2, t3, t4;

t1 = ((5.) * gh) / e26;
t2 = (2.) / e23 * (u2 + v2);
f_eq_0[index] = h[index] * ((1) - t1 - t2);

dxu = distribution_x[1] * u[index];
dyv = distribution_y[1] * v[index];
t1 = (gh) / e26;
t2 = (dxu + dyv) / e23;
t3 = (dxu * dxu + (2.) * dxu * dyv + dyv * dyv) / e42;
t4 = (u2 + v2) / e26;
f_eq_1[index] = h[index] * (t1 + t2 + t3 - t4);

dxu = distribution_x[3] * u[index];
dyv = distribution_y[3] * v[index];
t2 = (dxu + dyv) / e23;
t3 = (dxu * dxu + (2.) * dxu * dyv + dyv * dyv) / e42;
f_eq_3[index] = h[index] * (t1 + t2 + t3 - t4);

dxu = distribution_x[5] * u[index];
dyv = distribution_y[5] * v[index];
t2 = (dxu + dyv) / e23;
t3 = (dxu * dxu + (2.) * dxu * dyv + dyv * dyv) / e42;
f_eq_5[index] = h[index] * (t1 + t2 + t3 - t4);

dxu = distribution_x[7] * u[index];
dyv = distribution_y[7] * v[index];
t2 = (dxu + dyv) / e23;
t3 = (dxu * dxu + (2.) * dxu * dyv + dyv * dyv) / e42;
f_eq_7[index] = h[index] * (t1 + t2 + t3 - t4);

dxu = distribution_x[2] * u[index];
dyv = distribution_y[2] * v[index];
t1 = (gh) / e224;
t2 = (dxu + dyv) / e212;
t3 = (dxu * dxu + (2.) * dxu * dyv + dyv * dyv) / e48;
t4 = (u2 + v2) / e224;
f_eq_2[index] =  h[index] * (t1 + t2 + t3 - t4);

dxu = distribution_x[2] * u[index];
dyv = distribution_y[2] * v[index];
t2 = (dxu + dyv) / e212;
t3 = (dxu * dxu + (2.) * dxu * dyv + dyv * dyv) / e48;
f_eq_2[index] =  h[index] * (t1 + t2 + t3 - t4);

dxu = distribution_x[4] * u[index];
dyv = distribution_y[4] * v[index];
t2 = (dxu + dyv) / e212;
t3 = (dxu * dxu + (2.) * dxu * dyv + dyv * dyv) / e48;
f_eq_4[index] =  h[index] * (t1 + t2 + t3 - t4);

dxu = distribution_x[6] * u[index];
dyv = distribution_y[6] * v[index];
t2 = (dxu + dyv) / e212;
t3 = (dxu * dxu + (2.) * dxu * dyv + dyv * dyv) / e48;
f_eq_6[index] =  h[index] * (t1 + t2 + t3 - t4);

dxu = distribution_x[8] * u[index];
dyv = distribution_y[8] * v[index];
t2 = (dxu + dyv) / e212;
t3 = (dxu * dxu + (2.) * dxu * dyv + dyv * dyv) / e48;
f_eq_8[index] =  h[index] * (t1 + t2 + t3 - t4);
}
}

__kernel void eq_dist_grid_d(const unsigned long start, const unsigned long end,
                                      __global double * f_eq_0, __global double * f_eq_1, __global double * f_eq_2,
                                      __global double * f_eq_3, __global double * f_eq_4, __global double * f_eq_5,
                                      __global double * f_eq_6, __global double * f_eq_7, __global double * f_eq_8,
                                      __global double * h, __global double * u, __global double * v,
                                      __global double * distribution_x_data, __global double * distribution_y_data,
                                      const double g, const double e,
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
unsigned long index = idx;

double e2=(e);
double e23=(3. * e2);
double e26=(6. * e2);
double e42=(2. * e2 * e2);
double e48=(8. * e2 * e2);
double e212=(12. * e2);
double e224=(24. * e2);

double u2=(u[index] * u[index]);
double v2=(v[index] * v[index]);
double gh=(g * h[index]);

double dxu, dyv;
double t1, t2, t3, t4;

t1 = ((5.) * gh) / e26;
t2 = (2.) / e23 * (u2 + v2);
f_eq_0[index] = h[index] * ((1) - t1 - t2);

dxu = distribution_x[1] * u[index];
dyv = distribution_y[1] * v[index];
t1 = (gh) / e26;
t2 = (dxu + dyv) / e23;
t3 = (dxu * dxu + (2.) * dxu * dyv + dyv * dyv) / e42;
t4 = (u2 + v2) / e26;
f_eq_1[index] = h[index] * (t1 + t2 + t3 - t4);

dxu = distribution_x[3] * u[index];
dyv = distribution_y[3] * v[index];
t2 = (dxu + dyv) / e23;
t3 = (dxu * dxu + (2.) * dxu * dyv + dyv * dyv) / e42;
f_eq_3[index] = h[index] * (t1 + t2 + t3 - t4);

dxu = distribution_x[5] * u[index];
dyv = distribution_y[5] * v[index];
t2 = (dxu + dyv) / e23;
t3 = (dxu * dxu + (2.) * dxu * dyv + dyv * dyv) / e42;
f_eq_5[index] = h[index] * (t1 + t2 + t3 - t4);

dxu = distribution_x[7] * u[index];
dyv = distribution_y[7] * v[index];
t2 = (dxu + dyv) / e23;
t3 = (dxu * dxu + (2.) * dxu * dyv + dyv * dyv) / e42;
f_eq_7[index] = h[index] * (t1 + t2 + t3 - t4);

dxu = distribution_x[2] * u[index];
dyv = distribution_y[2] * v[index];
t1 = (gh) / e224;
t2 = (dxu + dyv) / e212;
t3 = (dxu * dxu + (2.) * dxu * dyv + dyv * dyv) / e48;
t4 = (u2 + v2) / e224;
f_eq_2[index] =  h[index] * (t1 + t2 + t3 - t4);

dxu = distribution_x[2] * u[index];
dyv = distribution_y[2] * v[index];
t2 = (dxu + dyv) / e212;
t3 = (dxu * dxu + (2.) * dxu * dyv + dyv * dyv) / e48;
f_eq_2[index] =  h[index] * (t1 + t2 + t3 - t4);

dxu = distribution_x[4] * u[index];
dyv = distribution_y[4] * v[index];
t2 = (dxu + dyv) / e212;
t3 = (dxu * dxu + (2.) * dxu * dyv + dyv * dyv) / e48;
f_eq_4[index] =  h[index] * (t1 + t2 + t3 - t4);

dxu = distribution_x[6] * u[index];
dyv = distribution_y[6] * v[index];
t2 = (dxu + dyv) / e212;
t3 = (dxu * dxu + (2.) * dxu * dyv + dyv * dyv) / e48;
f_eq_6[index] =  h[index] * (t1 + t2 + t3 - t4);

dxu = distribution_x[8] * u[index];
dyv = distribution_y[8] * v[index];
t2 = (dxu + dyv) / e212;
t3 = (dxu * dxu + (2.) * dxu * dyv + dyv * dyv) / e48;
f_eq_8[index] =  h[index] * (t1 + t2 + t3 - t4);
}
}
