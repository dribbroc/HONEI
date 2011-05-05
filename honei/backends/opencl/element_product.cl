/* vim: set sw=4 sts=4 et nofoldenable filetype=cpp : */

/*
* Copyright (c) 2011 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
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

__kernel void element_product_three_float(__global  float * output,
                                  __global const float * x,
                                  __global const float * y,
                                  const unsigned int size)
{
uint tid = get_global_id(0);

    if (tid < size) output[tid] = x[tid] * y[tid];
}

#pragma OPENCL EXTENSION cl_khr_fp64 : enable
__kernel void element_product_three_double(__global  double * output,
                                  __global const double * x,
                                  __global const double * y,
                                  const unsigned int size)
{
uint tid = get_global_id(0);

    if (tid < size) output[tid] = x[tid] * y[tid];
}
