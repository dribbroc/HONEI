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

#ifdef __CPU__

__kernel void dot_product_float(__global  float * x,
                                          __global float * y,
                                          __global float * tmp,
                                          const unsigned int size)
{
        uint tid = get_global_id(0);
        uint blocksize = get_local_size(0);

            // calculate how many elements each thread needs to calculate
            const unsigned int iter = size / (get_global_size(0));
            unsigned int pos =  tid;

            // clear the output
            tmp[pos] = 0;

            for (unsigned long i = 0 ; i < iter ; ++i)
            {
                tmp[tid] += x[pos] * y[pos];
                pos += get_global_size(0);
            }

            // for the last iteration, check if the elements are still available
            if (pos < size)
            {
                tmp[tid] += x[pos] * y[pos];
            }
}

#pragma OPENCL EXTENSION cl_amd_fp64 : enable
__kernel void dot_product_double(__global  double * x,
                                          __global double * y,
                                          __global double * tmp,
                                          const unsigned int size)
{
        uint tid = get_global_id(0);
        uint blocksize = get_local_size(0);

            // calculate how many elements each thread needs to calculate
            const unsigned int iter = size / (get_global_size(0));
            unsigned int pos =  tid;

            // clear the output
            tmp[pos] = 0;

            for (unsigned long i = 0 ; i < iter ; ++i)
            {
                tmp[tid] += x[pos] * y[pos];
                pos += get_global_size(0);
            }

            // for the last iteration, check if the elements are still available
            if (pos < size)
            {
                tmp[tid] += x[pos] * y[pos];
            }
}
//#endif

#else
//#ifdef __GPU__
__kernel void dot_product_float(__global  float * x,
                                          __global float * y,
                                          __global float * tmp,
                                          const unsigned int size)
{
        uint tid = get_global_id(0);
        uint blocksize = get_local_size(0);

            // calculate how many elements each thread needs to calculate
            const unsigned int iter = size / (get_global_size(0));
            unsigned int pos =  tid;

            // clear the output
            tmp[pos] = 0;

            for (unsigned long i = 0 ; i < iter ; ++i)
            {
                tmp[tid] += x[pos] * y[pos];
                pos += get_global_size(0);
            }

            // for the last iteration, check if the elements are still available
            if (pos < size)
            {
                tmp[tid] += x[pos] * y[pos];
            }
}


#ifdef HONEI_CUDA_DOUBLE
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
__kernel void dot_product_double(__global  double * x,
                                          __global double * y,
                                          __global double * tmp,
                                          const unsigned int size)
{
        uint tid = get_global_id(0);
        uint blocksize = get_local_size(0);

            // calculate how many elements each thread needs to calculate
            const unsigned int iter = size / (get_global_size(0));
            unsigned int pos =  tid;

            // clear the output
            tmp[pos] = 0;

            for (unsigned long i = 0 ; i < iter ; ++i)
            {
                tmp[tid] += x[pos] * y[pos];
                pos += get_global_size(0);
            }

            // for the last iteration, check if the elements are still available
            if (pos < size)
            {
                tmp[tid] += x[pos] * y[pos];
            }
}
#endif
#endif
