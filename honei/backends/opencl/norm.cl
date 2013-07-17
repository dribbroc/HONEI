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

#ifndef __CPU__
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif

#ifdef __CPU__
__kernel void norm_l2_false_f(__global  float * x,
                                          __global float * tmp,
                                          const unsigned long size)
{
        uint tid = get_global_id(0);
        uint blocksize = get_local_size(0);

            // calculate how many elements each thread needs to calculate
            const unsigned long iter = size / (get_global_size(0));
            unsigned long pos =  tid * iter;

            // clear the output
            float temp = 0;

            for (unsigned long i = 0 ; i < iter ; ++i)
            {
                temp += x[pos + i] * x[pos + i];
            }

            // for the last iteration, check if the elements are still available
            if (tid == get_global_size(0) - 1)
            {
                pos += iter;
                while (pos < size)
                {
                    temp += x[pos] * x[pos];
                    ++pos;
                }
            }
            tmp[tid] = temp;
}

__kernel void norm_l2_false_d(__global  double * x,
                                          __global double * tmp,
                                          const unsigned long size)
{
        uint tid = get_global_id(0);
        uint blocksize = get_local_size(0);

            // calculate how many elements each thread needs to calculate
            const unsigned long iter = size / (get_global_size(0));
            unsigned long pos =  tid * iter;

            // clear the output
            double temp = 0;

            for (unsigned long i = 0 ; i < iter ; ++i)
            {
                temp += x[pos + i] * x[pos + i];
            }

            // for the last iteration, check if the elements are still available
            if (tid == get_global_size(0) - 1)
            {
                pos += iter;
                while (pos < size)
                {
                    temp += x[pos] * x[pos];
                    ++pos;
                }
            }
            tmp[tid] = temp;
}
#else
__kernel void norm_l2_false_f(__global  float * x,
                                          __global float * tmp,
                                          const unsigned long size)
{
        uint tid = get_global_id(0);
        uint blocksize = get_local_size(0);

            // calculate how many elements each thread needs to calculate
            const unsigned long iter = size / (get_global_size(0));
            unsigned long pos =  tid;

            // clear the output
            float temp = 0;

            for (unsigned long i = 0 ; i < iter ; ++i)
            {
                temp += x[pos] * x[pos];
                pos += get_global_size(0);
            }

            // for the last iteration, check if the elements are still available
            if (pos < size)
            {
                temp += x[pos] * x[pos];
            }
            tmp[tid] = temp;
}

__kernel void norm_l2_false_d(__global  double * x,
                                          __global double * tmp,
                                          const unsigned long size)
{
        uint tid = get_global_id(0);
        uint blocksize = get_local_size(0);

            // calculate how many elements each thread needs to calculate
            const unsigned long iter = size / (get_global_size(0));
            unsigned long pos =  tid;

            // clear the output
            double temp = 0;

            for (unsigned long i = 0 ; i < iter ; ++i)
            {
                temp += x[pos] * x[pos];
                pos += get_global_size(0);
            }

            // for the last iteration, check if the elements are still available
            if (pos < size)
            {
                temp += x[pos] * x[pos];
            }
            tmp[tid] = temp;
}
#endif
