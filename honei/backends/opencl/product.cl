/* vim: set sw=4 sts=4 et nofoldenable filetype=cpp : */

/*
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#ifdef __CPU__
__kernel void product_smell_dv_float(__global float * x, __global float * y, __global unsigned long * Aj, __global float * Ax,
                                              __global unsigned long * Arl, unsigned long num_rows, unsigned long stride, unsigned long threads,
                                              __local float * shared_ell_float)
{
    uint idx = get_global_id(0);
    uint idb = get_local_id(0);
    uint row = idx;

            shared_ell_float[idb] = 0;
            if(row <  num_rows)
            {
                float sum = 0;

                uint max = Arl[row];
                //const unsigned long max = num_cols_per_row;
                for(uint k = 0; k < max ; ++k)
                {
                    for (uint thread = 0 ; thread < threads ; ++thread)
                    {
                        float value = Ax[k*stride+row*threads+ thread];
                        uint col = Aj[k*stride+row*threads + thread];
                        sum += value * x[col];
                    }
                }

                y[row] = sum;
            }
}

__kernel void product_smell_dv_double(__global double * x, __global double * y, __global unsigned long * Aj, __global double * Ax,
                                              __global unsigned long * Arl, unsigned long num_rows, unsigned long stride, unsigned long threads,
                                              __local double * shared_ell_double)
{
    uint idx = get_global_id(0);
    uint idb = get_local_id(0);
    uint row = idx;

            shared_ell_double[idb] = 0;
            if(row <  num_rows)
            {
                double sum = 0;

                uint max = Arl[row];
                //const unsigned long max = num_cols_per_row;
                for(uint k = 0; k < max ; ++k)
                {
                    for (uint thread = 0 ; thread < threads ; ++thread)
                    {
                        double value = Ax[k*stride+row*threads+ thread];
                        uint col = Aj[k*stride+row*threads + thread];
                        sum += value * x[col];
                    }
                }

                y[row] = sum;
            }
}
#else
__kernel void product_smell_dv_float(__global float * x, __global float * y, __global unsigned long * Aj, __global float * Ax,
                                              __global unsigned long * Arl, unsigned long num_rows, unsigned long stride, unsigned long threads,
                                              __local float * shared_ell_float)
{
    uint T = threads;
    uint idx = get_global_id(0);
    uint idb = get_local_id(0);
    uint idp = idb % T;
    uint row = idx / T;

            shared_ell_float[idb] = 0;
            if(row <  num_rows)
            {
                float sum = 0;

                uint max = Arl[row];
                //const unsigned long max = num_cols_per_row;
                for(uint k = 0; k < max ; ++k)
                {
                    float value = Ax[k*stride+row*T+idp];
                    if (value != 0)
                    {
                        uint col = Aj[k*stride+row*T+idp];
                        sum += value * x[col];
                    }
                }
                shared_ell_float[idb] = sum;
                mem_fence(CLK_LOCAL_MEM_FENCE);

                switch (threads)
                {
                    case 16:
                        if (idp < 8)
                            shared_ell_float[idb] += shared_ell_float[idb + 8];
                    case 8:
                        if (idp < 4)
                            shared_ell_float[idb] += shared_ell_float[idb + 4];
                    case 4:
                        if (idp < 2)
                            shared_ell_float[idb] += shared_ell_float[idb + 2];
                    case 2:
                        if (idp == 0)
                            y[row] = shared_ell_float[idb] + shared_ell_float[idb + 1];
                        break;
                    case 1:
                        y[row] = shared_ell_float[idb];
                        break;
                    default:
                        break;
                }
            }
}

__kernel void product_smell_dv_double(__global double * x, __global double * y, __global unsigned long * Aj, __global double * Ax,
                                              __global unsigned long * Arl, unsigned long num_rows, unsigned long stride, unsigned long threads,
                                              __local double * shared_ell_double)
{
    uint T = threads;
    uint idx = get_global_id(0);
    uint idb = get_local_id(0);
    uint idp = idb % T;
    uint row = idx / T;

            shared_ell_double[idb] = 0;
            if(row <  num_rows)
            {
                double sum = 0;

                const unsigned long max = Arl[row];
                //const unsigned long max = num_cols_per_row;
                for(unsigned long k = 0; k < max ; ++k)
                {
                    const double value = Ax[k*stride+row*T+idp];
                    if (value != 0)
                    {
                        const unsigned long col = Aj[k*stride+row*T+idp];
                        sum += value * x[col];
                    }
                }
                shared_ell_double[idb] = sum;

                switch (threads)
                {
                    case 16:
                        if (idp < 8)
                            shared_ell_double[idb] += shared_ell_double[idb + 8];
                    case 8:
                        if (idp < 4)
                            shared_ell_double[idb] += shared_ell_double[idb + 4];
                    case 4:
                        if (idp < 2)
                            shared_ell_double[idb] += shared_ell_double[idb + 2];
                    case 2:
                        if (idp == 0)
                            y[row] = shared_ell_double[idb] + shared_ell_double[idb + 1];
                        break;
                    case 1:
                        y[row] = shared_ell_double[idb];
                        break;
                    default:
                        break;
                }
            }
}
#endif
