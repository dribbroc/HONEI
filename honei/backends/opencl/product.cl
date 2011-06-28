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
__kernel void product_smell_dv_f(__global float * x, __global float * y, __global unsigned long * Aj, __global float * Ax,
                                              __global unsigned long * Arl, unsigned long num_rows, unsigned long stride, unsigned long threads,
                                              __local float * shared_ell_float)
{
    uint idx = get_global_id(0);
    uint idb = get_local_id(0);
    uint row = idx;

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

__kernel void product_smell_dv_d(__global double * x, __global double * y, __global unsigned long * Aj, __global double * Ax,
                                              __global unsigned long * Arl, unsigned long num_rows, unsigned long stride, unsigned long threads,
                                              __local double * shared_ell_double)
{
    uint idx = get_global_id(0);
    uint idb = get_local_id(0);
    uint row = idx;

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
__kernel void product_smell_dv_f(__global float * x, __global float * y, __global unsigned long * Aj, __global float * Ax,
                                              __global unsigned long * Arl, unsigned long num_rows, unsigned long stride, unsigned long threads,
                                              __local float * shared_ell_float)
{
    uint T = threads;
    uint idx = get_global_id(0);
    uint idb = get_local_id(0);
    uint idp = idb % T;
    uint row = idx / T;

            if(row >= num_rows){ return; }
            shared_ell_float[idb] = 0;
            float sum = 0;

            const unsigned long max = Arl[row];
            Ax += (row*T)+idp;
            Aj += (row*T)+idp;
            for(unsigned long k = 0; k < max ; ++k)
            {
                sum += *Ax * x[*Aj];
                Ax += stride;
                Aj += stride;
            }
            shared_ell_float[idb] = sum;

            switch (threads)
            {
                case 32:
                    if (idp < 16)
                        shared_ell_float[idb] += shared_ell_float[idb + 16];
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

__kernel void product_smell_dv_d(__global double * x, __global double * y, __global unsigned long * Aj, __global double * Ax,
                                              __global unsigned long * Arl, unsigned long num_rows, unsigned long stride, unsigned long threads,
                                              __local double * shared_ell_double)
{
    uint T = threads;
    uint idx = get_global_id(0);
    uint idb = get_local_id(0);
    uint idp = idb % T;
    uint row = idx / T;

            if(row >= num_rows){ return; }
            shared_ell_double[idb] = 0;
            double sum = 0;

            const unsigned long max = Arl[row];
            Ax += (row*T)+idp;
            Aj += (row*T)+idp;
            for(unsigned long k = 0; k < max ; ++k)
            {
                sum += *Ax * x[*Aj];
                Ax += stride;
                Aj += stride;
            }
            shared_ell_double[idb] = sum;

            switch (threads)
            {
                case 32:
                    if (idp < 16)
                        shared_ell_double[idb] += shared_ell_double[idb + 16];
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
#endif

#ifdef __CPU__
__kernel void product_bmdv_q1_f(__global float* ll, __global float* ld, __global float* lu,
                                    __global float* dl, __global float * dd, __global float* du,
                                    __global float* ul, __global float* ud, __global float* uu,
                                    __global float * x, __global float * y, unsigned long n, unsigned long m,
                                    __local float * smvf_cache)
{

unsigned long idx = get_global_id(0);

// runs from 0 to blockDim.x-1
unsigned long lindex = get_local_id(0);

__local float* Dcache = smvf_cache;
__local float* Lcache = smvf_cache + get_local_size(0) + 2;
__local float* Ucache = smvf_cache + 2 * (get_local_size(0) + 2);

// prefetch chunks from iteration vector
//
//
// data needed for DD, DU, DL: each thread loads one element, the first and last one load the border cases
// x_0 ... x_blockdim-1 into c_1...c_blockdim
Dcache[lindex + 1] = x[idx];
if (idx  >= m) Lcache[lindex + 1] = x[idx - m];
if (idx + m < n) Ucache[lindex + 1] = x[idx + m];
if (lindex == 0)
{
// x_-1 in c_0
if (get_local_size(0) * get_group_id(0) - 1 < n) Dcache[0] = x[get_local_size(0) * get_group_id(0) - 1];
if (get_local_size(0) * get_group_id(0) - m - 1 < n) Lcache[0] = x[get_local_size(0) * get_group_id(0) - m - 1];
if (get_local_size(0) * get_group_id(0) + m - 1 < n) Ucache[0] = x[get_local_size(0) * get_group_id(0) + m - 1];
}
if (lindex == get_local_size(0) - 1)
{
// x_blockdim in c_blockdim+1
if (get_local_size(0) * (get_group_id(0) + 1) < n) Dcache[get_local_size(0) + 1] = x[get_local_size(0) * (get_group_id(0) + 1)];
if (get_local_size(0) * (get_group_id(0) + 1) - m < n) Lcache[get_local_size(0) + 1] = x[get_local_size(0) * (get_group_id(0) + 1) - m];
if (get_local_size(0) * (get_group_id(0) + 1) + m  < n) Ucache[get_local_size(0) + 1] = x[get_local_size(0) * (get_group_id(0) + 1) + m];
}

// now, compute
if (idx < n)
{
float ytemp1 = dd[idx] * x[idx];
if (idx > 0) ytemp1 += dl[idx] * x[idx-1];
if (idx < n - 1) ytemp1 += du[idx] * x[idx+1];

if (idx > m) ytemp1 += ll[idx] * x[idx-m-1];
if (idx > m - 1) ytemp1 += ld[idx] * x[idx-m];
if (idx > m - 2) ytemp1 += lu[idx] * x[idx-m+1];

if (idx < n - m + 1) ytemp1 += ul[idx] * x[idx + m - 1];
if (idx < n - m) ytemp1 += ud[idx] * x[idx + m];;
if (idx < n - m - 1) ytemp1 += uu[idx] * x[idx + m + 1];
y[idx] = ytemp1;
}
}

__kernel void product_bmdv_q1_d(__global double* ll, __global double* ld, __global double* lu,
                                    __global double* dl, __global double * dd, __global double* du,
                                    __global double* ul, __global double* ud, __global double* uu,
                                    __global double * x, __global double * y, unsigned long n, unsigned long m,
                                    __local double * smvf_cache)
{

unsigned long idx = get_global_id(0);

// runs from 0 to blockDim.x-1
unsigned long lindex = get_local_id(0);

__local double* Dcache = smvf_cache;
__local double* Lcache = smvf_cache + get_local_size(0) + 2;
__local double* Ucache = smvf_cache + 2 * (get_local_size(0) + 2);

// prefetch chunks from iteration vector
//
//
// data needed for DD, DU, DL: each thread loads one element, the first and last one load the border cases
// x_0 ... x_blockdim-1 into c_1...c_blockdim
Dcache[lindex + 1] = x[idx];
if (idx  >= m) Lcache[lindex + 1] = x[idx - m];
if (idx + m < n) Ucache[lindex + 1] = x[idx + m];
if (lindex == 0)
{
// x_-1 in c_0
if (get_local_size(0) * get_group_id(0) - 1 < n) Dcache[0] = x[get_local_size(0) * get_group_id(0) - 1];
if (get_local_size(0) * get_group_id(0) - m - 1 < n) Lcache[0] = x[get_local_size(0) * get_group_id(0) - m - 1];
if (get_local_size(0) * get_group_id(0) + m - 1 < n) Ucache[0] = x[get_local_size(0) * get_group_id(0) + m - 1];
}
if (lindex == get_local_size(0) - 1)
{
// x_blockdim in c_blockdim+1
if (get_local_size(0) * (get_group_id(0) + 1) < n) Dcache[get_local_size(0) + 1] = x[get_local_size(0) * (get_group_id(0) + 1)];
if (get_local_size(0) * (get_group_id(0) + 1) - m < n) Lcache[get_local_size(0) + 1] = x[get_local_size(0) * (get_group_id(0) + 1) - m];
if (get_local_size(0) * (get_group_id(0) + 1) + m  < n) Ucache[get_local_size(0) + 1] = x[get_local_size(0) * (get_group_id(0) + 1) + m];
}

// now, compute
if (idx < n)
{
double ytemp1 = dd[idx] * x[idx];
if (idx > 0) ytemp1 += dl[idx] * x[idx-1];
if (idx < n - 1) ytemp1 += du[idx] * x[idx+1];

if (idx > m) ytemp1 += ll[idx] * x[idx-m-1];
if (idx > m - 1) ytemp1 += ld[idx] * x[idx-m];
if (idx > m - 2) ytemp1 += lu[idx] * x[idx-m+1];

if (idx < n - m + 1) ytemp1 += ul[idx] * x[idx + m - 1];
if (idx < n - m) ytemp1 += ud[idx] * x[idx + m];;
if (idx < n - m - 1) ytemp1 += uu[idx] * x[idx + m + 1];
y[idx] = ytemp1;
}
}

#else

__kernel void product_bmdv_q1_f(__global float* ll, __global float* ld, __global float* lu,
                                    __global float* dl, __global float * dd, __global float* du,
                                    __global float* ul, __global float* ud, __global float* uu,
                                    __global float * x, __global float * y, unsigned long n, unsigned long m,
                                    __local float * smvf_cache)
{

unsigned long idx = get_global_id(0);

// runs from 0 to blockDim.x-1
unsigned long lindex = get_local_id(0);

__local float* Dcache = smvf_cache;
__local float* Lcache = smvf_cache + get_local_size(0) + 2;
__local float* Ucache = smvf_cache + 2 * (get_local_size(0) + 2);

// prefetch chunks from iteration vector
//
//
// data needed for DD, DU, DL: each thread loads one element, the first and last one load the border cases
// x_0 ... x_blockdim-1 into c_1...c_blockdim
Dcache[lindex + 1] = x[idx];
if (idx  >= m) Lcache[lindex + 1] = x[idx - m];
if (idx + m < n) Ucache[lindex + 1] = x[idx + m];
if (lindex == 0)
{
// x_-1 in c_0
if (get_local_size(0) * get_group_id(0) - 1 < n) Dcache[0] = x[get_local_size(0) * get_group_id(0) - 1];
if (get_local_size(0) * get_group_id(0) - m - 1 < n) Lcache[0] = x[get_local_size(0) * get_group_id(0) - m - 1];
if (get_local_size(0) * get_group_id(0) + m - 1 < n) Ucache[0] = x[get_local_size(0) * get_group_id(0) + m - 1];
}
if (lindex == get_local_size(0) - 1)
{
// x_blockdim in c_blockdim+1
if (get_local_size(0) * (get_group_id(0) + 1) < n) Dcache[get_local_size(0) + 1] = x[get_local_size(0) * (get_group_id(0) + 1)];
if (get_local_size(0) * (get_group_id(0) + 1) - m < n) Lcache[get_local_size(0) + 1] = x[get_local_size(0) * (get_group_id(0) + 1) - m];
if (get_local_size(0) * (get_group_id(0) + 1) + m  < n) Ucache[get_local_size(0) + 1] = x[get_local_size(0) * (get_group_id(0) + 1) + m];
}

mem_fence(CLK_LOCAL_MEM_FENCE);

// now, compute
if (idx < n)
{
float ytemp1 = dd[idx] * Dcache[lindex + 1];
if (idx > 0) ytemp1 += dl[idx] * Dcache[lindex];
if (idx < n - 1) ytemp1 += du[idx] * Dcache[lindex + 2];

if (idx > m) ytemp1 += ll[idx] * Lcache[lindex];
if (idx > m - 1) ytemp1 += ld[idx] * Lcache[lindex + 1];
if (idx > m - 2) ytemp1 += lu[idx] * Lcache[lindex + 2];

if (idx < n - m + 1) ytemp1 += ul[idx] * Ucache[lindex];
if (idx < n - m) ytemp1 += ud[idx] * Ucache[lindex + 1];
if (idx < n - m - 1) ytemp1 += uu[idx] * Ucache[lindex + 2];
y[idx] = ytemp1;
}
}

__kernel void product_bmdv_q1_d(__global double* ll, __global double* ld, __global double* lu,
                                    __global double* dl, __global double * dd, __global double* du,
                                    __global double* ul, __global double* ud, __global double* uu,
                                    __global double * x, __global double * y, unsigned long n, unsigned long m,
                                    __local double * smvf_cache)
{

unsigned long idx = get_global_id(0);

// runs from 0 to blockDim.x-1
unsigned long lindex = get_local_id(0);

__local double* Dcache = smvf_cache;
__local double* Lcache = smvf_cache + get_local_size(0) + 2;
__local double* Ucache = smvf_cache + 2 * (get_local_size(0) + 2);

// prefetch chunks from iteration vector
//
//
// data needed for DD, DU, DL: each thread loads one element, the first and last one load the border cases
// x_0 ... x_blockdim-1 into c_1...c_blockdim
Dcache[lindex + 1] = x[idx];
if (idx  >= m) Lcache[lindex + 1] = x[idx - m];
if (idx + m < n) Ucache[lindex + 1] = x[idx + m];
if (lindex == 0)
{
// x_-1 in c_0
if (get_local_size(0) * get_group_id(0) - 1 < n) Dcache[0] = x[get_local_size(0) * get_group_id(0) - 1];
if (get_local_size(0) * get_group_id(0) - m - 1 < n) Lcache[0] = x[get_local_size(0) * get_group_id(0) - m - 1];
if (get_local_size(0) * get_group_id(0) + m - 1 < n) Ucache[0] = x[get_local_size(0) * get_group_id(0) + m - 1];
}
if (lindex == get_local_size(0) - 1)
{
// x_blockdim in c_blockdim+1
if (get_local_size(0) * (get_group_id(0) + 1) < n) Dcache[get_local_size(0) + 1] = x[get_local_size(0) * (get_group_id(0) + 1)];
if (get_local_size(0) * (get_group_id(0) + 1) - m < n) Lcache[get_local_size(0) + 1] = x[get_local_size(0) * (get_group_id(0) + 1) - m];
if (get_local_size(0) * (get_group_id(0) + 1) + m  < n) Ucache[get_local_size(0) + 1] = x[get_local_size(0) * (get_group_id(0) + 1) + m];
}

mem_fence(CLK_LOCAL_MEM_FENCE);

// now, compute
if (idx < n)
{
double ytemp1 = dd[idx] * Dcache[lindex + 1];
if (idx > 0) ytemp1 += dl[idx] * Dcache[lindex];
if (idx < n - 1) ytemp1 += du[idx] * Dcache[lindex + 2];

if (idx > m) ytemp1 += ll[idx] * Lcache[lindex];
if (idx > m - 1) ytemp1 += ld[idx] * Lcache[lindex + 1];
if (idx > m - 2) ytemp1 += lu[idx] * Lcache[lindex + 2];

if (idx < n - m + 1) ytemp1 += ul[idx] * Ucache[lindex];
if (idx < n - m) ytemp1 += ud[idx] * Ucache[lindex + 1];
if (idx < n - m - 1) ytemp1 += uu[idx] * Ucache[lindex + 2];
y[idx] = ytemp1;
}
}
#endif
