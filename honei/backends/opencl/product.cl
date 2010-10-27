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

__kernel void product_smell_dv_float(__global float * x, __global float * y, __global unsigned int * Aj, __global float * Ax,
                                              __global unsigned int * Arl, unsigned int num_rows, unsigned int stride)
{
    uint row = get_global_id(0);

    if(row >= num_rows){ return; }
    float sum = 0.f;

    Aj += row;
    Ax += row;

    const unsigned int max = Arl[row];
    for(unsigned int n = 0; n < max ; n++){
        const float A_ij = *Ax;

        //if (A_ij != 0)
        {
            const unsigned int col = *Aj;
            sum += A_ij * x[col];
        }

        Aj += stride;
        Ax += stride;
    }

    y[row] = sum;
}

__kernel void product_smell_dv_double(__global double * x, __global double * y, __global unsigned int * Aj, __global double * Ax,
                                               __global unsigned int * Arl, unsigned int num_rows, unsigned int stride)
{
    uint row = get_global_id(0);

    if(row >= num_rows){ return; }
    double sum = 0.;

    Aj += row;
    Ax += row;

    const unsigned int max = Arl[row];
    for(unsigned int n = 0; n < max ; n++){
        const double A_ij = *Ax;

        //if (A_ij != 0)
        {
            const unsigned int col = *Aj;
            sum += A_ij * x[col];
        }

        Aj += stride;
        Ax += stride;
    }

    y[row] = sum;
}
