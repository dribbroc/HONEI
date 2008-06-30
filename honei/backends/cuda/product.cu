/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
#include <honei/backends/cuda/cuda_util.hh>

namespace honei
{
    namespace cuda
    {
        //__global__ void product_bmdv_q1_gpu(float * x, float * y, unsigned long size)
        // further optimised version: don't compute on zeros in offdiagonals
        __global__ void product_bmdv_q1_gpu(
                float* ll, float* ld, float* lu,
                float* dl, float * dd, float* du,
                float* ul, float* ud, float* uu,
                float * x, float * y, int n, int m)
        {
            extern __shared__ float  smv_cache[];

            int idx = blockDim.x*blockIdx.x+threadIdx.x;

            // runs from 0 to blockDim.x-1
            int lindex = threadIdx.x;

            float* Dcache = smv_cache;
            float* Lcache = smv_cache + blockDim.x + 2;
            float* Ucache = smv_cache + 2 * (blockDim.x + 2);

            // prefetch chunks from iteration vector
            //
            //
            // data needed for DD, DU, DL: each thread loads one element, the first and last one load the border cases
            // x_0 ... x_blockdim-1 into c_1...c_blockdim
            Dcache[lindex + 1] = x[idx];
            if (idx - m >= 0) Lcache[lindex + 1] = x[idx - m];
            if (idx + m < n) Ucache[lindex + 1] = x[idx + m];
            if (lindex == 0)
            {
                // x_-1 in c_0
                if (blockDim.x * blockIdx.x - 1 >= 0 && blockDim.x * blockIdx.x - 1 < n) Dcache[0] = x[blockDim.x * blockIdx.x - 1];
                if (blockDim.x * blockIdx.x - m - 1 >= 0 && blockDim.x * blockIdx.x - m - 1 < n) Lcache[0] = x[blockDim.x * blockIdx.x - m - 1];
                if (blockDim.x * blockIdx.x + m - 1 >= 0 && blockDim.x * blockIdx.x + m - 1 < n) Ucache[0] = x[blockDim.x * blockIdx.x + m - 1];
            }
            if (lindex == blockDim.x - 1)
            {
                // x_blockdim in c_blockdim+1
                if (blockDim.x * (blockIdx.x +1) >= 0 && blockDim.x * (blockIdx.x +1) < n) Dcache[blockDim.x + 1] = x[blockDim.x * (blockIdx.x + 1)];
                if (blockDim.x * (blockIdx.x +1) - m >= 0 && blockDim.x * (blockIdx.x +1) - m < n) Lcache[blockDim.x + 1] = x[blockDim.x * (blockIdx.x + 1) - m];
                if (blockDim.x * (blockIdx.x +1) + m  >= 0 && blockDim.x * (blockIdx.x +1) + m  < n) Ucache[blockDim.x + 1] = x[blockDim.x * (blockIdx.x + 1) + m];
            }
            __syncthreads();
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
    }
}

extern "C" void cuda_product_bmdv_q1_float (void * ll, void * ld, void * lu,
        void * dl, void * dd, void *du,
        void * ul, void * ud, void *uu, void * x, void * y,
        unsigned long size, unsigned long blocksize, unsigned long m)
{
    dim3 grid;
    dim3 block;
    block.x = blocksize;
    grid.x = (int)ceil(size/(double)(block.x));
    float * x_gpu((float *)x);
    float * y_gpu((float *)y);
    float * ll_gpu((float *)ll);
    float * ld_gpu((float *)ld);
    float * lu_gpu((float *)lu);
    float * dl_gpu((float *)dl);
    float * dd_gpu((float *)dd);
    float * du_gpu((float *)du);
    float * ul_gpu((float *)ul);
    float * ud_gpu((float *)ud);
    float * uu_gpu((float *)uu);


    honei::cuda::product_bmdv_q1_gpu<<<grid, block, 3 * (block.x + 2 ) * sizeof(float)>>>(ll_gpu, ld_gpu, lu_gpu, dl_gpu, dd_gpu, du_gpu, ul_gpu, ud_gpu, uu_gpu, x_gpu, y_gpu, size, m);

    CUDA_ERROR();
}
