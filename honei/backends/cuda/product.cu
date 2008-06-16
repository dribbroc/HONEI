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
            Lcache[lindex + 1] = x[idx - m];
            Ucache[lindex + 1] = x[idx + m];
            if (lindex == 0)
            {
                // x_-1 in c_0
                Dcache[0] = x[blockDim.x * blockIdx.x - 1];
                Lcache[0] = x[blockDim.x * blockIdx.x - m - 1];
                Ucache[0] = x[blockDim.x * blockIdx.x + m - 1];
            }
            if (lindex == blockDim.x - 1)
            {
                // x_blockdim in c_blockdim+1
                Dcache[blockDim.x + 1] = x[blockDim.x * (blockIdx.x + 1)];
                Lcache[blockDim.x + 1] = x[blockDim.x * (blockIdx.x + 1) - m];
                Ucache[blockDim.x + 1] = x[blockDim.x * (blockIdx.x + 1) + m];
            }
            __syncthreads();
            // now, compute
            if (idx < n)
            {
                float ytemp1 = dd[idx] * Dcache[lindex + 1];
                if (idx > 0) ytemp1 += dl[idx - 1] * Dcache[lindex];
                if (idx < n - 1) ytemp1 += du[idx] * Dcache[lindex + 2];

                if (idx > m) ytemp1 += ll[idx - m - 1] * Lcache[lindex];
                if (idx > m - 1) ytemp1 += ld[idx - m] * Lcache[lindex + 1];
                if (idx > m - 2) ytemp1 += lu[idx - m + 1] * Lcache[lindex + 2];

                if (idx < n - m + 1) ytemp1 += ul[idx] * Ucache[lindex];
                if (idx < n - m) ytemp1 += ud[idx] * Ucache[lindex + 1];
                if (idx < n - m - 1) ytemp1 += uu[idx] * Ucache[lindex + 2];
                y[idx] = ytemp1;
            }
        }
    }
}

extern "C" void cuda_product_bmdv_q1_float (float * ll, float * ld, float * lu,
        float * dl, float * dd, float *du,
        float * ul, float * ud, float *uu, float * x, float * y,
        unsigned long size, unsigned long blocksize, unsigned long m)
{
    dim3 grid;
    dim3 block;
    block.x = blocksize;
    grid.x = (int)ceil(size/(double)(block.x));
    float * x_gpu(0);
    float * y_gpu(0);
    float * ll_gpu;
    float * ld_gpu;
    float * lu_gpu;
    float * dl_gpu;
    float * dd_gpu;
    float * du_gpu;
    float * ul_gpu;
    float * ud_gpu;
    float * uu_gpu;

    //todo: malloc and transfer only non zero padded parts
    cudaMalloc((void**)&x_gpu, size * sizeof(float));
    cudaMalloc((void**)&y_gpu, size * sizeof(float));
    cudaMalloc((void**)&ll_gpu, size * sizeof(float));
    cudaMalloc((void**)&ld_gpu, size * sizeof(float));
    cudaMalloc((void**)&lu_gpu, size * sizeof(float));
    cudaMalloc((void**)&dl_gpu, size * sizeof(float));
    cudaMalloc((void**)&dd_gpu, size * sizeof(float));
    cudaMalloc((void**)&du_gpu, size * sizeof(float));
    cudaMalloc((void**)&ul_gpu, size * sizeof(float));
    cudaMalloc((void**)&ud_gpu, size * sizeof(float));
    cudaMalloc((void**)&uu_gpu, size * sizeof(float));

    cudaMemcpy(x_gpu, x, size * sizeof(float), cudaMemcpyHostToDevice);
    //cudaMemcpy(y_gpu, y, size * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(ll_gpu, ll, (size - m - 1)* sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(ld_gpu, ld, (size - m ) * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(lu_gpu, lu, (size -m + 1) * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dl_gpu, dl, (size - 1) * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dd_gpu, dd, size * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(du_gpu, du, (size - 1) * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(ul_gpu, ul, (size - m + 1) * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(ud_gpu, ud, (size - m) * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(uu_gpu, uu, (size - m - 1) * sizeof(float), cudaMemcpyHostToDevice);

    honei::cuda::product_bmdv_q1_gpu<<<grid, block, 3 * (block.x + 2 ) * sizeof(float)>>>(ll_gpu, ld_gpu, lu_gpu, dl_gpu, dd_gpu, du_gpu, ul_gpu, ud_gpu, uu_gpu, x_gpu, y_gpu, size, m);
    cudaMemcpy(y, y_gpu, size * sizeof(float), cudaMemcpyDeviceToHost);
    cudaFree(x_gpu);
    cudaFree(y_gpu);
    cudaFree(ll_gpu);
    cudaFree(ld_gpu);
    cudaFree(lu_gpu);
    cudaFree(dl_gpu);
    cudaFree(dd_gpu);
    cudaFree(du_gpu);
    cudaFree(ul_gpu);
    cudaFree(ud_gpu);
    cudaFree(uu_gpu);

    CUDA_ERROR();
}
