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

// ceil(x/y) for integers, used to determine # of blocks/warps etc.
#define DIVIDE_INTO(x,y) ((x + y - 1)/y)
#define large_grid_thread_id(void) ((__umul24(blockDim.x,blockIdx.x + __umul24(blockIdx.y,gridDim.x)) + threadIdx.x))

dim3 make_large_grid_defect(const unsigned int num_threads, const unsigned int blocksize){
    const unsigned int num_blocks = DIVIDE_INTO(num_threads, blocksize);
    if (num_blocks <= 65535){
        //fits in a 1D grid
        return dim3(num_blocks);
    } else {
        //2D grid is required
        const unsigned int side = (unsigned int) ceil(sqrt((double)num_blocks));
        return dim3(side,side);
    }
}


namespace honei
{
    namespace cuda
    {
        // further optimised version: don't compute on zeros in offdiagonals
        __global__ void defect_q1_gpu(float * rhs,
                float* ll, float* ld, float* lu,
                float* dl, float * dd, float* du,
                float* ul, float* ud, float* uu,
                float * x, float * y, unsigned long n, unsigned long m)
        {
            extern __shared__ float  smv_cache[];

            unsigned long idx = blockDim.x*blockIdx.x+threadIdx.x;

            // runs from 0 to blockDim.x-1
            unsigned long lindex = threadIdx.x;

            float* Dcache = smv_cache;
            float* Lcache = smv_cache + blockDim.x + 2;
            float* Ucache = smv_cache + 2 * (blockDim.x + 2);

            // prefetch chunks from iteration vector
            //
            //
            // data needed for DD, DU, DL: each thread loads one element, the first and last one load the border cases
            // x_0 ... x_blockdim-1 into c_1...c_blockdim
            Dcache[lindex + 1] = x[idx];
            if (idx >= m) Lcache[lindex + 1] = x[idx - m];
            if (idx + m < n) Ucache[lindex + 1] = x[idx + m];
            if (lindex == 0)
            {
                // x_-1 in c_0
                if (blockDim.x * blockIdx.x - 1 < n) Dcache[0] = x[blockDim.x * blockIdx.x - 1];
                if (blockDim.x * blockIdx.x - m - 1 < n) Lcache[0] = x[blockDim.x * blockIdx.x - m - 1];
                if (blockDim.x * blockIdx.x + m - 1 < n) Ucache[0] = x[blockDim.x * blockIdx.x + m - 1];
            }
            if (lindex == blockDim.x - 1)
            {
                // x_blockdim in c_blockdim+1
                if (blockDim.x * (blockIdx.x + 1) < n) Dcache[blockDim.x + 1] = x[blockDim.x * (blockIdx.x + 1)];
                if (blockDim.x * (blockIdx.x + 1) - m < n) Lcache[blockDim.x + 1] = x[blockDim.x * (blockIdx.x + 1) - m];
                if (blockDim.x * (blockIdx.x + 1) + m  < n) Ucache[blockDim.x + 1] = x[blockDim.x * (blockIdx.x + 1) + m];
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
                y[idx] = rhs[idx] - ytemp1;
            }
        }

        __global__ void defect_smell_dv_gpu(float * rhs, float * y, unsigned long * Aj, float * Ax, float * x,
                unsigned long num_rows, unsigned long num_cols, unsigned long num_cols_per_row, unsigned long stride)
        {
            const unsigned long row = large_grid_thread_id();

            if(row >= num_rows){ return; }

            //float sum = y[row];
            float sum(float(0));

            Aj += row;
            Ax += row;

            for(unsigned long n = 0; n < num_cols_per_row; n++){
                const float A_ij = *Ax;

                if (A_ij != 0){
                    const unsigned long col = *Aj;
                    sum += A_ij * x[col];
                }

                Aj += stride;
                Ax += stride;
            }

            y[row] = rhs[row] - sum;
        }

        __global__ void defect_smell_dv_gpu(double * rhs, double * y, unsigned long * Aj, double * Ax, double * x,
                unsigned long num_rows, unsigned long num_cols, unsigned long num_cols_per_row, unsigned long stride)
        {
            const unsigned long row = large_grid_thread_id();

            if(row >= num_rows){ return; }

            //double sum = y[row];
            double sum(double(0));

            Aj += row;
            Ax += row;

            for(unsigned long n = 0; n < num_cols_per_row; n++){
                const double A_ij = *Ax;

                if (A_ij != 0){
                    const unsigned long col = *Aj;
                    sum += A_ij * x[col];
                }

                Aj += stride;
                Ax += stride;
            }

            y[row] = rhs[row] - sum;
        }
    }
}

extern "C" void cuda_defect_q1_float (void * rhs, void * ll, void * ld, void * lu,
        void * dl, void * dd, void *du,
        void * ul, void * ud, void *uu, void * x, void * y,
        unsigned long size, unsigned long blocksize, unsigned long m)
{
    dim3 grid;
    dim3 block;
    block.x = blocksize;
    grid.x = (unsigned)ceil(size/(double)(block.x));
    float * rhs_gpu((float *)rhs);
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


    honei::cuda::defect_q1_gpu<<<grid, block, 3 * (block.x + 2 ) * sizeof(float)>>>(rhs_gpu,
            ll_gpu, ld_gpu, lu_gpu, dl_gpu, dd_gpu, du_gpu, ul_gpu, ud_gpu, uu_gpu,
            x_gpu, y_gpu, size, m);

    CUDA_ERROR();
}

extern "C" void cuda_defect_smell_dv_float(void * rhs, void * result, void * Aj, void * Ax, void * b,
        unsigned long num_rows, unsigned long num_cols, unsigned long num_cols_per_row, unsigned long stride,
        unsigned long blocksize)
{
    const dim3 grid = make_large_grid_defect(num_rows, blocksize);

    float * rhs_gpu((float *)rhs);
    float * result_gpu((float *)result);
    float * b_gpu((float *)b);
    unsigned long * Aj_gpu((unsigned long *)Aj);
    float * Ax_gpu((float *)Ax);

    honei::cuda::defect_smell_dv_gpu<<<grid, blocksize>>>(rhs_gpu, result_gpu, Aj_gpu, Ax_gpu, b_gpu,
            num_rows, num_cols, num_cols_per_row, stride);

    CUDA_ERROR();
}

extern "C" void cuda_defect_smell_dv_double(void * rhs, void * result, void * Aj, void * Ax, void * b,
        unsigned long num_rows, unsigned long num_cols, unsigned long num_cols_per_row, unsigned long stride,
        unsigned long blocksize)
{
    const dim3 grid = make_large_grid_defect(num_rows, blocksize);

    double * rhs_gpu((double *)rhs);
    double * result_gpu((double *)result);
    double * b_gpu((double *)b);
    unsigned long * Aj_gpu((unsigned long *)Aj);
    double * Ax_gpu((double *)Ax);

    honei::cuda::defect_smell_dv_gpu<<<grid, blocksize>>>(rhs_gpu, result_gpu, Aj_gpu, Ax_gpu, b_gpu,
            num_rows, num_cols, num_cols_per_row, stride);

    CUDA_ERROR();
}
