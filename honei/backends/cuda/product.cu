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

dim3 make_large_grid_product(const unsigned int num_threads, const unsigned int blocksize){
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

texture<float,1> tex_x_float_product;
texture<int2,1>  tex_x_double_product;

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
                float * x, float * y, unsigned long n, unsigned long m)
        {
            extern __shared__ float  smvf_cache[];

            unsigned long idx = blockDim.x*blockIdx.x+threadIdx.x;

            // runs from 0 to blockDim.x-1
            unsigned long lindex = threadIdx.x;

            float* Dcache = smvf_cache;
            float* Lcache = smvf_cache + blockDim.x + 2;
            float* Ucache = smvf_cache + 2 * (blockDim.x + 2);

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
                y[idx] = ytemp1;
            }
        }

        #ifdef HONEI_CUDA_DOUBLE
        __global__ void product_bmdv_q1_gpu(
                double* ll, double* ld, double* lu,
                double* dl, double * dd, double* du,
                double* ul, double* ud, double* uu,
                double * x, double * y, unsigned long n, unsigned long m)
        {
            extern __shared__ double  smvd_cache[];

            unsigned long idx = blockDim.x*blockIdx.x+threadIdx.x;

            // runs from 0 to blockDim.x-1
            unsigned long lindex = threadIdx.x;

            double* Dcache = smvd_cache;
            double* Lcache = smvd_cache + blockDim.x + 2;
            double* Ucache = smvd_cache + 2 * (blockDim.x + 2);

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

        __global__ void product_smell_dv_gpu(float * x, float * y, unsigned long * Aj, float * Ax, unsigned long * Arl,
                unsigned long num_rows, unsigned long num_cols, unsigned long num_cols_per_row, unsigned long stride)
        {
            const unsigned long row = large_grid_thread_id();

            if(row >= num_rows){ return; }

            //float sum = y[row];
            float sum = float(0);

            Aj += row;
            Ax += row;

            const unsigned long max = Arl[row];
            for(unsigned long n = 0; n < max ; n++)
            {
                const float A_ij = *Ax;

                //if (A_ij != 0)
                {
                    const unsigned long col = *Aj;
                    //sum += A_ij * x[col];
                    sum += A_ij * tex1Dfetch(tex_x_float_product, col);
                }

                Aj += stride;
                Ax += stride;
            }

            y[row] = sum;
        }

#ifdef HONEI_CUDA_DOUBLE
        __global__ void product_smell_dv_gpu(double * x, double * y, unsigned long * Aj, double * Ax, unsigned long * Arl,
                unsigned long num_rows, unsigned long num_cols, unsigned long num_cols_per_row, unsigned long stride)
        {
            const unsigned long row = large_grid_thread_id();

            if(row >= num_rows){ return; }

            //double sum = y[row];
            double sum = double(0);

            Aj += row;
            Ax += row;

            const unsigned long max = Arl[row];
            for(unsigned long n = 0; n < max ; n++)
            {
                const double A_ij = *Ax;

                //if (A_ij != 0)
                {
                    const unsigned long col = *Aj;
                    //sum += A_ij * x[col];
                    int2 v = tex1Dfetch(tex_x_double_product, col);
                    sum += A_ij * __hiloint2double(v.y, v.x);
                }

                Aj += stride;
                Ax += stride;
            }

            y[row] = sum;
        }
#endif
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
    grid.x = (unsigned)ceil(size/(double)(block.x));
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

#ifdef HONEI_CUDA_DOUBLE
extern "C" void cuda_product_bmdv_q1_double (void * ll, void * ld, void * lu,
        void * dl, void * dd, void *du,
        void * ul, void * ud, void *uu, void * x, void * y,
        unsigned long size, unsigned long blocksize, unsigned long m)
{
    dim3 grid;
    dim3 block;
    block.x = blocksize;
    grid.x = (unsigned)ceil(size/(double)(block.x));
    double * x_gpu((double *)x);
    double * y_gpu((double *)y);
    double * ll_gpu((double *)ll);
    double * ld_gpu((double *)ld);
    double * lu_gpu((double *)lu);
    double * dl_gpu((double *)dl);
    double * dd_gpu((double *)dd);
    double * du_gpu((double *)du);
    double * ul_gpu((double *)ul);
    double * ud_gpu((double *)ud);
    double * uu_gpu((double *)uu);


    honei::cuda::product_bmdv_q1_gpu<<<grid, block, 3 * (block.x + 2 ) * sizeof(double)>>>(ll_gpu, ld_gpu, lu_gpu, dl_gpu, dd_gpu, du_gpu, ul_gpu, ud_gpu, uu_gpu, x_gpu, y_gpu, size, m);

    CUDA_ERROR();
}
#endif

extern "C" void cuda_product_smell_dv_float(void * x, void * y, void * Aj, void * Ax, void *Arl,
        unsigned long num_rows, unsigned long num_cols, unsigned long num_cols_per_row,
        unsigned long stride, unsigned long blocksize)
{
    const dim3 grid = make_large_grid_product(num_rows, blocksize);

    float * x_gpu((float *)x);
    float * y_gpu((float *)y);
    unsigned long * Aj_gpu((unsigned long *)Aj);
    float * Ax_gpu((float *)Ax);
    unsigned long * Arl_gpu((unsigned long *)Arl);

    cudaBindTexture(NULL, tex_x_float_product, x_gpu);
    honei::cuda::product_smell_dv_gpu<<<grid, blocksize>>>(x_gpu, y_gpu, Aj_gpu, Ax_gpu, Arl_gpu,
            num_rows, num_cols, num_cols_per_row, stride);
    cudaUnbindTexture(tex_x_float_product);

    CUDA_ERROR();
}

#ifdef HONEI_CUDA_DOUBLE
extern "C" void cuda_product_smell_dv_double(void * x, void * y, void * Aj, void * Ax, void * Arl,
        unsigned long num_rows, unsigned long num_cols, unsigned long num_cols_per_row,
        unsigned long stride, unsigned long blocksize)
{
    const dim3 grid = make_large_grid_product(num_rows, blocksize);

    double * x_gpu((double *)x);
    double * y_gpu((double *)y);
    unsigned long * Aj_gpu((unsigned long *)Aj);
    double * Ax_gpu((double *)Ax);
    unsigned long * Arl_gpu((unsigned long *)Arl);

    cudaBindTexture(NULL, tex_x_double_product, x_gpu);
    honei::cuda::product_smell_dv_gpu<<<grid, blocksize>>>(x_gpu, y_gpu, Aj_gpu, Ax_gpu, Arl_gpu,
            num_rows, num_cols, num_cols_per_row, stride);
    cudaUnbindTexture(tex_x_double_product);

    CUDA_ERROR();
}
#endif
