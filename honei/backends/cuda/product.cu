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
            if (idx < n) Dcache[lindex + 1] = x[idx];
            if (idx >= m && idx - m < n) Lcache[lindex + 1] = x[idx - m];
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
            if (idx < n) Dcache[lindex + 1] = x[idx];
            if (idx >= m && idx - m < n) Lcache[lindex + 1] = x[idx - m];
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

        __global__ void product_smell_dv_gpu(float * x, float * y, const unsigned long * Aj, const float * Ax, const unsigned long * Arl,
                unsigned long row_start, unsigned long row_end, unsigned long num_cols_per_row, unsigned long stride, unsigned long threads)
        {
            extern __shared__ float  shared_ell_float[];

            const unsigned long T = threads;
            const unsigned long idx = blockDim.x*blockIdx.x+threadIdx.x + (row_start * T);
            const unsigned long idb = threadIdx.x;
            const unsigned long idp = idb % T;
            const unsigned long row = idx / T;

            if(row >= row_end){ return; }
            shared_ell_float[idb] = 0;
            float sum = float(0);

            const unsigned long max = Arl[row];
            Ax += (row*T)+idp;
            Aj += (row*T)+idp;
            for(unsigned long k = 0; k < max ; ++k)
            {
                //sum += value * x[col];
                sum += *Ax * tex1Dfetch(tex_x_float_product, *Aj);
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
                        y[row - row_start] = shared_ell_float[idb] + shared_ell_float[idb + 1];
                    break;
                case 1:
                    y[row - row_start] = shared_ell_float[idb];
                    break;
                default:
                    break;
            }
        }

#ifdef HONEI_CUDA_DOUBLE
        __global__ void product_smell_dv_gpu(double * x, double * y, const unsigned long * Aj, const double * Ax, const unsigned long * Arl,
                unsigned long row_start, unsigned long row_end, unsigned long num_cols_per_row, unsigned long stride, unsigned long threads)
        {
            extern __shared__ double  shared_ell_double[];

            const unsigned long T = threads;
            const unsigned long idx = blockDim.x*blockIdx.x+threadIdx.x + (row_start * T);
            const unsigned long idb = threadIdx.x;
            const unsigned long idp = idb % T;
            const unsigned long row = idx / T;

            if(row >= row_end){ return; }
            shared_ell_double[idb] = 0;
            double sum = double(0);

            const unsigned long max = Arl[row];
            Ax += (row*T)+idp;
            Aj += (row*T)+idp;
            for(unsigned long k = 0; k < max ; ++k)
            {
                //sum += value * x[col];
                int2 v = tex1Dfetch(tex_x_double_product, *Aj);
                sum += *Ax * __hiloint2double(v.y, v.x);
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
                        y[row - row_start] = shared_ell_double[idb] + shared_ell_double[idb + 1];
                    break;
                case 1:
                    y[row - row_start] = shared_ell_double[idb];
                    break;
                default:
                    break;
            }
        }
#endif

        __global__ void product_csr_dv_gpu(float * x, float * y, const unsigned long * Aj, const float * Ax, const unsigned long * Ar,
                unsigned long row_start, unsigned long row_end, unsigned long atomicsize)
        {
            unsigned long idx = blockDim.x*blockIdx.x+threadIdx.x;
            unsigned long row = idx + row_start;

            if(row >= row_end){ return; }

            float sum(0);
            const unsigned long end(Ar[row+1]);
            for (unsigned long i(Ar[row]) ; i < end ; ++i)
            {
                for (unsigned long blocki(0) ; blocki < atomicsize ; ++blocki)
                    sum += Ax[(i * atomicsize) + blocki] * tex1Dfetch(tex_x_float_product, Aj[i] + blocki);
            }
            y[row] = sum;
        }

#ifdef HONEI_CUDA_DOUBLE
        __global__ void product_csr_dv_gpu(double * x, double * y, const unsigned long * Aj, const double * Ax, const unsigned long * Ar,
                unsigned long row_start, unsigned long row_end, unsigned long atomicsize)
        {
            unsigned long idx = blockDim.x*blockIdx.x+threadIdx.x;
            unsigned long row = idx + row_start;

            if(row >= row_end){ return; }

            double sum(0);
            const unsigned long end(Ar[row+1]);
            for (unsigned long i(Ar[row]) ; i < end ; ++i)
            {
                for (unsigned long blocki(0) ; blocki < atomicsize ; ++blocki)
                {
                    int2 v = tex1Dfetch(tex_x_double_product, Aj[i] + blocki);
                    sum += Ax[(i * atomicsize) + blocki] * __hiloint2double(v.y, v.x);
                }
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

extern "C" void cuda_product_smell_dv_float(void * x, void * y, void * Aj, void * Ax, void * Arl,
        unsigned long row_start, unsigned long row_end, unsigned long num_cols_per_row,
        unsigned long stride, unsigned long blocksize, unsigned long threads, cudaStream_t stream)
{
    dim3 grid;
    dim3 block;
    block.x = blocksize;
    grid.x = (unsigned)ceil((threads*(row_end - row_start))/(double)(block.x));

    float * x_gpu((float *)x);
    float * y_gpu((float *)y);
    unsigned long * Aj_gpu((unsigned long *)Aj);
    float * Ax_gpu((float *)Ax);
    unsigned long * Arl_gpu((unsigned long *)Arl);

    cudaBindTexture(NULL, tex_x_float_product, x_gpu);
    honei::cuda::product_smell_dv_gpu<<<grid, block, block.x * sizeof(float), stream>>>(x_gpu, y_gpu, Aj_gpu, Ax_gpu, Arl_gpu,
            row_start, row_end, num_cols_per_row, stride, threads);
    cudaUnbindTexture(tex_x_float_product);

    CUDA_ERROR();
}

#ifdef HONEI_CUDA_DOUBLE
extern "C" void cuda_product_smell_dv_double(void * x, void * y, void * Aj, void * Ax, void * Arl,
        unsigned long row_start, unsigned long row_end, unsigned long num_cols_per_row,
        unsigned long stride, unsigned long blocksize, unsigned long threads, cudaStream_t stream)
{
    dim3 grid;
    dim3 block;
    block.x = blocksize;
    grid.x = (unsigned)ceil((threads*(row_end - row_start))/(double)(block.x));

    double * x_gpu((double *)x);
    double * y_gpu((double *)y);
    unsigned long * Aj_gpu((unsigned long *)Aj);
    double * Ax_gpu((double *)Ax);
    unsigned long * Arl_gpu((unsigned long *)Arl);

    cudaBindTexture(NULL, tex_x_double_product, x_gpu);
    honei::cuda::product_smell_dv_gpu<<<grid, block, block.x * sizeof(double), stream>>>(x_gpu, y_gpu, Aj_gpu, Ax_gpu, Arl_gpu,
            row_start, row_end, num_cols_per_row, stride, threads);
    cudaUnbindTexture(tex_x_double_product);

    CUDA_ERROR();
}
#endif

extern "C" void cuda_product_csr_dv_float(void * x, void * y, void * Aj, void * Ax, void * Ar,
        unsigned long row_start, unsigned long row_end, unsigned long atomicsize, unsigned long blocksize, cudaStream_t stream)
{
    dim3 grid;
    dim3 block;
    block.x = blocksize;
    grid.x = (unsigned)ceil((row_end - row_start)/(double)(block.x));

    float * x_gpu((float *)x);
    float * y_gpu((float *)y);
    unsigned long * Aj_gpu((unsigned long *)Aj);
    float * Ax_gpu((float *)Ax);
    unsigned long * Ar_gpu((unsigned long *)Ar);

    cudaBindTexture(NULL, tex_x_float_product, x_gpu);
    honei::cuda::product_csr_dv_gpu<<<grid, block, block.x * sizeof(float), stream>>>(x_gpu, y_gpu, Aj_gpu, Ax_gpu, Ar_gpu,
            row_start, row_end, atomicsize);
    cudaUnbindTexture(tex_x_float_product);

    CUDA_ERROR();
}

#ifdef HONEI_CUDA_DOUBLE
extern "C" void cuda_product_csr_dv_double(void * x, void * y, void * Aj, void * Ax, void * Ar,
        unsigned long row_start, unsigned long row_end, unsigned long atomicsize, unsigned long blocksize, cudaStream_t stream)
{
    dim3 grid;
    dim3 block;
    block.x = blocksize;
    grid.x = (unsigned)ceil((row_end - row_start)/(double)(block.x));

    double * x_gpu((double *)x);
    double * y_gpu((double *)y);
    unsigned long * Aj_gpu((unsigned long *)Aj);
    double * Ax_gpu((double *)Ax);
    unsigned long * Ar_gpu((unsigned long *)Ar);

    cudaBindTexture(NULL, tex_x_double_product, x_gpu);
    honei::cuda::product_csr_dv_gpu<<<grid, block, block.x * sizeof(double), stream>>>(x_gpu, y_gpu, Aj_gpu, Ax_gpu, Ar_gpu,
            row_start, row_end, atomicsize);
    cudaUnbindTexture(tex_x_double_product);

    CUDA_ERROR();
}
#endif
