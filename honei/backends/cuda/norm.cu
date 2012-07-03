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
#ifdef HONEI_CUBLAS
#include <cublas.h>
#endif

namespace honei
{
    namespace cuda
    {
        __global__ void norm_l2_gpu(float * x, float * tmp, unsigned long size, unsigned long blocksize)
        {
            // calculate how many elements each thread needs to calculate
            const unsigned long iter = size / (blockDim.x * gridDim.x);
            unsigned long pos = blockIdx.x* blocksize + threadIdx.x;

            // clear the output
            tmp[blockIdx.x * blocksize + threadIdx.x] = 0;

            for (unsigned long i = 0 ; i < iter ; ++i)
            {
                tmp[blockIdx.x * blocksize + threadIdx.x] += x[pos] * x[pos];
                pos += blockDim.x * gridDim.x;
            }

            // for the last iteration, check if the elements are still available
            if (pos < size)
            {
                tmp[blockIdx.x * blocksize + threadIdx.x] += x[pos] * x[pos];
            }
        }

#ifdef HONEI_CUDA_DOUBLE
        __global__ void norm_l2_gpu(double * x, double * tmp, unsigned long size, unsigned long blocksize)
        {
            // calculate how many elements each thread needs to calculate
            const unsigned long iter = size / (blockDim.x * gridDim.x);
            unsigned long pos = blockIdx.x* blocksize + threadIdx.x;

            // clear the output
            tmp[blockIdx.x * blocksize + threadIdx.x] = 0;

            for (unsigned long i = 0 ; i < iter ; ++i)
            {
                tmp[blockIdx.x * blocksize + threadIdx.x] += x[pos] * x[pos];
                pos += blockDim.x * gridDim.x;
            }

            // for the last iteration, check if the elements are still available
            if (pos < size)
            {
                tmp[blockIdx.x * blocksize + threadIdx.x] += x[pos] * x[pos];
            }
        }
#endif
    }
}

extern "C" float cuda_norm_l2_one_float(const void * x, unsigned long size, unsigned long blocksize,
        unsigned long gridsize)
{
#ifdef HONEI_CUBLAS
    float * x_gpu((float *)x);
    cublasInit();
    float result = cublasSnrm2(size, x_gpu, 1);
    cublasShutdown();
    return result;
#else
    float result(0.);

    dim3 grid(gridsize);
    dim3 block(blocksize);
    float * x_gpu((float* )x);
    float * tmp_cpu(0);
    float * tmp_gpu(0);

    cudaMalloc((void**)&tmp_gpu, gridsize * blocksize * sizeof(float));
    //cudaMallocHost((void**)&tmp_cpu, gridsize * blocksize * sizeof(float));
    tmp_cpu = new float[gridsize * blocksize];

    honei::cuda::norm_l2_gpu<<<grid, block>>>(x_gpu, tmp_gpu, size, blocksize);

    cudaMemcpy(tmp_cpu, tmp_gpu, blocksize * gridsize * sizeof(float), cudaMemcpyDeviceToHost);
    for (unsigned long i(0) ; i < blocksize * gridsize ; ++i)
    {
        result += tmp_cpu[i];
    }

    cudaFree(tmp_gpu);
    //cudaFreeHost(tmp_cpu);
    delete[] tmp_cpu;

    CUDA_ERROR();
    return (result);
#endif
}

#ifdef HONEI_CUDA_DOUBLE
extern "C" double cuda_norm_l2_one_double(const void * x, unsigned long size, unsigned long blocksize,
        unsigned long gridsize)
{
#ifdef HONEI_CUBLAS
    double * x_gpu((double *)x);
    cublasInit();
    double result = cublasDnrm2(size, x_gpu, 1);
    cublasShutdown();
    return result;
#else
    double result(0.);

    dim3 grid(gridsize);
    dim3 block(blocksize);
    double * x_gpu((double* )x);
    double * tmp_cpu(0);
    double * tmp_gpu(0);

    cudaMalloc((void**)&tmp_gpu, gridsize * blocksize * sizeof(double));
    //cudaMallocHost((void**)&tmp_cpu, gridsize * blocksize * sizeof(double));
    tmp_cpu = new double[gridsize * blocksize];

    honei::cuda::norm_l2_gpu<<<grid, block>>>(x_gpu, tmp_gpu, size, blocksize);

    cudaMemcpy(tmp_cpu, tmp_gpu, blocksize * gridsize * sizeof(double), cudaMemcpyDeviceToHost);
    for (unsigned long i(0) ; i < blocksize * gridsize ; ++i)
    {
        result += tmp_cpu[i];
    }

    cudaFree(tmp_gpu);
    //cudaFreeHost(tmp_cpu);
    delete[] tmp_cpu;

    CUDA_ERROR();
    return (result);
#endif
}
#endif
