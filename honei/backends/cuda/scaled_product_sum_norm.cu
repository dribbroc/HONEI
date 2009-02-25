/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2009 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 * Copyright (c) 2009 Markus Geveler <apryde@gmx.de>
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
        __global__ void scaled_product_sum_norm_tut_gpu(float a, float * y, float b, float * A_x, float * tmp, unsigned long size, unsigned long blocksize)
        {
            // calculate how many elements each thread needs to calculate
            const unsigned long iter = size / (blockDim.x * gridDim.x);
            unsigned long pos = blockIdx.x* blocksize + threadIdx.x;

            // clear the output
            tmp[blockIdx.x * blocksize + threadIdx.x] = 0;

            for (unsigned long i = 0 ; i < iter ; ++i)
            {
                tmp[blockIdx.x * blocksize + threadIdx.x] += (A_x[pos] * b + y[pos] * a) * (A_x[pos] * b + y[pos] * a);
                pos += blockDim.x * gridDim.x;
            }

            // for the last iteration, check if the elements are still available
            if (pos < size)
            {
                tmp[blockIdx.x * blocksize + threadIdx.x] += (A_x[pos] * b + y[pos] * a) * (A_x[pos] * b + y[pos] * a);
            }
        }
    }
}

extern "C" float cuda_scaled_product_sum_norm_float_tut(unsigned long size, float a, void * y, float b, void * A_x, unsigned long blocksize, unsigned long gridsize)
{
    float result(0.);

    dim3 grid(gridsize);
    dim3 block(blocksize);
    float * A_x_gpu((float* )A_x);
    float * y_gpu((float* )y);
    float * tmp_cpu(0);
    float * tmp_gpu(0);

    cudaMalloc((void**)&tmp_gpu, gridsize * blocksize * sizeof(float));
    cudaMallocHost((void**)&tmp_cpu, gridsize * blocksize * sizeof(float));

    honei::cuda::scaled_product_sum_norm_tut_gpu<<<grid, block>>>(a, y_gpu, b, A_x_gpu, tmp_gpu, size, blocksize);

    cudaMemcpy(tmp_cpu, tmp_gpu, blocksize * gridsize * sizeof(float), cudaMemcpyDeviceToHost);
    for (unsigned long i(0) ; i < blocksize * gridsize ; ++i)
    {
        result += tmp_cpu[i];
    }

    cudaFree(tmp_gpu);
    cudaFreeHost(tmp_cpu);

    CUDA_ERROR();
    return result;
}
