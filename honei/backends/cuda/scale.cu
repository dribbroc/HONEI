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
        __global__ void scale_gpu(float a, float * x, unsigned long size)
        {
            int idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            if (idx < size)
            {
                x[idx] = x[idx] * a;
            }
        }
    }
}

extern "C" void cuda_scale_one_float(float a, float * x, unsigned long size, unsigned long blocksize)
{
    dim3 grid;
    dim3 block;
    block.x = blocksize;
    grid.x = ceil(sqrt(size/(double)block.x));
    grid.y = grid.x;
    float * x_gpu(0);

    cudaMalloc((void**)&x_gpu, size * sizeof(float));

    cudaMemcpy(x_gpu, x, size * sizeof(float), cudaMemcpyHostToDevice);

    honei::cuda::scale_gpu<<<grid, block>>>(a, x_gpu, size);

    cudaMemcpy(x, x_gpu, size * sizeof(float), cudaMemcpyDeviceToHost);
    cudaFree(x_gpu);

    CUDA_ERROR();
}
