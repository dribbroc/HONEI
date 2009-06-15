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
        __global__ void element_product_gpu(float * r, float * x, float * y, unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            if (idx < size)
            {
                r[idx] = x[idx] * y[idx];
            }
        }

#ifdef HONEI_CUDA_DOUBLE
        __global__ void element_product_gpu(double * r, double * x, double * y, unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            if (idx < size)
            {
                r[idx] = x[idx] * y[idx];
            }
        }
#endif
    }
}

extern "C" void cuda_element_product_three_float(void * r, void * x, void * y, unsigned long size, unsigned long blocksize)
{
    dim3 grid;
    dim3 block;
    block.x = blocksize;
    grid.x = (unsigned)ceil(sqrt(size/(double)block.x));
    grid.y = grid.x;
    float * r_gpu((float *)r);
    float * x_gpu((float *)x);
    float * y_gpu((float *)y);

    honei::cuda::element_product_gpu<<<grid, block>>>(r_gpu, x_gpu, y_gpu, size);

    CUDA_ERROR();
}

#ifdef HONEI_CUDA_DOUBLE
extern "C" void cuda_element_product_three_double(void * r, void * x, void * y, unsigned long size, unsigned long blocksize)
{
    dim3 grid;
    dim3 block;
    block.x = blocksize;
    grid.x = (unsigned)ceil(sqrt(size/(double)block.x));
    grid.y = grid.x;
    double * r_gpu((double *)r);
    double * x_gpu((double *)x);
    double * y_gpu((double *)y);

    honei::cuda::element_product_gpu<<<grid, block>>>(r_gpu, x_gpu, y_gpu, size);

    CUDA_ERROR();
}
#endif
