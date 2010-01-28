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
        __global__ void extraction_grid_dry_gpu(
                float * f_0, float * f_1, float * f_2,
                float * f_3, float * f_4, float * f_5,
                float * f_6, float * f_7, float * f_8,
                float * h, float * u, float * v,
                float * distribution_x_data, float * distribution_y_data, float epsilon,
                unsigned long offset, unsigned long size)
        {
            extern __shared__ float  distribution_cache[];
            float* distribution_x = distribution_cache;
            float* distribution_y = distribution_cache + 9;

            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;

            if (idx % blockDim.x < 9)
            {
                distribution_x[idx % blockDim.x] = distribution_x_data[idx % blockDim.x];
                distribution_y[idx % blockDim.x] = distribution_y_data[idx % blockDim.x];
            }
            __syncthreads();

            if (idx >= offset && idx < size)
            {
                unsigned long i(idx);

                h[i] = f_0[i] +
                    f_1[i] +
                    f_2[i] +
                    f_3[i] +
                    f_4[i] +
                    f_5[i] +
                    f_6[i] +
                    f_7[i] +
                    f_8[i];

                float lax_upper(epsilon);
                float lax_lower(-lax_upper);

                if(h[i] < lax_lower || h[i] > lax_upper)
                {
                    u[i] = (distribution_x[0] * f_0[i] +
                            distribution_x[1] * f_1[i] +
                            distribution_x[2] * f_2[i] +
                            distribution_x[3] * f_3[i] +
                            distribution_x[4] * f_4[i] +
                            distribution_x[5] * f_5[i] +
                            distribution_x[6] * f_6[i] +
                            distribution_x[7] * f_7[i] +
                            distribution_x[8] * f_8[i]) / h[i];

                    v[i] = (distribution_y[0] * f_0[i] +
                            distribution_y[1] * f_1[i] +
                            distribution_y[2] * f_2[i] +
                            distribution_y[3] * f_3[i] +
                            distribution_y[4] * f_4[i] +
                            distribution_y[5] * f_5[i] +
                            distribution_y[6] * f_6[i] +
                            distribution_y[7] * f_7[i] +
                            distribution_y[8] * f_8[i]) / h[i];
                }
                else
                {
                    h[i] = 0;
                    u[i] = 0;
                    v[i] = 0;
                }
                h[i] = max(float(0), min(float(1), h[i]));
            }
        }

        __global__ void extraction_grid_wet_gpu(
                float * f_0, float * f_1, float * f_2,
                float * f_3, float * f_4, float * f_5,
                float * f_6, float * f_7, float * f_8,
                float * h, float * u, float * v,
                float * distribution_x_data, float * distribution_y_data, float epsilon,
                unsigned long offset, unsigned long size)
        {
            extern __shared__ float  distribution_cache[];
            float* distribution_x = distribution_cache;
            float* distribution_y = distribution_cache + 9;

            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;

            if (idx % blockDim.x < 9)
            {
                distribution_x[idx % blockDim.x] = distribution_x_data[idx % blockDim.x];
                distribution_y[idx % blockDim.x] = distribution_y_data[idx % blockDim.x];
            }
            __syncthreads();

            if (idx >= offset && idx < size)
            {
                unsigned long i(idx);

                h[i] = f_0[i] +
                    f_1[i] +
                    f_2[i] +
                    f_3[i] +
                    f_4[i] +
                    f_5[i] +
                    f_6[i] +
                    f_7[i] +
                    f_8[i];

                u[i] = (distribution_x[0] * f_0[i] +
                        distribution_x[1] * f_1[i] +
                        distribution_x[2] * f_2[i] +
                        distribution_x[3] * f_3[i] +
                        distribution_x[4] * f_4[i] +
                        distribution_x[5] * f_5[i] +
                        distribution_x[6] * f_6[i] +
                        distribution_x[7] * f_7[i] +
                        distribution_x[8] * f_8[i]) / h[i];

                v[i] = (distribution_y[0] * f_0[i] +
                        distribution_y[1] * f_1[i] +
                        distribution_y[2] * f_2[i] +
                        distribution_y[3] * f_3[i] +
                        distribution_y[4] * f_4[i] +
                        distribution_y[5] * f_5[i] +
                        distribution_y[6] * f_6[i] +
                        distribution_y[7] * f_7[i] +
                        distribution_y[8] * f_8[i]) / h[i];
            }
        }
    }
}

extern "C" void cuda_extraction_grid_dry_float(
        unsigned long start, unsigned long end,
        void * f_0, void * f_1, void * f_2,
        void * f_3, void * f_4, void * f_5,
        void * f_6, void * f_7, void * f_8,
        void * h, void * u, void * v,
        void * distribution_x, void * distribution_y, float epsilon,
        unsigned long blocksize)
{
    unsigned long size(end);
    dim3 grid;
    dim3 block;
    block.x = blocksize;
    grid.x = (unsigned)ceil(sqrt(size/(double)block.x));
    grid.y = grid.x;

    float * f_0_gpu((float *)f_0);
    float * f_1_gpu((float *)f_1);
    float * f_2_gpu((float *)f_2);
    float * f_3_gpu((float *)f_3);
    float * f_4_gpu((float *)f_4);
    float * f_5_gpu((float *)f_5);
    float * f_6_gpu((float *)f_6);
    float * f_7_gpu((float *)f_7);
    float * f_8_gpu((float *)f_8);

    float * h_gpu((float *)h);
    float * u_gpu((float *)u);
    float * v_gpu((float *)v);
    float * distribution_x_gpu((float *)distribution_x);
    float * distribution_y_gpu((float *)distribution_y);

    honei::cuda::extraction_grid_dry_gpu<<<grid, block, 18 * sizeof(float)>>>(
            f_0_gpu, f_1_gpu, f_2_gpu, f_3_gpu, f_4_gpu,
            f_5_gpu, f_6_gpu, f_7_gpu, f_8_gpu,
            h_gpu, u_gpu, v_gpu,
            distribution_x_gpu, distribution_y_gpu, epsilon,
            start, size);

    CUDA_ERROR();
}

extern "C" void cuda_extraction_grid_wet_float(
        unsigned long start, unsigned long end,
        void * f_0, void * f_1, void * f_2,
        void * f_3, void * f_4, void * f_5,
        void * f_6, void * f_7, void * f_8,
        void * h, void * u, void * v,
        void * distribution_x, void * distribution_y, float epsilon,
        unsigned long blocksize)
{
    unsigned long size(end);
    dim3 grid;
    dim3 block;
    block.x = blocksize;
    grid.x = (unsigned)ceil(sqrt(size/(double)block.x));
    grid.y = grid.x;

    float * f_0_gpu((float *)f_0);
    float * f_1_gpu((float *)f_1);
    float * f_2_gpu((float *)f_2);
    float * f_3_gpu((float *)f_3);
    float * f_4_gpu((float *)f_4);
    float * f_5_gpu((float *)f_5);
    float * f_6_gpu((float *)f_6);
    float * f_7_gpu((float *)f_7);
    float * f_8_gpu((float *)f_8);

    float * h_gpu((float *)h);
    float * u_gpu((float *)u);
    float * v_gpu((float *)v);
    float * distribution_x_gpu((float *)distribution_x);
    float * distribution_y_gpu((float *)distribution_y);

    honei::cuda::extraction_grid_wet_gpu<<<grid, block, 18 * sizeof(float)>>>(
            f_0_gpu, f_1_gpu, f_2_gpu, f_3_gpu, f_4_gpu,
            f_5_gpu, f_6_gpu, f_7_gpu, f_8_gpu,
            h_gpu, u_gpu, v_gpu,
            distribution_x_gpu, distribution_y_gpu, epsilon,
            start, size);

    CUDA_ERROR();
}
