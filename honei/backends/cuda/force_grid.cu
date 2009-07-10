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
        __global__ void force_grid_gpu_x(
                unsigned long * dir_1,
                float * f_temp_1,
                float * h, float * b,
                float distribution_x,
                float g, float d_x, float d_y, float d_t,
                unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            float force_multiplier(d_t / (6 * d_x * d_x / (d_t * d_t)));
            float gravity_multiplier(-g);
            float force_times_gravity(force_multiplier * gravity_multiplier);
            float dx_by_two(float(2) * d_x);

            if (idx < size)
            {
                unsigned long i(idx);
                if (dir_1[i] < size)
                {
                    f_temp_1[i] += force_times_gravity * distribution_x * (h[i] + h[dir_1[i]]) * (b[dir_1[i]] - b[i]) / dx_by_two;
                }
            }
        }
        __global__ void force_grid_gpu_xy(
                unsigned long * dir_2,
                unsigned long * dir_1,
                unsigned long * dir_3,
                float * f_temp_2,
                float * h, float * b,
                float distribution_x, float distribution_y,
                float g, float d_x, float d_y, float d_t,
                unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            float force_multiplier(d_t / (6 * d_x * d_x / (d_t * d_t)));
            float gravity_multiplier(-g);
            float force_times_gravity(force_multiplier * gravity_multiplier);
            float dx_by_two(float(2) * d_x);
            float dy_by_two(float(2) * d_y);

            if (idx < size)
            {
                unsigned long i(idx);
                if (dir_2[i] < size && dir_3[i] < size && dir_1[i] < size)
                {
                    f_temp_2[i] += (force_times_gravity * distribution_x * (h[i] + h[dir_2[i]]) * (b[dir_1[i]] - b[i]) / dx_by_two) +
                                   (force_times_gravity * distribution_y * (h[i] + h[dir_2[i]]) * (b[dir_3[i]] - b[i]) / dy_by_two);
                }
            }
        }
        __global__ void force_grid_gpu_y(
                unsigned long * dir_3,
                float * f_temp_3,
                float * h, float * b,
                float distribution_y,
                float g, float d_x, float d_y, float d_t,
                unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            float force_multiplier(d_t / (6 * d_x * d_x / (d_t * d_t)));
            float gravity_multiplier(-g);
            float force_times_gravity(force_multiplier * gravity_multiplier);
            float dy_by_two(float(2) * d_y);

            if (idx < size)
            {
                unsigned long i(idx);
                if (dir_3[i] < size)
                {
                    f_temp_3[i] += force_times_gravity * distribution_y * (h[i] + h[dir_3[i]]) * (b[dir_3[i]] - b[i]) / dy_by_two;
                }
            }
        }
    }
}

extern "C" void cuda_force_grid_float(
        void * dir_1, void * dir_2, void * dir_3, void * dir_4,
        void * dir_5, void * dir_6, void * dir_7, void * dir_8,
        void * h, void * b,
        void * f_temp_1, void * f_temp_2,
        void * f_temp_3, void * f_temp_4, void * f_temp_5,
        void * f_temp_6, void * f_temp_7, void * f_temp_8,
        float distribution_x_1, float distribution_y_1,
        float distribution_x_2, float distribution_y_2,
        float distribution_x_3, float distribution_y_3,
        float distribution_x_4, float distribution_y_4,
        float distribution_x_5, float distribution_y_5,
        float distribution_x_6, float distribution_y_6,
        float distribution_x_7, float distribution_y_7,
        float distribution_x_8, float distribution_y_8,
        float g, float d_x, float d_y, float d_t,
        unsigned long size,
        unsigned long blocksize)
{
    dim3 grid;
    dim3 block;
    block.x = blocksize;
    grid.x = (unsigned)ceil(sqrt(size/(double)block.x));
    grid.y = grid.x;

    unsigned long * dir_1_gpu((unsigned long *)dir_1);
    unsigned long * dir_2_gpu((unsigned long *)dir_2);
    unsigned long * dir_3_gpu((unsigned long *)dir_3);
    unsigned long * dir_4_gpu((unsigned long *)dir_4);
    unsigned long * dir_5_gpu((unsigned long *)dir_5);
    unsigned long * dir_6_gpu((unsigned long *)dir_6);
    unsigned long * dir_7_gpu((unsigned long *)dir_7);
    unsigned long * dir_8_gpu((unsigned long *)dir_8);

    float * h_gpu((float *)h);
    float * b_gpu((float *)b);

    float * f_temp_1_gpu((float *)f_temp_1);
    float * f_temp_2_gpu((float *)f_temp_2);
    float * f_temp_3_gpu((float *)f_temp_3);
    float * f_temp_4_gpu((float *)f_temp_4);
    float * f_temp_5_gpu((float *)f_temp_5);
    float * f_temp_6_gpu((float *)f_temp_6);
    float * f_temp_7_gpu((float *)f_temp_7);
    float * f_temp_8_gpu((float *)f_temp_8);

    honei::cuda::force_grid_gpu_x<<<grid, block>>>(
            dir_1_gpu, f_temp_1_gpu,
            h_gpu, b_gpu,
            distribution_x_1,
            g, d_x, d_y, d_t,
            size);
    honei::cuda::force_grid_gpu_xy<<<grid, block>>>(
            dir_2_gpu, dir_1_gpu, dir_3_gpu, f_temp_2_gpu,
            h_gpu, b_gpu,
            distribution_x_2, distribution_y_2,
            g, d_x, d_y, d_t,
            size);
    honei::cuda::force_grid_gpu_y<<<grid, block>>>(
            dir_3_gpu, f_temp_3_gpu,
            h_gpu, b_gpu,
            distribution_y_3,
            g, d_x, d_y, d_t,
            size);
    honei::cuda::force_grid_gpu_xy<<<grid, block>>>(
            dir_4_gpu, dir_1_gpu, dir_3_gpu, f_temp_4_gpu,
            h_gpu, b_gpu,
            distribution_x_4, distribution_y_4,
            g, d_x, d_y, d_t,
            size);
    honei::cuda::force_grid_gpu_x<<<grid, block>>>(
            dir_5_gpu, f_temp_5_gpu,
            h_gpu, b_gpu,
            distribution_x_5,
            g, d_x, d_y, d_t,
            size);
    honei::cuda::force_grid_gpu_xy<<<grid, block>>>(
            dir_6_gpu, dir_1_gpu, dir_3_gpu, f_temp_6_gpu,
            h_gpu, b_gpu,
            distribution_x_6, distribution_y_6,
            g, d_x, d_y, d_t,
            size);
    honei::cuda::force_grid_gpu_y<<<grid, block>>>(
            dir_7_gpu, f_temp_7_gpu,
            h_gpu, b_gpu,
            distribution_y_7,
            g, d_x, d_y, d_t,
            size);
    honei::cuda::force_grid_gpu_xy<<<grid, block>>>(
            dir_8_gpu, dir_1_gpu, dir_3_gpu, f_temp_8_gpu,
            h_gpu, b_gpu,
            distribution_x_8, distribution_y_8,
            g, d_x, d_y, d_t,
            size);

    CUDA_ERROR();
}
