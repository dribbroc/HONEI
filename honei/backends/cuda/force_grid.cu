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
        __global__ void force_grid_gpu(
                unsigned long ** directions,
                float ** fs,
                float * h, float * b_x, float * b_y,
                float * distribution_x, float * distribution_y,
                float g, float d_x, float d_y, float d_t,
                unsigned long size)
        {
            int idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            if (idx < size)
            {
                /*(*data.f_temp_1)[i] += d_t / (6 * d_x / d_t) * ((*data.distribution_x)[1]) *
                    (-g * (((*data.h)[(*info.dir_1)[half] + offset]) - ((*data.b_x)[i]))/ DT1_(2.) *
                     ((((*data.b_x)[(*info.dir_1)[half] + offset]) - ((*data.b_x)[i]))/ DT1_(2.)));

                (*data.f_temp_1)[i] += d_t / (6 * d_y / d_t) * ((*data.distribution_y)[1]) *
                    (- g * (((*data.h)[(*info.dir_1)[half] + offset]) - ((*data.b_y)[i]))/ DT1_(2.) *
                     ((((*data.b_y)[(*info.dir_1)[half] + offset]) - ((*data.b_y)[i]))/ DT1_(2.)));*/
                unsigned long i(idx);
                if (directions[0][i] < size)
                {
                    fs[0][i] += d_t / (6 * d_x / d_t) * distribution_x[1] *
                        (-g * (h[directions[0][i]] - b_x[i])/ float(2.) *
                         ((b_x[directions[0][i]] - b_x[i])/ float(2.)));

                    fs[0][i] += d_t / (6 * d_y / d_t) * distribution_y[1] *
                        (-g * (h[directions[0][i]] - b_y[i])/ float(2.) *
                         ((b_y[directions[0][i]] - b_y[i])/ float(2.)));
                }
                if (directions[1][i] < size)
                {
                    fs[1][i] += d_t / (6 * d_x / d_t) * distribution_x[2] *
                        (-g * (h[directions[1][i]] - b_x[i])/ float(2.) *
                         ((b_x[directions[1][i]] - b_x[i])/ float(2.)));

                    fs[1][i] += d_t / (6 * d_y / d_t) * distribution_y[2] *
                        (-g * (h[directions[1][i]] - b_y[i])/ float(2.) *
                         ((b_y[directions[1][i]] - b_y[i])/ float(2.)));
                }
                if (directions[2][i] < size)
                {
                    fs[2][i] += d_t / (6 * d_x / d_t) * distribution_x[3] *
                        (-g * (h[directions[2][i]] - b_x[i])/ float(2.) *
                         ((b_x[directions[2][i]] - b_x[i])/ float(2.)));

                    fs[2][i] += d_t / (6 * d_y / d_t) * distribution_y[3] *
                        (-g * (h[directions[2][i]] - b_y[i])/ float(2.) *
                         ((b_y[directions[2][i]] - b_y[i])/ float(2.)));
                }
                if (directions[3][i] < size)
                {
                    fs[3][i] += d_t / (6 * d_x / d_t) * distribution_x[4] *
                        (-g * (h[directions[3][i]] - b_x[i])/ float(2.) *
                         ((b_x[directions[3][i]] - b_x[i])/ float(2.)));

                    fs[3][i] += d_t / (6 * d_y / d_t) * distribution_y[4] *
                        (-g * (h[directions[3][i]] - b_y[i])/ float(2.) *
                         ((b_y[directions[3][i]] - b_y[i])/ float(2.)));
                }
                if (directions[4][i] < size)
                {
                    fs[4][i] += d_t / (6 * d_x / d_t) * distribution_x[5] *
                        (-g * (h[directions[4][i]] - b_x[i])/ float(2.) *
                         ((b_x[directions[4][i]] - b_x[i])/ float(2.)));

                    fs[4][i] += d_t / (6 * d_y / d_t) * distribution_y[5] *
                        (-g * (h[directions[4][i]] - b_y[i])/ float(2.) *
                         ((b_y[directions[4][i]] - b_y[i])/ float(2.)));
                }
                if (directions[5][i] < size)
                {
                    fs[5][i] += d_t / (6 * d_x / d_t) * distribution_x[6] *
                        (-g * (h[directions[5][i]] - b_x[i])/ float(2.) *
                         ((b_x[directions[5][i]] - b_x[i])/ float(2.)));

                    fs[5][i] += d_t / (6 * d_y / d_t) * distribution_y[6] *
                        (-g * (h[directions[5][i]] - b_y[i])/ float(2.) *
                         ((b_y[directions[5][i]] - b_y[i])/ float(2.)));
                }
                if (directions[6][i] < size)
                {
                    fs[6][i] += d_t / (6 * d_x / d_t) * distribution_x[7] *
                        (-g * (h[directions[6][i]] - b_x[i])/ float(2.) *
                         ((b_x[directions[6][i]] - b_x[i])/ float(2.)));

                    fs[6][i] +=  d_t / (6 * d_y / d_t) * distribution_y[7] *
                        (-g * (h[directions[6][i]] - b_y[i])/ float(2.) *
                         ((b_y[directions[6][i]] - b_y[i])/ float(2.)));
                }
                if (directions[7][i] < size)
                {
                    fs[7][i] += d_t / (6 * d_x / d_t) * distribution_x[8] *
                        (-g * (h[directions[7][i]] - b_x[i])/ float(2.) *
                         ((b_x[directions[7][i]] - b_x[i])/ float(2.)));

                    fs[7][i] += d_t / (6 * d_y / d_t) * distribution_y[8] *
                        (-g * (h[directions[7][i]] - b_y[i])/ float(2.) *
                         ((b_y[directions[7][i]] - b_y[i])/ float(2.)));
                }
            }
        }
    }
}

extern "C" void cuda_force_grid_float(
        void * dir_1, void * dir_2, void * dir_3, void * dir_4,
        void * dir_5, void * dir_6, void * dir_7, void * dir_8,
        void * h, void * b_x, void * b_y,
        void * distribution_x, void * distribution_y,
        void * f_temp_1, void * f_temp_2,
        void * f_temp_3, void * f_temp_4, void * f_temp_5,
        void * f_temp_6, void * f_temp_7, void * f_temp_8,
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
    float * b_x_gpu((float *)b_x);
    float * b_y_gpu((float *)b_y);
    float * distribution_x_gpu((float *)distribution_x);
    float * distribution_y_gpu((float *)distribution_y);

    float * f_temp_1_gpu((float *)f_temp_1);
    float * f_temp_2_gpu((float *)f_temp_2);
    float * f_temp_3_gpu((float *)f_temp_3);
    float * f_temp_4_gpu((float *)f_temp_4);
    float * f_temp_5_gpu((float *)f_temp_5);
    float * f_temp_6_gpu((float *)f_temp_6);
    float * f_temp_7_gpu((float *)f_temp_7);
    float * f_temp_8_gpu((float *)f_temp_8);

    unsigned long * directions[8];
    directions[0] = dir_1_gpu;
    directions[1] = dir_2_gpu;
    directions[2] = dir_3_gpu;
    directions[3] = dir_4_gpu;
    directions[4] = dir_5_gpu;
    directions[5] = dir_6_gpu;
    directions[6] = dir_7_gpu;
    directions[7] = dir_8_gpu;

    float * fs[8];
    fs[0] = f_temp_1_gpu;
    fs[1] = f_temp_2_gpu;
    fs[2] = f_temp_3_gpu;
    fs[3] = f_temp_4_gpu;
    fs[4] = f_temp_5_gpu;
    fs[5] = f_temp_6_gpu;
    fs[6] = f_temp_7_gpu;
    fs[7] = f_temp_8_gpu;

    unsigned long ** directions_gpu;
    float ** fs_gpu;
    cudaMalloc((void **) &directions_gpu, sizeof(directions));
    cudaMalloc((void **) &fs_gpu, sizeof(fs));
    cudaMemcpy(directions_gpu, directions, sizeof(directions), cudaMemcpyHostToDevice);
    cudaMemcpy(fs_gpu, fs, sizeof(fs), cudaMemcpyHostToDevice);

    honei::cuda::force_grid_gpu<<<grid, block>>>(
            directions_gpu, fs_gpu,
            h_gpu, b_x_gpu, b_y_gpu,
            distribution_x_gpu, distribution_y_gpu,
            g, d_x, d_y, d_t,
            size);

    cudaFree(directions_gpu);
    cudaFree(fs_gpu);

    CUDA_ERROR();
}
