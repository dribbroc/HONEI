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
        __global__ void extraction_grid_gpu(
                float ** fs,
                float * h, float * u, float * v,
                float * distribution_x, float * distribution_y,
                unsigned long offset, unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            if (idx < size)
            {
                unsigned long i(idx + offset);

                fs[0][i] = fs[9][i];
                fs[1][i] = fs[10][i];
                fs[2][i] = fs[11][i];
                fs[3][i] = fs[12][i];
                fs[4][i] = fs[13][i];
                fs[5][i] = fs[14][i];
                fs[6][i] = fs[15][i];
                fs[7][i] = fs[16][i];
                fs[8][i] = fs[17][i];

                h[i] = fs[0][i] +
                    fs[1][i] +
                    fs[2][i] +
                    fs[3][i] +
                    fs[4][i] +
                    fs[5][i] +
                    fs[6][i] +
                    fs[7][i] +
                    fs[8][i];

                u[i] = (distribution_x[0] * fs[0][i] +
                        distribution_x[1] * fs[1][i] +
                        distribution_x[2] * fs[2][i] +
                        distribution_x[3] * fs[3][i] +
                        distribution_x[4] * fs[4][i] +
                        distribution_x[5] * fs[5][i] +
                        distribution_x[6] * fs[6][i] +
                        distribution_x[7] * fs[7][i] +
                        distribution_x[8] * fs[8][i]) / h[i];

                v[i] = (distribution_y[0] * fs[0][i] +
                        distribution_y[1] * fs[1][i] +
                        distribution_y[2] * fs[2][i] +
                        distribution_y[3] * fs[3][i] +
                        distribution_y[4] * fs[4][i] +
                        distribution_y[5] * fs[5][i] +
                        distribution_y[6] * fs[6][i] +
                        distribution_y[7] * fs[7][i] +
                        distribution_y[8] * fs[8][i]) / h[i];
            }
        }
    }
}

extern "C" void cuda_extraction_grid_float(
        unsigned long start, unsigned long end,
        void * f_0, void * f_1, void * f_2,
        void * f_3, void * f_4, void * f_5,
        void * f_6, void * f_7, void * f_8,
        void * f_temp_0, void * f_temp_1, void * f_temp_2,
        void * f_temp_3, void * f_temp_4, void * f_temp_5,
        void * f_temp_6, void * f_temp_7, void * f_temp_8,
        void * h, void * u, void * v,
        void * distribution_x, void * distribution_y,
        unsigned long blocksize)
{
    unsigned long size(end - start);
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

    float * f_temp_0_gpu((float *)f_temp_0);
    float * f_temp_1_gpu((float *)f_temp_1);
    float * f_temp_2_gpu((float *)f_temp_2);
    float * f_temp_3_gpu((float *)f_temp_3);
    float * f_temp_4_gpu((float *)f_temp_4);
    float * f_temp_5_gpu((float *)f_temp_5);
    float * f_temp_6_gpu((float *)f_temp_6);
    float * f_temp_7_gpu((float *)f_temp_7);
    float * f_temp_8_gpu((float *)f_temp_8);

    float * h_gpu((float *)h);
    float * u_gpu((float *)u);
    float * v_gpu((float *)v);
    float * distribution_x_gpu((float *)distribution_x);
    float * distribution_y_gpu((float *)distribution_y);

    float * fs[18];
    fs[0] = f_0_gpu;
    fs[1] = f_1_gpu;
    fs[2] = f_2_gpu;
    fs[3] = f_3_gpu;
    fs[4] = f_4_gpu;
    fs[5] = f_5_gpu;
    fs[6] = f_6_gpu;
    fs[7] = f_7_gpu;
    fs[8] = f_8_gpu;
    fs[9] = f_temp_0_gpu;
    fs[10] = f_temp_1_gpu;
    fs[11] = f_temp_2_gpu;
    fs[12] = f_temp_3_gpu;
    fs[13] = f_temp_4_gpu;
    fs[14] = f_temp_5_gpu;
    fs[15] = f_temp_6_gpu;
    fs[16] = f_temp_7_gpu;
    fs[17] = f_temp_8_gpu;

    float ** fs_gpu;
    cudaMalloc((void **) &fs_gpu, sizeof(fs));
    cudaMemcpy(fs_gpu, fs, sizeof(fs), cudaMemcpyHostToDevice);

    honei::cuda::extraction_grid_gpu<<<grid, block>>>(
            fs_gpu,
            h_gpu, u_gpu, v_gpu,
            distribution_x_gpu, distribution_y_gpu,
            start, size);

    cudaFree(fs_gpu);

    CUDA_ERROR();
}
