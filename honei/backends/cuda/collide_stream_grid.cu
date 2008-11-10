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
        __global__ void collide_stream_grid_gpu(
                unsigned long ** directions,
                float ** fs,
                float tau,
                unsigned long offset_0, unsigned long size_0, unsigned long size)
        {
            int idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            if (idx < size_0)
            {
                unsigned long i(idx + offset_0);
                //(*data.f_temp_0)[i] = (*data.f_0)[i] - ((*data.f_0)[i] - (*data.f_eq_0)[i])/tau;
                fs[18][i] = fs[9][i] - (fs[9][i] - fs[0][i]) / tau;
            }
            if (idx < size)
            {
                unsigned long i(idx);
                if (directions[0][i] < size)
                    fs[19][directions[0][i]] = fs[10][i] - (fs[10][i] - fs[1][i])/tau;
                if (directions[1][i] < size)
                    fs[20][directions[1][i]] = fs[11][i] - (fs[11][i] - fs[2][i])/tau;
                if (directions[2][i] < size)
                    fs[21][directions[2][i]] = fs[12][i] - (fs[12][i] - fs[3][i])/tau;
                if (directions[3][i] < size)
                    fs[22][directions[3][i]] = fs[13][i] - (fs[13][i] - fs[4][i])/tau;
                if (directions[4][i] < size)
                    fs[23][directions[4][i]] = fs[14][i] - (fs[14][i] - fs[5][i])/tau;
                if (directions[5][i] < size)
                    fs[24][directions[5][i]] = fs[15][i] - (fs[15][i] - fs[6][i])/tau;
                if (directions[6][i] < size)
                    fs[25][directions[6][i]] = fs[16][i] - (fs[16][i] - fs[7][i])/tau;
                if (directions[7][i] < size)
                    fs[26][directions[7][i]] = fs[17][i] - (fs[17][i] - fs[8][i])/tau;
            }
        }
    }
}

extern "C" void cuda_collide_stream_grid_float(unsigned long start, unsigned long end,
            void * dir_1, void * dir_2, void * dir_3, void * dir_4,
            void * dir_5, void * dir_6, void * dir_7, void * dir_8,
            void * f_eq_0, void * f_eq_1, void * f_eq_2,
            void * f_eq_3, void * f_eq_4, void * f_eq_5,
            void * f_eq_6, void * f_eq_7, void * f_eq_8,
            void * f_0, void * f_1, void * f_2,
            void * f_3, void * f_4, void * f_5,
            void * f_6, void * f_7, void * f_8,
            void * f_temp_0, void * f_temp_1, void * f_temp_2,
            void * f_temp_3, void * f_temp_4, void * f_temp_5,
            void * f_temp_6, void * f_temp_7, void * f_temp_8,
            float tau, unsigned long size,
            unsigned long blocksize)
{
    unsigned long size_0(end - start);
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

    float * f_eq_0_gpu((float *)f_eq_0);
    float * f_eq_1_gpu((float *)f_eq_1);
    float * f_eq_2_gpu((float *)f_eq_2);
    float * f_eq_3_gpu((float *)f_eq_3);
    float * f_eq_4_gpu((float *)f_eq_4);
    float * f_eq_5_gpu((float *)f_eq_5);
    float * f_eq_6_gpu((float *)f_eq_6);
    float * f_eq_7_gpu((float *)f_eq_7);
    float * f_eq_8_gpu((float *)f_eq_8);

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

    unsigned long * directions[8];
    directions[0] = dir_1_gpu;
    directions[1] = dir_2_gpu;
    directions[2] = dir_3_gpu;
    directions[3] = dir_4_gpu;
    directions[4] = dir_5_gpu;
    directions[5] = dir_6_gpu;
    directions[6] = dir_7_gpu;
    directions[7] = dir_8_gpu;

    float * fs[27];
    fs[0] = f_eq_0_gpu;
    fs[1] = f_eq_1_gpu;
    fs[2] = f_eq_2_gpu;
    fs[3] = f_eq_3_gpu;
    fs[4] = f_eq_4_gpu;
    fs[5] = f_eq_5_gpu;
    fs[6] = f_eq_6_gpu;
    fs[7] = f_eq_7_gpu;
    fs[8] = f_eq_8_gpu;
    fs[9] = f_0_gpu;
    fs[10] = f_1_gpu;
    fs[11] = f_2_gpu;
    fs[12] = f_3_gpu;
    fs[13] = f_4_gpu;
    fs[14] = f_5_gpu;
    fs[15] = f_6_gpu;
    fs[16] = f_7_gpu;
    fs[17] = f_8_gpu;
    fs[18] = f_temp_0_gpu;
    fs[19] = f_temp_1_gpu;
    fs[20] = f_temp_2_gpu;
    fs[21] = f_temp_3_gpu;
    fs[22] = f_temp_4_gpu;
    fs[23] = f_temp_5_gpu;
    fs[24] = f_temp_6_gpu;
    fs[25] = f_temp_7_gpu;
    fs[26] = f_temp_8_gpu;

    unsigned long ** directions_gpu;
    float ** fs_gpu;
    cudaMalloc((void **) &directions_gpu, sizeof(directions));
    cudaMalloc((void **) &fs_gpu, sizeof(fs));
    cudaMemcpy(directions_gpu, directions, sizeof(directions), cudaMemcpyHostToDevice);
    cudaMemcpy(fs_gpu, fs, sizeof(fs), cudaMemcpyHostToDevice);

    honei::cuda::collide_stream_grid_gpu<<<grid, block>>>(
            directions_gpu, fs_gpu,
            tau,
            start, size_0, size);

    cudaFree(directions_gpu);
    cudaFree(fs_gpu);

    CUDA_ERROR();
}
