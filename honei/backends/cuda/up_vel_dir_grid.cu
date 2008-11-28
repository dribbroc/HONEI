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
        __global__ void up_vel_dir_grid_gpu(
                float ** fs,
                unsigned long * types,
                unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            if (idx < size)
            {
                unsigned long i(idx);
                if((types[i] & 1<<0) == 1<<0)
                    fs[4][i] = fs[0][i];
                if((types[i] & 1<<1) == 1<<1)
                    fs[5][i] = fs[1][i];
                if((types[i] & 1<<2) == 1<<2)
                    fs[6][i] = fs[2][i];
                if((types[i] & 1<<3) == 1<<3)
                    fs[7][i] = fs[3][i];
                if((types[i] & 1<<4) == 1<<4)
                    fs[0][i] = fs[4][i];
                if((types[i] & 1<<5) == 1<<5)
                    fs[1][i] = fs[5][i];
                if((types[i] & 1<<6) == 1<<6)
                    fs[2][i] = fs[6][i];
                if((types[i] & 1<<7) == 1<<7)
                    fs[3][i] = fs[7][i];

                // Corners
                if((types[i] & 1<<2) == 1<<2 && (types[i] & 1<<4) == 1<<4)
                {
                    fs[1][i] = fs[7][i];
                    fs[5][i] = fs[7][i];
                }
                if((types[i] & 1<<4) == 1<<4 && (types[i] & 1<<6) == 1<<6)
                {
                    fs[3][i] = fs[1][i];
                    fs[7][i] = fs[1][i];
                }
                if((types[i] & 1<<0) == 1<<0 && (types[i] & 1<<6) == 1<<6)
                {
                    fs[1][i] = fs[3][i];
                    fs[5][i] = fs[3][i];
                }
                if((types[i] & 1<<0) == 1<<0 && (types[i] & 1<<2) == 1<<2)
                {
                    fs[3][i] = fs[5][i];
                    fs[7][i] = fs[5][i];
                }
            }
        }
    }
}

extern "C" void cuda_up_vel_dir_grid_float(void * types,
        void * f_temp_1, void * f_temp_2,
        void * f_temp_3, void * f_temp_4, void * f_temp_5,
        void * f_temp_6, void * f_temp_7, void * f_temp_8,
        unsigned long size,
        unsigned long blocksize)
{
    dim3 grid;
    dim3 block;
    block.x = blocksize;
    grid.x = (unsigned)ceil(sqrt(size/(double)block.x));
    grid.y = grid.x;


    unsigned long * types_gpu((unsigned long *)types);

    float * f_temp_1_gpu((float *)f_temp_1);
    float * f_temp_2_gpu((float *)f_temp_2);
    float * f_temp_3_gpu((float *)f_temp_3);
    float * f_temp_4_gpu((float *)f_temp_4);
    float * f_temp_5_gpu((float *)f_temp_5);
    float * f_temp_6_gpu((float *)f_temp_6);
    float * f_temp_7_gpu((float *)f_temp_7);
    float * f_temp_8_gpu((float *)f_temp_8);

    float * fs[8];
    fs[0] = f_temp_1_gpu;
    fs[1] = f_temp_2_gpu;
    fs[2] = f_temp_3_gpu;
    fs[3] = f_temp_4_gpu;
    fs[4] = f_temp_5_gpu;
    fs[5] = f_temp_6_gpu;
    fs[6] = f_temp_7_gpu;
    fs[7] = f_temp_8_gpu;

    float ** fs_gpu;
    cudaMalloc((void **) &fs_gpu, sizeof(fs));
    cudaMemcpy(fs_gpu, fs, sizeof(fs), cudaMemcpyHostToDevice);

    honei::cuda::up_vel_dir_grid_gpu<<<grid, block>>>(
            fs_gpu,
            types_gpu,
            size);

    cudaFree(fs_gpu);

    CUDA_ERROR();
}
