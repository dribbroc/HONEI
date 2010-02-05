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
#include <iostream>
namespace honei
{
    namespace cuda
    {
        __global__ void up_vel_dir_grid_gpu_ordinary(
                float * f_temp_backward,
                float * f_forward,
                float * f_eq_forward,
                unsigned long * types,
                unsigned long type,
                float tau,
                unsigned long offset,
                unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            if (idx >= offset && idx < size)
            {
                unsigned long i(idx);
                if((types[i] & 1<<(type)) == 1<<(type))
                    f_temp_backward[i] = f_forward[i] - (f_forward[i] - f_eq_forward[i])/tau;
            }
        }
        __global__ void up_vel_dir_grid_gpu_corner(
                float * f_temp_old_1,
                float * f_temp_old_2,
                float * f_temp_new,
                unsigned long * types,
                unsigned long type_1,
                unsigned long type_2,
                unsigned long type_3,
                unsigned long type_4,
                float tau,
                unsigned long offset,
                unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            if (idx >= offset && idx < size)
            {
                unsigned long i(idx);
                if((types[i] & 1<<(type_1)) == 1<<(type_1) &&
                        (types[i] & 1<<(type_2)) == 1<<(type_2) &&
                        (types[i] & 1<<(type_3)) == 1<<(type_3) &&
                        (types[i] & 1<<(type_4)) == 1<<(type_4)
                  )
                {
                    f_temp_old_1[i] = f_temp_new[i];
                    f_temp_old_2[i] = f_temp_new[i];
                }
            }
        }
        /*__global__ void up_vel_dir_grid_gpu(
                float * f_temp_1, float * f_temp_2,
                float * f_temp_3, float * f_temp_4, float * f_temp_5,
                float * f_temp_6, float * f_temp_7, float * f_temp_8,
                unsigned long * types,
                float tau,
                unsigned long offset,
                unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            if (idx >= offset && idx < size)
            {
                unsigned long i(idx);
                if((types[i] & 1<<0) == 1<<0)
                    f_temp_5[i] = f_temp_1[i];
                if((types[i] & 1<<1) == 1<<1)
                    f_temp_6[i] = f_temp_2[i];
                if((types[i] & 1<<2) == 1<<2)
                    f_temp_7[i] = f_temp_3[i];
                if((types[i] & 1<<3) == 1<<3)
                    f_temp_8[i] = f_temp_4[i];
                if((types[i] & 1<<4) == 1<<4)
                    f_temp_1[i] = f_temp_5[i];
                if((types[i] & 1<<5) == 1<<5)
                    f_temp_2[i] = f_temp_6[i];
                if((types[i] & 1<<6) == 1<<6)
                    f_temp_3[i] = f_temp_7[i];
                if((types[i] & 1<<7) == 1<<7)
                    f_temp_4[i] = f_temp_8[i];

                if(((types)[i] & 1<<2) == 1<<2 &&
                        ((types)[i] & 1<<4) == 1<<4 &&
                        ((types)[i] & 1<<1) == 1<<1 &&
                        ((types)[i] & 1<<5) == 1<<5)
                {
                    (f_temp_2)[i] = (f_temp_8)[i];
                    (f_temp_6)[i] = (f_temp_8)[i];
                }
                if(((types)[i] & 1<<4) == 1<<4 &&
                        ((types)[i] & 1<<6) == 1<<6 &&
                        ((types)[i] & 1<<7) == 1<<7 &&
                        ((types)[i] & 1<<3) == 1<<3)
                {
                    (f_temp_4)[i] = (f_temp_2)[i];
                    (f_temp_8)[i] = (f_temp_2)[i];
                }
                if(((types)[i] & 1<<0) == 1<<0 &&
                        ((types)[i] & 1<<6) == 1<<6 &&
                        ((types)[i] & 1<<1) == 1<<1 &&
                        ((types)[i] & 1<<5) == 1<<5)
                {
                    (f_temp_2)[i] = (f_temp_4)[i];
                    (f_temp_6)[i] = (f_temp_4)[i];
                }
                if(((types)[i] & 1<<0) == 1<<0 &&
                        ((types)[i] & 1<<2) == 1<<2 &&
                        ((types)[i] & 1<<7) == 1<<7 &&
                        ((types)[i] & 1<<3) == 1<<3)
                {
                    (f_temp_4)[i] = (f_temp_6)[i];
                    (f_temp_8)[i] = (f_temp_6)[i];
                }
            }
        }*/
    }
}

extern "C" void cuda_up_vel_dir_grid_float(unsigned long start, unsigned long end,
        void * types, void * f_temp_1, void * f_temp_2,
        void * f_temp_3, void * f_temp_4, void * f_temp_5,
        void * f_temp_6, void * f_temp_7, void * f_temp_8,
        void * f_1, void * f_2, void * f_3, void * f_4, void * f_5, void * f_6, void * f_7, void * f_8,
        void * f_eq_1, void * f_eq_2, void * f_eq_3, void * f_eq_4, void * f_eq_5, void * f_eq_6, void * f_eq_7, void * f_eq_8,
        float tau,
        unsigned long blocksize)
{
    unsigned long size(end);
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

    float * f_1_gpu((float *)f_1);
    float * f_2_gpu((float *)f_2);
    float * f_3_gpu((float *)f_3);
    float * f_4_gpu((float *)f_4);
    float * f_5_gpu((float *)f_5);
    float * f_6_gpu((float *)f_6);
    float * f_7_gpu((float *)f_7);
    float * f_8_gpu((float *)f_8);

    float * f_eq_1_gpu((float *)f_eq_1);
    float * f_eq_2_gpu((float *)f_eq_2);
    float * f_eq_3_gpu((float *)f_eq_3);
    float * f_eq_4_gpu((float *)f_eq_4);
    float * f_eq_5_gpu((float *)f_eq_5);
    float * f_eq_6_gpu((float *)f_eq_6);
    float * f_eq_7_gpu((float *)f_eq_7);
    float * f_eq_8_gpu((float *)f_eq_8);
    /*honei::cuda::up_vel_dir_grid_gpu<<<grid, block>>>(
            f_temp_1_gpu, f_temp_2_gpu, f_temp_3_gpu, f_temp_4_gpu,
            f_temp_5_gpu, f_temp_6_gpu, f_temp_7_gpu, f_temp_8_gpu,
            types_gpu,
            tau,
            start, size);*/
    honei::cuda::up_vel_dir_grid_gpu_ordinary<<<grid, block>>>(
            f_temp_5_gpu,
            f_1_gpu,
            f_eq_1_gpu,
            types_gpu,
            0ul,
            tau,
            start, size);
    honei::cuda::up_vel_dir_grid_gpu_ordinary<<<grid, block>>>(
            f_temp_6_gpu,
            f_2_gpu,
            f_eq_2_gpu,
            types_gpu,
            1ul,
            tau,
            start, size);
    honei::cuda::up_vel_dir_grid_gpu_ordinary<<<grid, block>>>(
            f_temp_7_gpu,
            f_3_gpu,
            f_eq_3_gpu,
            types_gpu,
            2ul,
            tau,
            start, size);
    honei::cuda::up_vel_dir_grid_gpu_ordinary<<<grid, block>>>(
            f_temp_8_gpu,
            f_4_gpu,
            f_eq_4_gpu,
            types_gpu,
            3ul,
            tau,
            start, size);
    honei::cuda::up_vel_dir_grid_gpu_ordinary<<<grid, block>>>(
            f_temp_1_gpu,
            f_5_gpu,
            f_eq_5_gpu,
            types_gpu,
            4ul,
            tau,
            start, size);
    honei::cuda::up_vel_dir_grid_gpu_ordinary<<<grid, block>>>(
            f_temp_2_gpu,
            f_6_gpu,
            f_eq_6_gpu,
            types_gpu,
            5ul,
            tau,
            start, size);
    honei::cuda::up_vel_dir_grid_gpu_ordinary<<<grid, block>>>(
            f_temp_3_gpu,
            f_7_gpu,
            f_eq_7_gpu,
            types_gpu,
            6ul,
            tau,
            start, size);
    honei::cuda::up_vel_dir_grid_gpu_ordinary<<<grid, block>>>(
            f_temp_4_gpu,
            f_8_gpu,
            f_eq_8_gpu,
            types_gpu,
            7ul,
            tau,
            start, size);

    honei::cuda::up_vel_dir_grid_gpu_corner<<<grid, block>>>(
            f_temp_2_gpu,
            f_temp_6_gpu,
            f_temp_8_gpu,
            types_gpu,
            1ul,
            2ul,
            4ul,
            5ul,
            tau,
            start, size);
    honei::cuda::up_vel_dir_grid_gpu_corner<<<grid, block>>>(
            f_temp_4_gpu,
            f_temp_8_gpu,
            f_temp_2_gpu,
            types_gpu,
            3ul,
            4ul,
            6ul,
            7ul,
            tau,
            start, size);
    honei::cuda::up_vel_dir_grid_gpu_corner<<<grid, block>>>(
            f_temp_2_gpu,
            f_temp_6_gpu,
            f_temp_4_gpu,
            types_gpu,
            0ul,
            1ul,
            5ul,
            6ul,
            tau,
            start, size);
    honei::cuda::up_vel_dir_grid_gpu_corner<<<grid, block>>>(
            f_temp_4_gpu,
            f_temp_8_gpu,
            f_temp_6_gpu,
            types_gpu,
            0ul,
            2ul,
            3ul,
            7ul,
            tau,
            start, size);

    CUDA_ERROR();
}
