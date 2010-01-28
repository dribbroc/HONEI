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
        __global__ void collide_stream_grid_0_gpu(
                float * f_eq_0,
                float * f_0,
                float * f_temp_0,
                float tau,
                unsigned long offset_0, unsigned long size_0)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            if (idx >= offset_0 && idx < size_0)
            {
                unsigned long i(idx);
                //(*data.f_temp_0)[i] = (*data.f_0)[i] - ((*data.f_0)[i] - (*data.f_eq_0)[i])/tau;
                f_temp_0[i] = f_0[i] - (f_0[i] - f_eq_0[i]) / tau;
            }
        }
        __global__ void collide_stream_grid_n_gpu(
                unsigned long * dir,
                float * f_eq,
                float * f,
                float * f_temp,
                float tau,
                unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            if (idx < size)
            {
                unsigned long i(idx);
                if (dir[i] < size)
                    f_temp[dir[i]] = f[i] - (f[i] - f_eq[i])/tau;
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
    unsigned long size_0(end);
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

    honei::cuda::collide_stream_grid_0_gpu<<<grid, block>>>(
            f_eq_0_gpu,
            f_0_gpu,
            f_temp_0_gpu,
            tau,
            start, size_0);

    honei::cuda::collide_stream_grid_n_gpu<<<grid, block>>>(dir_1_gpu, f_eq_1_gpu, f_1_gpu, f_temp_1_gpu, tau, size);
    honei::cuda::collide_stream_grid_n_gpu<<<grid, block>>>(dir_2_gpu, f_eq_2_gpu, f_2_gpu, f_temp_2_gpu, tau, size);
    honei::cuda::collide_stream_grid_n_gpu<<<grid, block>>>(dir_3_gpu, f_eq_3_gpu, f_3_gpu, f_temp_3_gpu, tau, size);
    honei::cuda::collide_stream_grid_n_gpu<<<grid, block>>>(dir_4_gpu, f_eq_4_gpu, f_4_gpu, f_temp_4_gpu, tau, size);
    honei::cuda::collide_stream_grid_n_gpu<<<grid, block>>>(dir_5_gpu, f_eq_5_gpu, f_5_gpu, f_temp_5_gpu, tau, size);
    honei::cuda::collide_stream_grid_n_gpu<<<grid, block>>>(dir_6_gpu, f_eq_6_gpu, f_6_gpu, f_temp_6_gpu, tau, size);
    honei::cuda::collide_stream_grid_n_gpu<<<grid, block>>>(dir_7_gpu, f_eq_7_gpu, f_7_gpu, f_temp_7_gpu, tau, size);
    honei::cuda::collide_stream_grid_n_gpu<<<grid, block>>>(dir_8_gpu, f_eq_8_gpu, f_8_gpu, f_temp_8_gpu, tau, size);

    CUDA_ERROR();
}
