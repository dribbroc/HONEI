/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2009 Markus Geveler <apryde@gmx.de>
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
        __global__ void collide_stream_fsi_1_gpu(
                unsigned long * dir,
                float * f_temp_i,
                float * f_temp_m_i,
                float * f_mea_m_i,
                bool * line_flags,
                float * dist_x,
                float * dist_y,
                float d_xu,
                float d_yv,
                unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            if (idx < size)
            {
                unsigned long i(idx);
                if(f_temp_i[i] != 0 && line_flags[i])
                {
                    f_temp_m_i[dir[i]] = f_temp_i[dir[i]] + 2./3. * (d_xu * dist_x[5] + d_yv * dist_y[5]);
                    f_mea_m_i[i] = (dist_x[5] + dist_y[5]) * (f_temp_m_i[dir[i]] + f_temp_i[i]);
                    f_temp_i[i] = 0.;
                }
            }
        }

        __global__ void collide_stream_fsi_2_gpu(
                unsigned long * dir,
                float * f_temp_i,
                float * f_temp_m_i,
                float * f_mea_m_i,
                bool * line_flags,
                float * dist_x,
                float * dist_y,
                float d_xu,
                float d_yv,
                unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            if (idx < size)
            {
                unsigned long i(idx);
                if(f_temp_i[i] != 0 && line_flags[i])
                {
                    f_temp_m_i[dir[i]] = f_temp_i[dir[i]] + 1./6. * (d_xu * dist_x[6] + d_yv * dist_y[6]);
                    f_mea_m_i[i] = (dist_x[6] + dist_y[6]) * (f_temp_m_i[dir[i]] + f_temp_i[i]);
                    f_temp_i[i] = 0.;
                }
            }
        }
        __global__ void collide_stream_fsi_3_gpu(
                unsigned long * dir,
                float * f_temp_i,
                float * f_temp_m_i,
                float * f_mea_m_i,
                bool * line_flags,
                float * dist_x,
                float * dist_y,
                float d_xu,
                float d_yv,
                unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            if (idx < size)
            {
                unsigned long i(idx);
                if(f_temp_i[i] != 0 && line_flags[i])
                {
                    f_temp_m_i[dir[i]] = f_temp_i[dir[i]] + 2./3. * (d_xu * dist_x[7] + d_yv * dist_y[7]);
                    f_mea_m_i[i] = (dist_x[7] + dist_y[7]) * (f_temp_m_i[dir[i]] + f_temp_i[i]);
                    f_temp_i[i] = 0.;
                }
            }
        }
        __global__ void collide_stream_fsi_4_gpu(
                unsigned long * dir,
                float * f_temp_i,
                float * f_temp_m_i,
                float * f_mea_m_i,
                bool * line_flags,
                float * dist_x,
                float * dist_y,
                float d_xu,
                float d_yv,
                unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            if (idx < size)
            {
                unsigned long i(idx);
                if(f_temp_i[i] != 0 && line_flags[i])
                {
                    f_temp_m_i[dir[i]] = f_temp_i[dir[i]] + 1./6. * (d_xu * dist_x[8] + d_yv * dist_y[8]);
                    f_mea_m_i[i] = (dist_x[8] + dist_y[8]) * (f_temp_m_i[dir[i]] + f_temp_i[i]);
                    f_temp_i[i] = 0.;
                }
            }
        }
        __global__ void collide_stream_fsi_5_gpu(
                unsigned long * dir,
                float * f_temp_i,
                float * f_temp_m_i,
                float * f_mea_m_i,
                bool * line_flags,
                float * dist_x,
                float * dist_y,
                float d_xu,
                float d_yv,
                unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            if (idx < size)
            {
                unsigned long i(idx);
                if(f_temp_i[i] != 0 && line_flags[i])
                {
                    f_temp_m_i[dir[i]] = f_temp_i[dir[i]] + 2./3. * (d_xu * dist_x[1] + d_yv * dist_y[1]);
                    f_mea_m_i[i] = (dist_x[1] + dist_y[1]) * (f_temp_m_i[dir[i]] + f_temp_i[i]);
                    f_temp_i[i] = 0.;
                }
            }
        }
        __global__ void collide_stream_fsi_6_gpu(
                unsigned long * dir,
                float * f_temp_i,
                float * f_temp_m_i,
                float * f_mea_m_i,
                bool * line_flags,
                float * dist_x,
                float * dist_y,
                float d_xu,
                float d_yv,
                unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            if (idx < size)
            {
                unsigned long i(idx);
                if(f_temp_i[i] != 0 && line_flags[i])
                {
                    f_temp_m_i[dir[i]] = f_temp_i[dir[i]] + 1./6. * (d_xu * dist_x[2] + d_yv * dist_y[2]);
                    f_mea_m_i[i] = (dist_x[2] + dist_y[2]) * (f_temp_m_i[dir[i]] + f_temp_i[i]);
                    f_temp_i[i] = 0.;
                }
            }
        }
        __global__ void collide_stream_fsi_7_gpu(
                unsigned long * dir,
                float * f_temp_i,
                float * f_temp_m_i,
                float * f_mea_m_i,
                bool * line_flags,
                float * dist_x,
                float * dist_y,
                float d_xu,
                float d_yv,
                unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            if (idx < size)
            {
                unsigned long i(idx);
                if(f_temp_i[i] != 0 && line_flags[i])
                {
                    f_temp_m_i[dir[i]] = f_temp_i[dir[i]] + 2./3. * (d_xu * dist_x[3] + d_yv * dist_y[3]);
                    f_mea_m_i[i] = (dist_x[3] + dist_y[3]) * (f_temp_m_i[dir[i]] + f_temp_i[i]);
                    f_temp_i[i] = 0.;
                }
            }
        }
        __global__ void collide_stream_fsi_8_gpu(
                unsigned long * dir,
                float * f_temp_i,
                float * f_temp_m_i,
                float * f_mea_m_i,
                bool * line_flags,
                float * dist_x,
                float * dist_y,
                float d_xu,
                float d_yv,
                unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            if (idx < size)
            {
                unsigned long i(idx);
                if(f_temp_i[i] != 0 && line_flags[i])
                {
                    f_temp_m_i[dir[i]] = f_temp_i[dir[i]] + 1./6. * (d_xu * dist_x[4] + d_yv * dist_y[4]);
                    f_mea_m_i[i] = (dist_x[4] + dist_y[4]) * (f_temp_m_i[dir[i]] + f_temp_i[i]);
                    f_temp_i[i] = 0.;
                }
            }
        }
    }
}

extern "C" void cuda_collide_stream_fsi_float(unsigned long start, unsigned long end,
            void * dir_1, void * dir_2, void * dir_3, void * dir_4,
            void * dir_5, void * dir_6, void * dir_7, void * dir_8,
            void * f_temp_1, void * f_temp_2,
            void * f_temp_3, void * f_temp_4, void * f_temp_5,
            void * f_temp_6, void * f_temp_7, void * f_temp_8,
            void * f_mea_1, void * f_mea_2,
            void * f_mea_3, void * f_mea_4, void * f_mea_5,
            void * f_mea_6, void * f_mea_7, void * f_mea_8,
            void * line_flags, void * dist_x, void * dist_y,
            float d_xu, float d_yv,
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

    float * f_temp_1_gpu((float *)f_temp_1);
    float * f_temp_2_gpu((float *)f_temp_2);
    float * f_temp_3_gpu((float *)f_temp_3);
    float * f_temp_4_gpu((float *)f_temp_4);
    float * f_temp_5_gpu((float *)f_temp_5);
    float * f_temp_6_gpu((float *)f_temp_6);
    float * f_temp_7_gpu((float *)f_temp_7);
    float * f_temp_8_gpu((float *)f_temp_8);

    float * f_mea_1_gpu((float *)f_mea_1);
    float * f_mea_2_gpu((float *)f_mea_2);
    float * f_mea_3_gpu((float *)f_mea_3);
    float * f_mea_4_gpu((float *)f_mea_4);
    float * f_mea_5_gpu((float *)f_mea_5);
    float * f_mea_6_gpu((float *)f_mea_6);
    float * f_mea_7_gpu((float *)f_mea_7);
    float * f_mea_8_gpu((float *)f_mea_8);

    bool * line_flags_gpu((bool *)line_flags);
    float * dist_x_gpu((float *)dist_x);
    float * dist_y_gpu((float *)dist_y);

    honei::cuda::collide_stream_fsi_1_gpu<<<grid, block>>>(dir_1_gpu,
                                                           f_temp_1_gpu,
                                                           f_temp_5_gpu,
                                                           f_mea_5_gpu,
                                                           line_flags_gpu,
                                                           dist_x_gpu,
                                                           dist_y_gpu,
                                                           d_xu,
                                                           d_yv,
                                                           size);

    honei::cuda::collide_stream_fsi_2_gpu<<<grid, block>>>(dir_2_gpu,
                                                           f_temp_2_gpu,
                                                           f_temp_6_gpu,
                                                           f_mea_6_gpu,
                                                           line_flags_gpu,
                                                           dist_x_gpu,
                                                           dist_y_gpu,
                                                           d_xu,
                                                           d_yv,
                                                           size);

    honei::cuda::collide_stream_fsi_3_gpu<<<grid, block>>>(dir_3_gpu,
                                                           f_temp_3_gpu,
                                                           f_temp_7_gpu,
                                                           f_mea_7_gpu,
                                                           line_flags_gpu,
                                                           dist_x_gpu,
                                                           dist_y_gpu,
                                                           d_xu,
                                                           d_yv,
                                                           size);

    honei::cuda::collide_stream_fsi_4_gpu<<<grid, block>>>(dir_4_gpu,
                                                           f_temp_4_gpu,
                                                           f_temp_8_gpu,
                                                           f_mea_8_gpu,
                                                           line_flags_gpu,
                                                           dist_x_gpu,
                                                           dist_y_gpu,
                                                           d_xu,
                                                           d_yv,
                                                           size);

    honei::cuda::collide_stream_fsi_5_gpu<<<grid, block>>>(dir_5_gpu,
                                                           f_temp_5_gpu,
                                                           f_temp_1_gpu,
                                                           f_mea_1_gpu,
                                                           line_flags_gpu,
                                                           dist_x_gpu,
                                                           dist_y_gpu,
                                                           d_xu,
                                                           d_yv,
                                                           size);

    honei::cuda::collide_stream_fsi_6_gpu<<<grid, block>>>(dir_6_gpu,
                                                           f_temp_6_gpu,
                                                           f_temp_2_gpu,
                                                           f_mea_2_gpu,
                                                           line_flags_gpu,
                                                           dist_x_gpu,
                                                           dist_y_gpu,
                                                           d_xu,
                                                           d_yv,
                                                           size);

    honei::cuda::collide_stream_fsi_7_gpu<<<grid, block>>>(dir_7_gpu,
                                                           f_temp_7_gpu,
                                                           f_temp_3_gpu,
                                                           f_mea_3_gpu,
                                                           line_flags_gpu,
                                                           dist_x_gpu,
                                                           dist_y_gpu,
                                                           d_xu,
                                                           d_yv,
                                                           size);

    honei::cuda::collide_stream_fsi_8_gpu<<<grid, block>>>(dir_8_gpu,
                                                           f_temp_8_gpu,
                                                           f_temp_4_gpu,
                                                           f_mea_4_gpu,
                                                           line_flags_gpu,
                                                           dist_x_gpu,
                                                           dist_y_gpu,
                                                           d_xu,
                                                           d_yv,
                                                           size);
    CUDA_ERROR();
}
