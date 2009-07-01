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
        __global__ void boundary_init_fsi_1_set_gpu(
                unsigned long * dir_5,
                bool * boundary_flags,
                bool * solid_flags,
                bool * solid_old_flags,
                float * h,
                float * u,
                float * v,
                float c_u,
                float c_v,
                unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            if (idx < size)
            {
                unsigned long i(idx);

                bool stf(solid_old_flags[i] & !solid_flags[i]);
                bool fts(!solid_old_flags[i] & solid_flags[i]);

                bool prev((dir_5[i] < size) ? (!solid_flags[dir_5[i]] ? true : false) : false);
                unsigned long prev_index(prev ? dir_5[i] : i);

                bool pre_prev(prev ? ((dir_5[prev_index] < size) ?
                            (!solid_flags[dir_5[prev_index]] ? true : false) : false) : false);
                unsigned long pre_prev_index(pre_prev ? dir_5[prev_index] : prev_index);

                bool pre_pre_prev(pre_prev ? (dir_5[pre_prev_index] < size ?
                            (!solid_flags[dir_5[pre_prev_index]] ? true : false) : false) : false);
                unsigned long pre_pre_prev_index(pre_pre_prev ? dir_5[pre_prev_index] : pre_prev_index);

                h[i] = stf ? (float(3.) * (h[prev_index] - h[pre_prev_index]) + h[pre_pre_prev_index]) : (fts ? float(0) : h[i]);

                u[i] = (boundary_flags[i] & !stf) ? c_u  :
                    (stf ? (float(8./15.) * c_u) + (float(2./3.) * u[prev_index]) - (float(2./5.) * u[pre_prev_index]) :(fts ? float(0) : u[i]));

                v[i] = (boundary_flags[i] & !stf) ? c_v  :
                    (stf ? (float(8./15.) * c_v) + (float(2./3.) * v[prev_index]) - (float(2./5.) * v[pre_prev_index]) :(fts ? float(0) : v[i]));

            }
        }
        __global__ void boundary_init_fsi_1_reset_f_gpu(
                float * f_0,
                float * f_1,
                float * f_2,
                float * f_3,
                float * f_4,
                float * f_5,
                float * f_6,
                float * f_7,
                float * f_8,
                bool * solid_flags,
                bool * solid_old_flags,
                unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            if (idx < size)
            {
                unsigned long i(idx);

                bool fts(!solid_old_flags[i] & solid_flags[i]);
                f_0[i] = fts ? float(0) : f_0[i];
                f_1[i] = fts ? float(0) : f_1[i];
                f_2[i] = fts ? float(0) : f_2[i];
                f_3[i] = fts ? float(0) : f_3[i];
                f_4[i] = fts ? float(0) : f_4[i];
                f_5[i] = fts ? float(0) : f_5[i];
                f_6[i] = fts ? float(0) : f_6[i];
                f_7[i] = fts ? float(0) : f_7[i];
                f_8[i] = fts ? float(0) : f_8[i];
            }
        }
        __global__ void boundary_init_fsi_1_reset_f_eq_gpu(
                float * f_eq_0,
                float * f_eq_1,
                float * f_eq_2,
                float * f_eq_3,
                float * f_eq_4,
                float * f_eq_5,
                float * f_eq_6,
                float * f_eq_7,
                float * f_eq_8,
                bool * solid_flags,
                bool * solid_old_flags,
                unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            if (idx < size)
            {
                unsigned long i(idx);

                bool fts(!solid_old_flags[i] & solid_flags[i]);
                f_eq_0[i] = fts ? float(0) : f_eq_0[i];
                f_eq_1[i] = fts ? float(0) : f_eq_1[i];
                f_eq_2[i] = fts ? float(0) : f_eq_2[i];
                f_eq_3[i] = fts ? float(0) : f_eq_3[i];
                f_eq_4[i] = fts ? float(0) : f_eq_4[i];
                f_eq_5[i] = fts ? float(0) : f_eq_5[i];
                f_eq_6[i] = fts ? float(0) : f_eq_6[i];
                f_eq_7[i] = fts ? float(0) : f_eq_7[i];
                f_eq_8[i] = fts ? float(0) : f_eq_8[i];
            }
        }
        __global__ void boundary_init_fsi_1_reset_f_temp_gpu(
                float * f_temp_0,
                float * f_temp_1,
                float * f_temp_2,
                float * f_temp_3,
                float * f_temp_4,
                float * f_temp_5,
                float * f_temp_6,
                float * f_temp_7,
                float * f_temp_8,
                bool * solid_flags,
                bool * solid_old_flags,
                unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            if (idx < size)
            {
                unsigned long i(idx);

                bool fts(!solid_old_flags[i] & solid_flags[i]);
                f_temp_0[i] = fts ? float(0) : f_temp_0[i];
                f_temp_1[i] = fts ? float(0) : f_temp_1[i];
                f_temp_2[i] = fts ? float(0) : f_temp_2[i];
                f_temp_3[i] = fts ? float(0) : f_temp_3[i];
                f_temp_4[i] = fts ? float(0) : f_temp_4[i];
                f_temp_5[i] = fts ? float(0) : f_temp_5[i];
                f_temp_6[i] = fts ? float(0) : f_temp_6[i];
                f_temp_7[i] = fts ? float(0) : f_temp_7[i];
                f_temp_8[i] = fts ? float(0) : f_temp_8[i];

                solid_old_flags[i] = solid_flags[i];
            }
        }
    }
}

extern "C" void cuda_boundary_init_fsi_dir_1_float(
        void * dir_5,
        void * boundary_flags ,
        void * solid_flags ,
        void * solid_old_flags ,
        void * h ,
        void * u ,
        void * v ,
        void * f_0 ,
        void * f_1 ,
        void * f_2 ,
        void * f_3 ,
        void * f_4 ,
        void * f_5 ,
        void * f_6 ,
        void * f_7 ,
        void * f_8 ,
        void * f_eq_0 ,
        void * f_eq_1 ,
        void * f_eq_2 ,
        void * f_eq_3 ,
        void * f_eq_4 ,
        void * f_eq_5 ,
        void * f_eq_6 ,
        void * f_eq_7 ,
        void * f_eq_8 ,
        void * f_temp_0 ,
        void * f_temp_1 ,
        void * f_temp_2 ,
        void * f_temp_3 ,
        void * f_temp_4 ,
        void * f_temp_5 ,
        void * f_temp_6 ,
        void * f_temp_7 ,
        void * f_temp_8 ,
        float c_u ,
        float c_v ,
        unsigned long size,
        unsigned long blocksize)
{
    dim3 grid;
    dim3 block;
    block.x = blocksize;
    grid.x = (unsigned)ceil(sqrt(size/(double)block.x));
    grid.y = grid.x;

    unsigned long * dir_5_gpu((unsigned long *)dir_5);

    bool * boundary_flags_gpu((bool *)boundary_flags);
    bool * solid_flags_gpu((bool *)solid_flags);
    bool * solid_old_flags_gpu((bool *)solid_old_flags);

    float * h_gpu((float *)h);
    float * u_gpu((float *)u);
    float * v_gpu((float *)v);

    float * f_0_gpu((float *)f_0);
    float * f_1_gpu((float *)f_1);
    float * f_2_gpu((float *)f_2);
    float * f_3_gpu((float *)f_3);
    float * f_4_gpu((float *)f_4);
    float * f_5_gpu((float *)f_5);
    float * f_6_gpu((float *)f_6);
    float * f_7_gpu((float *)f_7);
    float * f_8_gpu((float *)f_8);
    float * f_eq_0_gpu((float *)f_eq_0);
    float * f_eq_1_gpu((float *)f_eq_1);
    float * f_eq_2_gpu((float *)f_eq_2);
    float * f_eq_3_gpu((float *)f_eq_3);
    float * f_eq_4_gpu((float *)f_eq_4);
    float * f_eq_5_gpu((float *)f_eq_5);
    float * f_eq_6_gpu((float *)f_eq_6);
    float * f_eq_7_gpu((float *)f_eq_7);
    float * f_eq_8_gpu((float *)f_eq_8);
    float * f_temp_0_gpu((float *)f_temp_0);
    float * f_temp_1_gpu((float *)f_temp_1);
    float * f_temp_2_gpu((float *)f_temp_2);
    float * f_temp_3_gpu((float *)f_temp_3);
    float * f_temp_4_gpu((float *)f_temp_4);
    float * f_temp_5_gpu((float *)f_temp_5);
    float * f_temp_6_gpu((float *)f_temp_6);
    float * f_temp_7_gpu((float *)f_temp_7);
    float * f_temp_8_gpu((float *)f_temp_8);



    honei::cuda::boundary_init_fsi_1_set_gpu<<<grid, block>>>(dir_5_gpu,
                                                           boundary_flags_gpu,
                                                           solid_flags_gpu,
                                                           solid_old_flags_gpu,
                                                           h_gpu,
                                                           u_gpu,
                                                           v_gpu,
                                                           c_u,
                                                           c_v,
                                                           size);

    honei::cuda::boundary_init_fsi_1_reset_f_gpu<<<grid, block>>>(f_0_gpu,
                                                                  f_1_gpu,
                                                                  f_2_gpu,
                                                                  f_3_gpu,
                                                                  f_4_gpu,
                                                                  f_5_gpu,
                                                                  f_6_gpu,
                                                                  f_7_gpu,
                                                                  f_8_gpu,
                                                                  solid_flags_gpu,
                                                                  solid_old_flags_gpu,
                                                                  size);

    honei::cuda::boundary_init_fsi_1_reset_f_eq_gpu<<<grid, block>>>(f_eq_0_gpu,
                                                                  f_eq_1_gpu,
                                                                  f_eq_2_gpu,
                                                                  f_eq_3_gpu,
                                                                  f_eq_4_gpu,
                                                                  f_eq_5_gpu,
                                                                  f_eq_6_gpu,
                                                                  f_eq_7_gpu,
                                                                  f_eq_8_gpu,
                                                                  solid_flags_gpu,
                                                                  solid_old_flags_gpu,
                                                                  size);

    honei::cuda::boundary_init_fsi_1_reset_f_temp_gpu<<<grid, block>>>(f_temp_0_gpu,
                                                                  f_temp_1_gpu,
                                                                  f_temp_2_gpu,
                                                                  f_temp_3_gpu,
                                                                  f_temp_4_gpu,
                                                                  f_temp_5_gpu,
                                                                  f_temp_6_gpu,
                                                                  f_temp_7_gpu,
                                                                  f_temp_8_gpu,
                                                                  solid_flags_gpu,
                                                                  solid_old_flags_gpu,
                                                                  size);
    CUDA_ERROR();
}
