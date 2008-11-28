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
        __global__ void eq_dist_grid_gpu(float * u, float * v, float * h,
                float * distribution_x, float * distribution_y,
                float * f_eq_0, float * f_eq_1, float * f_eq_2,
                float * f_eq_3, float * f_eq_4, float * f_eq_5,
                float * f_eq_6, float * f_eq_7, float * f_eq_8,
                float g, float e,
                unsigned long offset, unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            if (idx < size)
            {
                unsigned long index(idx + offset);
                float us(u[index]);
                float vs(v[index]);
                float hs(h[index]);
                f_eq_0[index] = hs -
                    ((float(5.) * g * hs * hs) / (float(6.) * e * e)) -
                    ((float(2.) * hs) /(float(3.) * e * e) * (us * us + vs * vs));

                f_eq_1[index] = ((g * hs * hs) /(float(6.) * e * e)) +
                    ((hs / (float(3.) * e * e)) * (distribution_x[1] * us + distribution_y[1] * vs)) +
                    ((hs / (float(2.) * e * e)) * (distribution_x[1] * us * distribution_x[1] * us + float(2.) * distribution_x[1] * us * distribution_y[1] * vs + distribution_y[1] * vs * distribution_y[1] * vs)) -
                    ((hs / (float(6.) * e * e)) * (us * us + vs * vs));

                f_eq_3[index] = ((g * hs * hs) /(float(6.) * e * e)) +
                    ((hs / (float(3.) * e * e)) * (distribution_x[3] * us + distribution_y[3] * vs)) +
                    ((hs / (float(2.) * e * e)) * (distribution_x[3] * us * distribution_x[3] * us + float(2.) * distribution_x[3] * us * distribution_y[3] * vs + distribution_y[3] * vs * distribution_y[3] * vs)) -
                    ((hs / (float(6.) * e * e)) * (us * us + vs * vs));

                f_eq_5[index] = ((g * hs * hs) /(float(6.) * e * e)) +
                    ((hs / (float(3.) * e * e)) * (distribution_x[5] * us + distribution_y[5] * vs)) +
                    ((hs / (float(2.) * e * e)) * (distribution_x[5] * us * distribution_x[5] * us + float(2.) * distribution_x[5] * us * distribution_y[5] * vs + distribution_y[5] * vs * distribution_y[5] * vs)) -
                    ((hs / (float(6.) * e * e)) * (us * us + vs * vs));

                f_eq_7[index] = ((g * hs * hs) /(float(6.) * e * e)) +
                    ((hs / (float(3.) * e * e)) * (distribution_x[7] * us + distribution_y[7] * vs)) +
                    ((hs / (float(2.) * e * e)) * (distribution_x[7] * us * distribution_x[7] * us + float(2.) * distribution_x[7] * us * distribution_y[7] * vs + distribution_y[7] * vs * distribution_y[7] * vs)) -
                    ((hs / (float(6.) * e * e)) * (us * us + vs * vs));

                f_eq_2[index] = ((g * hs * hs) /(float(24.) * e * e)) +
                    ((hs / (float(12.) * e * e)) * (distribution_x[2] * us + distribution_y[2] * vs)) +
                    ((hs / (float(8.) * e * e)) * (distribution_x[2] * us * distribution_x[2] * us + float(2.) * distribution_x[2] * us * distribution_y[2] * vs + distribution_y[2] * vs * distribution_y[2] * vs)) -
                    ((hs / (float(24.) * e * e)) * (us * us + vs * vs));

                f_eq_4[index] = ((g * hs * hs) /(float(24.) * e * e)) +
                    ((hs / (float(12.) * e * e)) * (distribution_x[4] * us + distribution_y[4] * vs)) +
                    ((hs / (float(8.) * e * e)) * (distribution_x[4] * us * distribution_x[4] * us + float(2.) * distribution_x[4] * us * distribution_y[4] * vs + distribution_y[4] * vs * distribution_y[4] * vs)) -
                    ((hs / (float(24.) * e * e)) * (us * us + vs * vs));

                f_eq_6[index] = ((g * hs * hs) /(float(24.) * e * e)) +
                    ((hs / (float(12.) * e * e)) * (distribution_x[6] * us + distribution_y[6] * vs)) +
                    ((hs / (float(8.) * e * e)) * (distribution_x[6] * us * distribution_x[6] * us + float(2.) * distribution_x[6] * us * distribution_y[6] * vs + distribution_y[6] * vs * distribution_y[6] * vs)) -
                    ((hs / (float(24.) * e * e)) * (us * us + vs * vs));

                f_eq_8[index] = ((g * hs * hs) /(float(24.) * e * e)) +
                    ((hs / (float(12.) * e * e)) * (distribution_x[8] * us + distribution_y[8] * vs)) +
                    ((hs / (float(8.) * e * e)) * (distribution_x[8] * us * distribution_x[8] * us + float(2.) * distribution_x[8] * us * distribution_y[8] * vs + distribution_y[8] * vs * distribution_y[8] * vs)) -
                    ((hs / (float(24.) * e * e)) * (us * us + vs * vs));
            }
        }
    }
}

extern "C" void cuda_eq_dist_grid_float(unsigned long start, unsigned long end, void * u, void * v, void * h,
        void * distribution_x, void * distribution_y,
        void * f_eq_0, void * f_eq_1, void * f_eq_2,
        void * f_eq_3, void * f_eq_4, void * f_eq_5,
        void * f_eq_6, void * f_eq_7, void * f_eq_8,
        float g, float e,
        unsigned long blocksize)
{
    unsigned long size(end - start);
    dim3 grid;
    dim3 block;
    block.x = blocksize;
    grid.x = (unsigned)ceil(sqrt(size/(double)block.x));
    grid.y = grid.x;
    float * u_gpu((float *)u);
    float * v_gpu((float *)v);
    float * h_gpu((float *)h);
    float * distribution_x_gpu((float *)distribution_x);
    float * distribution_y_gpu((float *)distribution_y);
    float * f_eq_0_gpu((float *)f_eq_0);
    float * f_eq_1_gpu((float *)f_eq_1);
    float * f_eq_2_gpu((float *)f_eq_2);
    float * f_eq_3_gpu((float *)f_eq_3);
    float * f_eq_4_gpu((float *)f_eq_4);
    float * f_eq_5_gpu((float *)f_eq_5);
    float * f_eq_6_gpu((float *)f_eq_6);
    float * f_eq_7_gpu((float *)f_eq_7);
    float * f_eq_8_gpu((float *)f_eq_8);

    honei::cuda::eq_dist_grid_gpu<<<grid, block>>>(u_gpu, v_gpu, h_gpu,
            distribution_x_gpu, distribution_y_gpu,
            f_eq_0_gpu, f_eq_1_gpu, f_eq_2_gpu,
            f_eq_3_gpu, f_eq_4_gpu, f_eq_5_gpu,
            f_eq_6_gpu, f_eq_7_gpu, f_eq_8_gpu,
            g, e,
            start, size);

    CUDA_ERROR();
}
