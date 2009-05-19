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
                float * distribution_x_data, float * distribution_y_data,
                float * f_eq_0, float * f_eq_1, float * f_eq_2,
                float * f_eq_3, float * f_eq_4, float * f_eq_5,
                float * f_eq_6, float * f_eq_7, float * f_eq_8,
                float g, float e,
                unsigned long offset, unsigned long size)
        {
            extern __shared__ float  distribution_cache[];
            float* distribution_x = distribution_cache;
            float* distribution_y = distribution_cache + 9;

            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;

            if (idx % blockDim.x < 9)
            {
                distribution_x[idx % blockDim.x] = distribution_x_data[idx % blockDim.x];
                distribution_y[idx % blockDim.x] = distribution_y_data[idx % blockDim.x];
            }
            __syncthreads();

            if (idx < size)
            {
                unsigned long index(idx + offset);

                float e2(e);
                float e23(float(3.) * e2);
                float e26(float(6.) * e2);
                float e42(float(2.) * e2 * e2);
                float e48(float(8.) * e2 * e2);
                float e212(float(12.) * e2);
                float e224(float(24.) * e2);

                float u2(u[index] * u[index]);
                float v2(v[index] * v[index]);
                float gh(g * h[index]);

                float dxu, dyv;
                float t1, t2, t3, t4;

                t1 = (float(5.) * gh) / e26;
                t2 = float(2.) / e23 * (u2 + v2);
                f_eq_0[index] = h[index] * (float(1) - t1 - t2);

                dxu = distribution_x[1] * u[index];
                dyv = distribution_y[1] * v[index];
                t1 = (gh) / e26;
                t2 = (dxu + dyv) / e23;
                t3 = (dxu * dxu + float(2.) * dxu * dyv + dyv * dyv) / e42;
                t4 = (u2 + v2) / e26;
                f_eq_1[index] = h[index] * (t1 + t2 + t3 - t4);

                dxu = distribution_x[3] * u[index];
                dyv = distribution_y[3] * v[index];
                t2 = (dxu + dyv) / e23;
                t3 = (dxu * dxu + float(2.) * dxu * dyv + dyv * dyv) / e42;
                f_eq_3[index] = h[index] * (t1 + t2 + t3 - t4);

                dxu = distribution_x[5] * u[index];
                dyv = distribution_y[5] * v[index];
                t2 = (dxu + dyv) / e23;
                t3 = (dxu * dxu + float(2.) * dxu * dyv + dyv * dyv) / e42;
                f_eq_5[index] = h[index] * (t1 + t2 + t3 - t4);

                dxu = distribution_x[7] * u[index];
                dyv = distribution_y[7] * v[index];
                t2 = (dxu + dyv) / e23;
                t3 = (dxu * dxu + float(2.) * dxu * dyv + dyv * dyv) / e42;
                f_eq_7[index] = h[index] * (t1 + t2 + t3 - t4);

                dxu = distribution_x[2] * u[index];
                dyv = distribution_y[2] * v[index];
                t1 = (gh) / e224;
                t2 = (dxu + dyv) / e212;
                t3 = (dxu * dxu + float(2.) * dxu * dyv + dyv * dyv) / e48;
                t4 = (u2 + v2) / e224;
                f_eq_2[index] =  h[index] * (t1 + t2 + t3 - t4);

                dxu = distribution_x[2] * u[index];
                dyv = distribution_y[2] * v[index];
                t2 = (dxu + dyv) / e212;
                t3 = (dxu * dxu + float(2.) * dxu * dyv + dyv * dyv) / e48;
                f_eq_2[index] =  h[index] * (t1 + t2 + t3 - t4);

                dxu = distribution_x[4] * u[index];
                dyv = distribution_y[4] * v[index];
                t2 = (dxu + dyv) / e212;
                t3 = (dxu * dxu + float(2.) * dxu * dyv + dyv * dyv) / e48;
                f_eq_4[index] =  h[index] * (t1 + t2 + t3 - t4);

                dxu = distribution_x[6] * u[index];
                dyv = distribution_y[6] * v[index];
                t2 = (dxu + dyv) / e212;
                t3 = (dxu * dxu + float(2.) * dxu * dyv + dyv * dyv) / e48;
                f_eq_6[index] =  h[index] * (t1 + t2 + t3 - t4);

                dxu = distribution_x[8] * u[index];
                dyv = distribution_y[8] * v[index];
                t2 = (dxu + dyv) / e212;
                t3 = (dxu * dxu + float(2.) * dxu * dyv + dyv * dyv) / e48;
                f_eq_8[index] =  h[index] * (t1 + t2 + t3 - t4);
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

    honei::cuda::eq_dist_grid_gpu<<<grid, block, 18 * sizeof(float)>>>(u_gpu, v_gpu, h_gpu,
            distribution_x_gpu, distribution_y_gpu,
            f_eq_0_gpu, f_eq_1_gpu, f_eq_2_gpu,
            f_eq_3_gpu, f_eq_4_gpu, f_eq_5_gpu,
            f_eq_6_gpu, f_eq_7_gpu, f_eq_8_gpu,
            g, e,
            start, size);

    CUDA_ERROR();
}
