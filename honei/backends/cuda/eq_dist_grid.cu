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
        template <typename DT_>
        class SharedMem
        {
        };

        template <>
        class SharedMem<float>
        {
            public:
                __device__ float * get_pointer(){extern __shared__ float s_float[]; return s_float;}
        };

        template <>
        class SharedMem<double>
        {
            public:
                __device__ double * get_pointer(){extern __shared__ double s_double[]; return s_double;}
        };

        template <typename DT_>
        __global__ void eq_dist_grid_gpu(DT_ * u, DT_ * v, DT_ * h,
                DT_ * distribution_x_data, DT_ * distribution_y_data,
                DT_ * f_eq_0, DT_ * f_eq_1, DT_ * f_eq_2,
                DT_ * f_eq_3, DT_ * f_eq_4, DT_ * f_eq_5,
                DT_ * f_eq_6, DT_ * f_eq_7, DT_ * f_eq_8,
                DT_ g, DT_ e,
                unsigned long offset, unsigned long size)
        {
            //extern __shared__ DT_  distribution_cache[];
            SharedMem<DT_> shared;
            DT_ * distribution_cache = shared.get_pointer();

            DT_* distribution_x = distribution_cache;
            DT_* distribution_y = distribution_cache + 9;

            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;

            if (idx % blockDim.x < 9)
            {
                distribution_x[idx % blockDim.x] = distribution_x_data[idx % blockDim.x];
                distribution_y[idx % blockDim.x] = distribution_y_data[idx % blockDim.x];
            }
            __syncthreads();

            if (idx >= offset && idx < size)
            {
                unsigned long index(idx);

                DT_ e2(e);
                DT_ e23(DT_(3.) * e2);
                DT_ e26(DT_(6.) * e2);
                DT_ e42(DT_(2.) * e2 * e2);
                DT_ e48(DT_(8.) * e2 * e2);
                DT_ e212(DT_(12.) * e2);
                DT_ e224(DT_(24.) * e2);

                DT_ u2(u[index] * u[index]);
                DT_ v2(v[index] * v[index]);
                DT_ gh(g * h[index]);

                DT_ dxu, dyv;
                DT_ t1, t2, t3, t4;

                t1 = (DT_(5.) * gh) / e26;
                t2 = DT_(2.) / e23 * (u2 + v2);
                f_eq_0[index] = h[index] * (DT_(1) - t1 - t2);

                dxu = distribution_x[1] * u[index];
                dyv = distribution_y[1] * v[index];
                t1 = (gh) / e26;
                t2 = (dxu + dyv) / e23;
                t3 = (dxu * dxu + DT_(2.) * dxu * dyv + dyv * dyv) / e42;
                t4 = (u2 + v2) / e26;
                f_eq_1[index] = h[index] * (t1 + t2 + t3 - t4);

                dxu = distribution_x[3] * u[index];
                dyv = distribution_y[3] * v[index];
                t2 = (dxu + dyv) / e23;
                t3 = (dxu * dxu + DT_(2.) * dxu * dyv + dyv * dyv) / e42;
                f_eq_3[index] = h[index] * (t1 + t2 + t3 - t4);

                dxu = distribution_x[5] * u[index];
                dyv = distribution_y[5] * v[index];
                t2 = (dxu + dyv) / e23;
                t3 = (dxu * dxu + DT_(2.) * dxu * dyv + dyv * dyv) / e42;
                f_eq_5[index] = h[index] * (t1 + t2 + t3 - t4);

                dxu = distribution_x[7] * u[index];
                dyv = distribution_y[7] * v[index];
                t2 = (dxu + dyv) / e23;
                t3 = (dxu * dxu + DT_(2.) * dxu * dyv + dyv * dyv) / e42;
                f_eq_7[index] = h[index] * (t1 + t2 + t3 - t4);

                dxu = distribution_x[2] * u[index];
                dyv = distribution_y[2] * v[index];
                t1 = (gh) / e224;
                t2 = (dxu + dyv) / e212;
                t3 = (dxu * dxu + DT_(2.) * dxu * dyv + dyv * dyv) / e48;
                t4 = (u2 + v2) / e224;
                f_eq_2[index] =  h[index] * (t1 + t2 + t3 - t4);

                dxu = distribution_x[2] * u[index];
                dyv = distribution_y[2] * v[index];
                t2 = (dxu + dyv) / e212;
                t3 = (dxu * dxu + DT_(2.) * dxu * dyv + dyv * dyv) / e48;
                f_eq_2[index] =  h[index] * (t1 + t2 + t3 - t4);

                dxu = distribution_x[4] * u[index];
                dyv = distribution_y[4] * v[index];
                t2 = (dxu + dyv) / e212;
                t3 = (dxu * dxu + DT_(2.) * dxu * dyv + dyv * dyv) / e48;
                f_eq_4[index] =  h[index] * (t1 + t2 + t3 - t4);

                dxu = distribution_x[6] * u[index];
                dyv = distribution_y[6] * v[index];
                t2 = (dxu + dyv) / e212;
                t3 = (dxu * dxu + DT_(2.) * dxu * dyv + dyv * dyv) / e48;
                f_eq_6[index] =  h[index] * (t1 + t2 + t3 - t4);

                dxu = distribution_x[8] * u[index];
                dyv = distribution_y[8] * v[index];
                t2 = (dxu + dyv) / e212;
                t3 = (dxu * dxu + DT_(2.) * dxu * dyv + dyv * dyv) / e48;
                f_eq_8[index] =  h[index] * (t1 + t2 + t3 - t4);
            }
        }
    }

    template <typename DT_>
    void cuda_eq_dist_grid(unsigned long start, unsigned long end, void * u, void * v, void * h,
            void * distribution_x, void * distribution_y,
            void * f_eq_0, void * f_eq_1, void * f_eq_2,
            void * f_eq_3, void * f_eq_4, void * f_eq_5,
            void * f_eq_6, void * f_eq_7, void * f_eq_8,
            DT_ g, DT_ e,
            unsigned long blocksize)
    {
        unsigned long size(end);
        dim3 grid;
        dim3 block;
        block.x = blocksize;
        grid.x = (unsigned)ceil(sqrt(size/(double)block.x));
        grid.y = grid.x;
        DT_ * u_gpu((DT_ *)u);
        DT_ * v_gpu((DT_ *)v);
        DT_ * h_gpu((DT_ *)h);
        DT_ * distribution_x_gpu((DT_ *)distribution_x);
        DT_ * distribution_y_gpu((DT_ *)distribution_y);
        DT_ * f_eq_0_gpu((DT_ *)f_eq_0);
        DT_ * f_eq_1_gpu((DT_ *)f_eq_1);
        DT_ * f_eq_2_gpu((DT_ *)f_eq_2);
        DT_ * f_eq_3_gpu((DT_ *)f_eq_3);
        DT_ * f_eq_4_gpu((DT_ *)f_eq_4);
        DT_ * f_eq_5_gpu((DT_ *)f_eq_5);
        DT_ * f_eq_6_gpu((DT_ *)f_eq_6);
        DT_ * f_eq_7_gpu((DT_ *)f_eq_7);
        DT_ * f_eq_8_gpu((DT_ *)f_eq_8);

        honei::cuda::eq_dist_grid_gpu<DT_><<<grid, block, 18 * sizeof(DT_)>>>(u_gpu, v_gpu, h_gpu,
                distribution_x_gpu, distribution_y_gpu,
                f_eq_0_gpu, f_eq_1_gpu, f_eq_2_gpu,
                f_eq_3_gpu, f_eq_4_gpu, f_eq_5_gpu,
                f_eq_6_gpu, f_eq_7_gpu, f_eq_8_gpu,
                g, e,
                start, size);

        CUDA_ERROR();
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
    honei::cuda_eq_dist_grid<float>(start, end, u, v, h, distribution_x, distribution_y, f_eq_0, f_eq_1, f_eq_2, f_eq_3, f_eq_4, f_eq_5,
            f_eq_6, f_eq_7, f_eq_8, g, e, blocksize);
}

#ifdef HONEI_CUDA_DOUBLE
extern "C" void cuda_eq_dist_grid_double(unsigned long start, unsigned long end, void * u, void * v, void * h,
        void * distribution_x, void * distribution_y,
        void * f_eq_0, void * f_eq_1, void * f_eq_2,
        void * f_eq_3, void * f_eq_4, void * f_eq_5,
        void * f_eq_6, void * f_eq_7, void * f_eq_8,
        double g, double e,
        unsigned long blocksize)
{
    honei::cuda_eq_dist_grid<double>(start, end, u, v, h, distribution_x, distribution_y, f_eq_0, f_eq_1, f_eq_2, f_eq_3, f_eq_4, f_eq_5,
            f_eq_6, f_eq_7, f_eq_8, g, e, blocksize);
}
#endif
