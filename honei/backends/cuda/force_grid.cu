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

#include <limits>
#include <honei/backends/cuda/cuda_util.hh>

namespace honei
{
    namespace cuda
    {
        template <typename DT_>
        __global__ void force_grid_gpu_x(
                unsigned long * dir_1,
                unsigned long * interdir_1,
                DT_ * f_temp_1,
                DT_ * h, DT_ * b,
                DT_ distribution_x,
                DT_ g, DT_ d_x, DT_ d_y, DT_ d_t,
                unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            DT_ force_multiplier(d_t / (6 * d_x * d_x / (d_t * d_t)));
            DT_ gravity_multiplier(-g);
            DT_ force_times_gravity(force_multiplier * gravity_multiplier);

            if (idx < size)
            {
                unsigned long i(idx);
                DT_ x(0);
                if (interdir_1[i] < size)
                {
                    x = force_times_gravity * distribution_x * (h[i] + h[interdir_1[i]]) / DT_(2);
                }
                if (dir_1[i] < size)
                {
                    x *= (b[dir_1[i]] - b[i]) / d_x;
                }
                f_temp_1[i] += x;
            }
        }

        template <typename DT_>
        __global__ void force_grid_gpu_xy(
                unsigned long * dir_2,
                unsigned long * dir_1,
                unsigned long * dir_3,
                DT_ * f_temp_2,
                DT_ * h, DT_ * b,
                DT_ distribution_x, DT_ distribution_y,
                DT_ g, DT_ d_x, DT_ d_y, DT_ d_t,
                unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            DT_ force_multiplier(d_t / (6 * d_x * d_x / (d_t * d_t)));
            DT_ gravity_multiplier(-g);
            DT_ force_times_gravity(force_multiplier * gravity_multiplier);

            if (idx < size)
            {
                unsigned long i(idx);
                DT_ x(0);
                if (dir_2[i] < size)
                {
                    x = force_times_gravity * distribution_x * (h[i] + h[dir_2[i]]) / DT_(2);
                }
                if (dir_1[i] < size)
                {
                    x *= (b[dir_1[i]] - b[i]) / d_x;
                }
                f_temp_2[i] += x;
                DT_ y(0);
                if (dir_2[i] < size)
                {
                    y = force_times_gravity * distribution_y * (h[i] + h[dir_2[i]]) / DT_(2);
                }
                if (dir_3[i] < size)
                {
                    y *= (b[dir_3[i]] - b[i]) / d_y;
                }
                f_temp_2[i] += y;
            }
        }

        template <typename DT_>
        __global__ void force_grid_gpu_y(
                unsigned long * dir_3,
                unsigned long * interdir_3,
                DT_ * f_temp_3,
                DT_ * h, DT_ * b,
                DT_ distribution_y,
                DT_ g, DT_ d_x, DT_ d_y, DT_ d_t,
                unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            DT_ force_multiplier(d_t / (6 * d_x * d_x / (d_t * d_t)));
            DT_ gravity_multiplier(-g);
            DT_ force_times_gravity(force_multiplier * gravity_multiplier);

            if (idx < size)
            {
                unsigned long i(idx);
                DT_ y(0);
                if (interdir_3[i] < size)
                {
                    y = force_times_gravity * distribution_y * (h[i] + h[interdir_3[i]]) / DT_(2);
                }
                if (dir_3[i] < size)
                {
                    y *= (b[dir_3[i]] - b[i]) / d_y;
                }
                f_temp_3[i] += y;
            }
        }
    }

    template <typename DT_>
    void cuda_force_grid(
            void * dir_1, void * dir_2, void * dir_3, void * dir_4,
            void * dir_5, void * dir_6, void * dir_7, void * dir_8,
            void * h, void * b,
            void * f_temp_1, void * f_temp_2,
            void * f_temp_3, void * f_temp_4, void * f_temp_5,
            void * f_temp_6, void * f_temp_7, void * f_temp_8,
            DT_ distribution_x_1, DT_ distribution_y_1,
            DT_ distribution_x_2, DT_ distribution_y_2,
            DT_ distribution_x_3, DT_ distribution_y_3,
            DT_ distribution_x_4, DT_ distribution_y_4,
            DT_ distribution_x_5, DT_ distribution_y_5,
            DT_ distribution_x_6, DT_ distribution_y_6,
            DT_ distribution_x_7, DT_ distribution_y_7,
            DT_ distribution_x_8, DT_ distribution_y_8,
            DT_ g, DT_ d_x, DT_ d_y, DT_ d_t,
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

        DT_ * h_gpu((DT_ *)h);
        DT_ * b_gpu((DT_ *)b);

        DT_ * f_temp_1_gpu((DT_ *)f_temp_1);
        DT_ * f_temp_2_gpu((DT_ *)f_temp_2);
        DT_ * f_temp_3_gpu((DT_ *)f_temp_3);
        DT_ * f_temp_4_gpu((DT_ *)f_temp_4);
        DT_ * f_temp_5_gpu((DT_ *)f_temp_5);
        DT_ * f_temp_6_gpu((DT_ *)f_temp_6);
        DT_ * f_temp_7_gpu((DT_ *)f_temp_7);
        DT_ * f_temp_8_gpu((DT_ *)f_temp_8);

        honei::cuda::force_grid_gpu_x<DT_><<<grid, block>>>(
                dir_1_gpu, dir_1_gpu, f_temp_1_gpu,
                h_gpu, b_gpu,
                distribution_x_1,
                g, d_x, d_y, d_t,
                size);
        honei::cuda::force_grid_gpu_xy<DT_><<<grid, block>>>(
                dir_2_gpu, dir_1_gpu, dir_3_gpu, f_temp_2_gpu,
                h_gpu, b_gpu,
                distribution_x_2, distribution_y_2,
                g, d_x, d_y, d_t,
                size);
        honei::cuda::force_grid_gpu_y<DT_><<<grid, block>>>(
                dir_3_gpu, dir_3_gpu, f_temp_3_gpu,
                h_gpu, b_gpu,
                distribution_y_3,
                g, d_x, d_y, d_t,
                size);
        honei::cuda::force_grid_gpu_xy<DT_><<<grid, block>>>(
                dir_4_gpu, dir_1_gpu, dir_3_gpu, f_temp_4_gpu,
                h_gpu, b_gpu,
                distribution_x_4, distribution_y_4,
                g, d_x, d_y, d_t,
                size);
        honei::cuda::force_grid_gpu_x<DT_><<<grid, block>>>(
                dir_1_gpu, dir_5_gpu, f_temp_5_gpu,
                h_gpu, b_gpu,
                distribution_x_5,
                g, d_x, d_y, d_t,
                size);
        honei::cuda::force_grid_gpu_xy<DT_><<<grid, block>>>(
                dir_6_gpu, dir_1_gpu, dir_3_gpu, f_temp_6_gpu,
                h_gpu, b_gpu,
                distribution_x_6, distribution_y_6,
                g, d_x, d_y, d_t,
                size);
        honei::cuda::force_grid_gpu_y<DT_><<<grid, block>>>(
                dir_3_gpu, dir_7_gpu, f_temp_7_gpu,
                h_gpu, b_gpu,
                distribution_y_7,
                g, d_x, d_y, d_t,
                size);
        honei::cuda::force_grid_gpu_xy<DT_><<<grid, block>>>(
                dir_8_gpu, dir_1_gpu, dir_3_gpu, f_temp_8_gpu,
                h_gpu, b_gpu,
                distribution_x_8, distribution_y_8,
                g, d_x, d_y, d_t,
                size);

        CUDA_ERROR();
    }
}

extern "C" void cuda_force_grid_float(
            void * dir_1, void * dir_2, void * dir_3, void * dir_4,
            void * dir_5, void * dir_6, void * dir_7, void * dir_8,
            void * h, void * b,
            void * f_temp_1, void * f_temp_2,
            void * f_temp_3, void * f_temp_4, void * f_temp_5,
            void * f_temp_6, void * f_temp_7, void * f_temp_8,
            float distribution_x_1, float distribution_y_1,
            float distribution_x_2, float distribution_y_2,
            float distribution_x_3, float distribution_y_3,
            float distribution_x_4, float distribution_y_4,
            float distribution_x_5, float distribution_y_5,
            float distribution_x_6, float distribution_y_6,
            float distribution_x_7, float distribution_y_7,
            float distribution_x_8, float distribution_y_8,
            float g, float d_x, float d_y, float d_t,
            unsigned long size,
            unsigned long blocksize)
{
    honei::cuda_force_grid<float>(
            dir_1, dir_2, dir_3, dir_4,
            dir_5, dir_6, dir_7, dir_8,
            h, b,
            f_temp_1, f_temp_2,
            f_temp_3, f_temp_4, f_temp_5,
            f_temp_6, f_temp_7, f_temp_8,
            distribution_x_1, distribution_y_1,
            distribution_x_2, distribution_y_2,
            distribution_x_3, distribution_y_3,
            distribution_x_4, distribution_y_4,
            distribution_x_5, distribution_y_5,
            distribution_x_6, distribution_y_6,
            distribution_x_7, distribution_y_7,
            distribution_x_8, distribution_y_8,
            g, d_x, d_y, d_t,
            size, blocksize);
}

#ifdef HONEI_CUDA_DOUBLE
extern "C" void cuda_force_grid_double(
        void * dir_1, void * dir_2, void * dir_3, void * dir_4,
        void * dir_5, void * dir_6, void * dir_7, void * dir_8,
        void * h, void * b,
        void * f_temp_1, void * f_temp_2,
        void * f_temp_3, void * f_temp_4, void * f_temp_5,
        void * f_temp_6, void * f_temp_7, void * f_temp_8,
        double distribution_x_1, double distribution_y_1,
        double distribution_x_2, double distribution_y_2,
        double distribution_x_3, double distribution_y_3,
        double distribution_x_4, double distribution_y_4,
        double distribution_x_5, double distribution_y_5,
        double distribution_x_6, double distribution_y_6,
        double distribution_x_7, double distribution_y_7,
        double distribution_x_8, double distribution_y_8,
        double g, double d_x, double d_y, double d_t,
        unsigned long size,
        unsigned long blocksize)
{
    honei::cuda_force_grid<double>(
            dir_1, dir_2, dir_3, dir_4,
            dir_5, dir_6, dir_7, dir_8,
            h, b,
            f_temp_1, f_temp_2,
            f_temp_3, f_temp_4, f_temp_5,
            f_temp_6, f_temp_7, f_temp_8,
            distribution_x_1, distribution_y_1,
            distribution_x_2, distribution_y_2,
            distribution_x_3, distribution_y_3,
            distribution_x_4, distribution_y_4,
            distribution_x_5, distribution_y_5,
            distribution_x_6, distribution_y_6,
            distribution_x_7, distribution_y_7,
            distribution_x_8, distribution_y_8,
            g, d_x, d_y, d_t,
            size, blocksize);
}
#endif

//-----------------------------------------------------------------------------------------------
namespace honei
{
    namespace cuda
    {
        template <typename DT_>
        __global__ void force_grid_gpu_x_2(
                DT_ * f_temp,
                DT_ * h, DT_ * u, DT_ * v,
                DT_ distribution_x,
                DT_ g, DT_ d_x, DT_ d_t, DT_ manning, DT_ epsilon,
                unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;

            if (idx < size)
            {
                unsigned long i(idx);
                DT_ force_multiplier(d_t / (DT_(6) * d_x * d_x / (d_t * d_t)));
                force_multiplier *= g;
                if ( (pow(h[i], DT_(1./3.))) > epsilon || (pow(h[i], DT_(1./3.))) < -epsilon )
                {
                    f_temp[i] -= force_multiplier * distribution_x * manning * manning *
                        u[i] * sqrt(u[i] * u[i] + v[i] * v[i]) / (pow(h[i], DT_(1./3.)));
                }
            }
        }

        template <typename DT_>
        __global__ void force_grid_gpu_xy_2(
                DT_ * f_temp,
                DT_ * h, DT_ * u, DT_ * v,
                DT_ distribution_x,
                DT_ distribution_y,
                DT_ g, DT_ d_x, DT_ d_t, DT_ manning, DT_ epsilon,
                unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;

            if (idx < size)
            {
                unsigned long i(idx);
                DT_ force_multiplier(d_t / (DT_(6) * d_x * d_x / (d_t * d_t)));
                force_multiplier *= g;
                if ( (pow(h[i], DT_(1./3.))) > epsilon || (pow(h[i], DT_(1./3.))) < -epsilon )
                {
                    f_temp[i] -= force_multiplier * distribution_x * manning * manning *
                        u[i] * sqrt(u[i] * u[i] + v[i] * v[i]) / (pow(h[i], DT_(1./3.)));
                    f_temp[i] -= force_multiplier * distribution_y * manning * manning *
                        u[i] * sqrt(u[i] * u[i] + v[i] * v[i]) / (pow(h[i], DT_(1./3.)));
                }
            }
        }
    }

    template <typename DT_>
    void cuda_force_grid_2(
            void * h, void * u, void * v,
            void * f_temp_1, void * f_temp_2,
            void * f_temp_3, void * f_temp_4, void * f_temp_5,
            void * f_temp_6, void * f_temp_7, void * f_temp_8,
            DT_ distribution_x_1, DT_ distribution_y_1,
            DT_ distribution_x_2, DT_ distribution_y_2,
            DT_ distribution_x_3, DT_ distribution_y_3,
            DT_ distribution_x_4, DT_ distribution_y_4,
            DT_ distribution_x_5, DT_ distribution_y_5,
            DT_ distribution_x_6, DT_ distribution_y_6,
            DT_ distribution_x_7, DT_ distribution_y_7,
            DT_ distribution_x_8, DT_ distribution_y_8,
            DT_ g, DT_ d_x, DT_ d_y, DT_ d_t, DT_ manning,
            unsigned long size,
            unsigned long blocksize)
    {
        dim3 grid;
        dim3 block;
        block.x = blocksize;
        grid.x = (unsigned)ceil(sqrt(size/(double)block.x));
        grid.y = grid.x;

        DT_ * h_gpu((DT_ *)h);
        DT_ * u_gpu((DT_ *)u);
        DT_ * v_gpu((DT_ *)v);

        DT_ * f_temp_1_gpu((DT_ *)f_temp_1);
        DT_ * f_temp_2_gpu((DT_ *)f_temp_2);
        DT_ * f_temp_3_gpu((DT_ *)f_temp_3);
        DT_ * f_temp_4_gpu((DT_ *)f_temp_4);
        DT_ * f_temp_5_gpu((DT_ *)f_temp_5);
        DT_ * f_temp_6_gpu((DT_ *)f_temp_6);
        DT_ * f_temp_7_gpu((DT_ *)f_temp_7);
        DT_ * f_temp_8_gpu((DT_ *)f_temp_8);

        honei::cuda::force_grid_gpu_x_2<DT_><<<grid, block>>>(
                f_temp_1_gpu,
                h_gpu, u_gpu, v_gpu,
                distribution_x_1,
                g, d_x, d_t, manning, std::numeric_limits<DT_>::epsilon(),
                size);
        honei::cuda::force_grid_gpu_xy_2<DT_><<<grid, block>>>(
                f_temp_2_gpu,
                h_gpu, u_gpu, v_gpu,
                distribution_x_2,
                distribution_y_2,
                g, d_x, d_t, manning, std::numeric_limits<DT_>::epsilon(),
                size);
        honei::cuda::force_grid_gpu_x_2<DT_><<<grid, block>>>(
                f_temp_3_gpu,
                h_gpu, u_gpu, v_gpu,
                distribution_y_3,
                g, d_x, d_t, manning, std::numeric_limits<DT_>::epsilon(),
                size);
        honei::cuda::force_grid_gpu_xy_2<DT_><<<grid, block>>>(
                f_temp_4_gpu,
                h_gpu, u_gpu, v_gpu,
                distribution_x_4,
                distribution_y_4,
                g, d_x, d_t, manning, std::numeric_limits<DT_>::epsilon(),
                size);
        honei::cuda::force_grid_gpu_x_2<DT_><<<grid, block>>>(
                f_temp_5_gpu,
                h_gpu, u_gpu, v_gpu,
                distribution_x_5,
                g, d_x, d_t, manning, std::numeric_limits<DT_>::epsilon(),
                size);
        honei::cuda::force_grid_gpu_xy_2<DT_><<<grid, block>>>(
                f_temp_6_gpu,
                h_gpu, u_gpu, v_gpu,
                distribution_x_6,
                distribution_y_6,
                g, d_x, d_t, manning, std::numeric_limits<DT_>::epsilon(),
                size);
        honei::cuda::force_grid_gpu_x_2<DT_><<<grid, block>>>(
                f_temp_7_gpu,
                h_gpu, u_gpu, v_gpu,
                distribution_y_7,
                g, d_x, d_t, manning, std::numeric_limits<DT_>::epsilon(),
                size);
        honei::cuda::force_grid_gpu_xy_2<DT_><<<grid, block>>>(
                f_temp_8_gpu,
                h_gpu, u_gpu, v_gpu,
                distribution_x_8,
                distribution_y_8,
                g, d_x, d_t, manning, std::numeric_limits<DT_>::epsilon(),
                size);

        CUDA_ERROR();
    }
}

extern "C" void cuda_force_grid_float_2(
        void * h, void * u, void * v,
        void * f_temp_1, void * f_temp_2,
        void * f_temp_3, void * f_temp_4, void * f_temp_5,
        void * f_temp_6, void * f_temp_7, void * f_temp_8,
        float distribution_x_1, float distribution_y_1,
        float distribution_x_2, float distribution_y_2,
        float distribution_x_3, float distribution_y_3,
        float distribution_x_4, float distribution_y_4,
        float distribution_x_5, float distribution_y_5,
        float distribution_x_6, float distribution_y_6,
        float distribution_x_7, float distribution_y_7,
        float distribution_x_8, float distribution_y_8,
        float g, float d_x, float d_y, float d_t, float manning,
        unsigned long size,
        unsigned long blocksize)
{
    honei::cuda_force_grid_2<float>(
            h, u, v,
            f_temp_1, f_temp_2,
            f_temp_3, f_temp_4, f_temp_5,
            f_temp_6, f_temp_7, f_temp_8,
            distribution_x_1, distribution_y_1,
            distribution_x_2, distribution_y_2,
            distribution_x_3, distribution_y_3,
            distribution_x_4, distribution_y_4,
            distribution_x_5, distribution_y_5,
            distribution_x_6, distribution_y_6,
            distribution_x_7, distribution_y_7,
            distribution_x_8, distribution_y_8,
            g, d_x, d_y, d_t, manning,
            size, blocksize);
}

#ifdef HONEI_CUDA_DOUBLE
extern "C" void cuda_force_grid_double_2(
        void * h, void * u, void * v,
        void * f_temp_1, void * f_temp_2,
        void * f_temp_3, void * f_temp_4, void * f_temp_5,
        void * f_temp_6, void * f_temp_7, void * f_temp_8,
        double distribution_x_1, double distribution_y_1,
        double distribution_x_2, double distribution_y_2,
        double distribution_x_3, double distribution_y_3,
        double distribution_x_4, double distribution_y_4,
        double distribution_x_5, double distribution_y_5,
        double distribution_x_6, double distribution_y_6,
        double distribution_x_7, double distribution_y_7,
        double distribution_x_8, double distribution_y_8,
        double g, double d_x, double d_y, double d_t, double manning,
        unsigned long size,
        unsigned long blocksize)
{
    honei::cuda_force_grid_2<double>(
            h, u, v,
            f_temp_1, f_temp_2,
            f_temp_3, f_temp_4, f_temp_5,
            f_temp_6, f_temp_7, f_temp_8,
            distribution_x_1, distribution_y_1,
            distribution_x_2, distribution_y_2,
            distribution_x_3, distribution_y_3,
            distribution_x_4, distribution_y_4,
            distribution_x_5, distribution_y_5,
            distribution_x_6, distribution_y_6,
            distribution_x_7, distribution_y_7,
            distribution_x_8, distribution_y_8,
            g, d_x, d_y, d_t, manning,
            size, blocksize);
}
#endif
