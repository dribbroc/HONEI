/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2012 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
        __global__ void collide_stream_grid_pollutant_0_gpu(
                DT_ * h,
                DT_ * f_eq_0,
                DT_ * f_0,
                DT_ * f_temp_0,
                DT_ tau,
                unsigned long offset_0, unsigned long size_0)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            if (idx >= offset_0 && idx < size_0)
            {
                DT_ one_over_2(1./2.);
                DT_ tau_local = one_over_2 + h[idx] * (tau - one_over_2);
                (f_temp_0)[idx] = (f_0)[idx] - ((f_0)[idx] - (f_eq_0)[idx])/tau_local;
            }
        }

        template <typename DT_>
        __global__ void collide_stream_grid_pollutant_n_gpu(
                unsigned long * dir,
                DT_ * h,
                DT_ * f_eq,
                DT_ * f,
                DT_ * f_temp,
                DT_ tau,
                unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;
            if (idx < size)
            {
                if (dir[idx] < size)
                {
                    DT_ one_over_2(1./2.);
                    DT_ tau_local = one_over_2 + h[dir[idx]] * (tau - one_over_2);
                    (f_temp)[dir[idx]] = (f)[idx] - ((f)[idx] - (f_eq)[idx])/tau_local;
                }
            }
        }
    }

    template <typename DT_>
        void cuda_collide_stream_grid_pollutant(unsigned long start, unsigned long end,
                void * h,
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
                DT_ tau, unsigned long size,
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

            DT_ * h_gpu((DT_ *)h);

            DT_ * f_eq_0_gpu((DT_ *)f_eq_0);
            DT_ * f_eq_1_gpu((DT_ *)f_eq_1);
            DT_ * f_eq_2_gpu((DT_ *)f_eq_2);
            DT_ * f_eq_3_gpu((DT_ *)f_eq_3);
            DT_ * f_eq_4_gpu((DT_ *)f_eq_4);
            DT_ * f_eq_5_gpu((DT_ *)f_eq_5);
            DT_ * f_eq_6_gpu((DT_ *)f_eq_6);
            DT_ * f_eq_7_gpu((DT_ *)f_eq_7);
            DT_ * f_eq_8_gpu((DT_ *)f_eq_8);

            DT_ * f_0_gpu((DT_ *)f_0);
            DT_ * f_1_gpu((DT_ *)f_1);
            DT_ * f_2_gpu((DT_ *)f_2);
            DT_ * f_3_gpu((DT_ *)f_3);
            DT_ * f_4_gpu((DT_ *)f_4);
            DT_ * f_5_gpu((DT_ *)f_5);
            DT_ * f_6_gpu((DT_ *)f_6);
            DT_ * f_7_gpu((DT_ *)f_7);
            DT_ * f_8_gpu((DT_ *)f_8);

            DT_ * f_temp_0_gpu((DT_ *)f_temp_0);
            DT_ * f_temp_1_gpu((DT_ *)f_temp_1);
            DT_ * f_temp_2_gpu((DT_ *)f_temp_2);
            DT_ * f_temp_3_gpu((DT_ *)f_temp_3);
            DT_ * f_temp_4_gpu((DT_ *)f_temp_4);
            DT_ * f_temp_5_gpu((DT_ *)f_temp_5);
            DT_ * f_temp_6_gpu((DT_ *)f_temp_6);
            DT_ * f_temp_7_gpu((DT_ *)f_temp_7);
            DT_ * f_temp_8_gpu((DT_ *)f_temp_8);

            honei::cuda::collide_stream_grid_pollutant_0_gpu<DT_><<<grid, block>>>(
                    h_gpu,
                    f_eq_0_gpu,
                    f_0_gpu,
                    f_temp_0_gpu,
                    tau,
                    start, size_0);

            honei::cuda::collide_stream_grid_pollutant_n_gpu<DT_><<<grid, block>>>(dir_1_gpu, h_gpu, f_eq_1_gpu, f_1_gpu, f_temp_1_gpu, tau, size);
            honei::cuda::collide_stream_grid_pollutant_n_gpu<DT_><<<grid, block>>>(dir_2_gpu, h_gpu, f_eq_2_gpu, f_2_gpu, f_temp_2_gpu, tau, size);
            honei::cuda::collide_stream_grid_pollutant_n_gpu<DT_><<<grid, block>>>(dir_3_gpu, h_gpu, f_eq_3_gpu, f_3_gpu, f_temp_3_gpu, tau, size);
            honei::cuda::collide_stream_grid_pollutant_n_gpu<DT_><<<grid, block>>>(dir_4_gpu, h_gpu, f_eq_4_gpu, f_4_gpu, f_temp_4_gpu, tau, size);
            honei::cuda::collide_stream_grid_pollutant_n_gpu<DT_><<<grid, block>>>(dir_5_gpu, h_gpu, f_eq_5_gpu, f_5_gpu, f_temp_5_gpu, tau, size);
            honei::cuda::collide_stream_grid_pollutant_n_gpu<DT_><<<grid, block>>>(dir_6_gpu, h_gpu, f_eq_6_gpu, f_6_gpu, f_temp_6_gpu, tau, size);
            honei::cuda::collide_stream_grid_pollutant_n_gpu<DT_><<<grid, block>>>(dir_7_gpu, h_gpu, f_eq_7_gpu, f_7_gpu, f_temp_7_gpu, tau, size);
            honei::cuda::collide_stream_grid_pollutant_n_gpu<DT_><<<grid, block>>>(dir_8_gpu, h_gpu, f_eq_8_gpu, f_8_gpu, f_temp_8_gpu, tau, size);

            CUDA_ERROR();
        }
}

extern "C" void cuda_collide_stream_grid_pollutant_float(unsigned long start, unsigned long end,
        void * h, void * dir_1, void * dir_2, void * dir_3, void * dir_4,
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
    honei::cuda_collide_stream_grid_pollutant<float>(start, end,
            h, dir_1, dir_2, dir_3, dir_4,
            dir_5, dir_6, dir_7, dir_8,
            f_eq_0, f_eq_1, f_eq_2,
            f_eq_3, f_eq_4, f_eq_5,
            f_eq_6, f_eq_7, f_eq_8,
            f_0, f_1, f_2,
            f_3, f_4, f_5,
            f_6, f_7, f_8,
            f_temp_0, f_temp_1, f_temp_2,
            f_temp_3, f_temp_4, f_temp_5,
            f_temp_6, f_temp_7, f_temp_8,
            tau, size, blocksize);
}

#ifdef HONEI_CUDA_DOUBLE
extern "C" void cuda_collide_stream_grid_pollutant_double(unsigned long start, unsigned long end,
        void * h, void * dir_1, void * dir_2, void * dir_3, void * dir_4,
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
        double tau, unsigned long size,
        unsigned long blocksize)
{
    honei::cuda_collide_stream_grid_pollutant<double>(start, end,
            h, dir_1, dir_2, dir_3, dir_4,
            dir_5, dir_6, dir_7, dir_8,
            f_eq_0, f_eq_1, f_eq_2,
            f_eq_3, f_eq_4, f_eq_5,
            f_eq_6, f_eq_7, f_eq_8,
            f_0, f_1, f_2,
            f_3, f_4, f_5,
            f_6, f_7, f_8,
            f_temp_0, f_temp_1, f_temp_2,
            f_temp_3, f_temp_4, f_temp_5,
            f_temp_6, f_temp_7, f_temp_8,
            tau, size, blocksize);
}
#endif
