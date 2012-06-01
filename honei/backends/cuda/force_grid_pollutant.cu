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
        __global__ void force_grid_pollutant_gpu(
                DT_ * h, DT_ * c,
                DT_ * f_temp_0, DT_ * f_temp_1, DT_ * f_temp_2,
                DT_ * f_temp_3, DT_ * f_temp_4, DT_ * f_temp_5,
                DT_ * f_temp_6, DT_ * f_temp_7, DT_ * f_temp_8,
                DT_ dt, DT_ k, DT_ s_0,
                unsigned long offset, unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;

            if (idx >= offset && idx < size)
            {
                DT_ t1, t2;
                DT_ w_0(4./9.);
                DT_ w_odd(1./9.);
                DT_ w_even(1./36.);
                DT_ const_multiplier_0(dt * w_0);
                DT_ const_multiplier_odd(dt * w_odd);
                DT_ const_multiplier_even(dt * w_even);

                t1 = -k * h[idx] * c[idx];
                t2 = s_0 * h[idx];
                t1 = t1 + t2;

                (f_temp_0)[idx] += const_multiplier_0 * (t1);

                t1 = const_multiplier_odd * (t1);
                t2 = const_multiplier_even * (t2);

                (f_temp_1)[idx] += t1;
                (f_temp_3)[idx] += t1;
                (f_temp_5)[idx] += t1;
                (f_temp_7)[idx] += t1;
                (f_temp_2)[idx] += t2;
                (f_temp_4)[idx] += t2;
                (f_temp_6)[idx] += t2;
                (f_temp_8)[idx] += t2;
            }
        }
    }

    template <typename DT_>
        void cuda_force_grid_pollutant(
                unsigned long start, unsigned long end,
                void * h, void * c,
                void * f_temp_0, void * f_temp_1, void * f_temp_2,
                void * f_temp_3, void * f_temp_4, void * f_temp_5,
                void * f_temp_6, void * f_temp_7, void * f_temp_8,
                DT_ dt, DT_ k, DT_ s_0,
                unsigned long blocksize)
        {
            unsigned long size(end);
            dim3 grid;
            dim3 block;
            block.x = blocksize;
            grid.x = (unsigned)ceil(sqrt(size/(double)block.x));
            grid.y = grid.x;

            DT_ * f_temp_0_gpu((DT_ *)f_temp_0);
            DT_ * f_temp_1_gpu((DT_ *)f_temp_1);
            DT_ * f_temp_2_gpu((DT_ *)f_temp_2);
            DT_ * f_temp_3_gpu((DT_ *)f_temp_3);
            DT_ * f_temp_4_gpu((DT_ *)f_temp_4);
            DT_ * f_temp_5_gpu((DT_ *)f_temp_5);
            DT_ * f_temp_6_gpu((DT_ *)f_temp_6);
            DT_ * f_temp_7_gpu((DT_ *)f_temp_7);
            DT_ * f_temp_8_gpu((DT_ *)f_temp_8);

            DT_ * h_gpu((DT_ *)h);
            DT_ * c_gpu((DT_ *)c);

            honei::cuda::force_grid_pollutant_gpu<DT_><<<grid, block>>>(h_gpu, c_gpu,
                    f_temp_0_gpu, f_temp_1_gpu, f_temp_2_gpu, f_temp_3_gpu, f_temp_4_gpu,
                    f_temp_5_gpu, f_temp_6_gpu, f_temp_7_gpu, f_temp_8_gpu,
                    dt, k, s_0,
                    start, size);

            CUDA_ERROR();
        }
}

extern "C" void cuda_force_grid_pollutant_float(
        unsigned long start, unsigned long end,
        void * h, void * c,
        void * f_temp_0, void * f_temp_1, void * f_temp_2,
        void * f_temp_3, void * f_temp_4, void * f_temp_5,
        void * f_temp_6, void * f_temp_7, void * f_temp_8,
        float dt, float k, float s_0,
        unsigned long blocksize)
{
    honei::cuda_force_grid_pollutant<float>(start, end, h, c, f_temp_0, f_temp_1, f_temp_2, f_temp_3, f_temp_4, f_temp_5, f_temp_6, f_temp_7, f_temp_8,
            dt, k, s_0, blocksize);
}

#ifdef HONEI_CUDA_DOUBLE
extern "C" void cuda_force_grid_pollutant_double(
        unsigned long start, unsigned long end,
        void * h, void * c,
        void * f_temp_0, void * f_temp_1, void * f_temp_2,
        void * f_temp_3, void * f_temp_4, void * f_temp_5,
        void * f_temp_6, void * f_temp_7, void * f_temp_8,
        double dt, double k, double s_0,
        unsigned long blocksize)
{
    honei::cuda_force_grid_pollutant<double>(start, end, h, c, f_temp_0, f_temp_1, f_temp_2, f_temp_3, f_temp_4, f_temp_5, f_temp_6, f_temp_7, f_temp_8,
            dt, k, s_0, blocksize);
}
#endif
