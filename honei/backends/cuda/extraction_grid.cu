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
        __global__ void extraction_grid_dry_gpu(
                DT_ * f_0, DT_ * f_1, DT_ * f_2,
                DT_ * f_3, DT_ * f_4, DT_ * f_5,
                DT_ * f_6, DT_ * f_7, DT_ * f_8,
                DT_ * h, DT_ * u, DT_ * v,
                DT_ * distribution_x_data, DT_ * distribution_y_data, DT_ epsilon,
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
                unsigned long i(idx);

                h[i] = f_0[i] +
                    f_1[i] +
                    f_2[i] +
                    f_3[i] +
                    f_4[i] +
                    f_5[i] +
                    f_6[i] +
                    f_7[i] +
                    f_8[i];

                DT_ lax_upper(epsilon);
                DT_ lax_lower(-lax_upper);

                if(h[i] < lax_lower || h[i] > lax_upper)
                {
                    u[i] = (distribution_x[0] * f_0[i] +
                            distribution_x[1] * f_1[i] +
                            distribution_x[2] * f_2[i] +
                            distribution_x[3] * f_3[i] +
                            distribution_x[4] * f_4[i] +
                            distribution_x[5] * f_5[i] +
                            distribution_x[6] * f_6[i] +
                            distribution_x[7] * f_7[i] +
                            distribution_x[8] * f_8[i]) / h[i];

                    v[i] = (distribution_y[0] * f_0[i] +
                            distribution_y[1] * f_1[i] +
                            distribution_y[2] * f_2[i] +
                            distribution_y[3] * f_3[i] +
                            distribution_y[4] * f_4[i] +
                            distribution_y[5] * f_5[i] +
                            distribution_y[6] * f_6[i] +
                            distribution_y[7] * f_7[i] +
                            distribution_y[8] * f_8[i]) / h[i];
                }
                else
                {
                    h[i] = 0;
                    u[i] = 0;
                    v[i] = 0;
                }
                h[i] = max(DT_(0), min(DT_(1), h[i]));
            }
        }

        template <typename DT_>
        __global__ void extraction_grid_wet_gpu(
                DT_ * f_0, DT_ * f_1, DT_ * f_2,
                DT_ * f_3, DT_ * f_4, DT_ * f_5,
                DT_ * f_6, DT_ * f_7, DT_ * f_8,
                DT_ * h, DT_ * u, DT_ * v,
                DT_ * distribution_x_data, DT_ * distribution_y_data, DT_ epsilon,
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
                unsigned long i(idx);

                h[i] = f_0[i] +
                    f_1[i] +
                    f_2[i] +
                    f_3[i] +
                    f_4[i] +
                    f_5[i] +
                    f_6[i] +
                    f_7[i] +
                    f_8[i];

                u[i] = (distribution_x[0] * f_0[i] +
                        distribution_x[1] * f_1[i] +
                        distribution_x[2] * f_2[i] +
                        distribution_x[3] * f_3[i] +
                        distribution_x[4] * f_4[i] +
                        distribution_x[5] * f_5[i] +
                        distribution_x[6] * f_6[i] +
                        distribution_x[7] * f_7[i] +
                        distribution_x[8] * f_8[i]) / h[i];

                v[i] = (distribution_y[0] * f_0[i] +
                        distribution_y[1] * f_1[i] +
                        distribution_y[2] * f_2[i] +
                        distribution_y[3] * f_3[i] +
                        distribution_y[4] * f_4[i] +
                        distribution_y[5] * f_5[i] +
                        distribution_y[6] * f_6[i] +
                        distribution_y[7] * f_7[i] +
                        distribution_y[8] * f_8[i]) / h[i];
            }
        }
    }

    template <typename DT_>
        void cuda_extraction_grid_dry(
                unsigned long start, unsigned long end,
                void * f_0, void * f_1, void * f_2,
                void * f_3, void * f_4, void * f_5,
                void * f_6, void * f_7, void * f_8,
                void * h, void * u, void * v,
                void * distribution_x, void * distribution_y, DT_ epsilon,
                unsigned long blocksize)
        {
            unsigned long size(end);
            dim3 grid;
            dim3 block;
            block.x = blocksize;
            grid.x = (unsigned)ceil(sqrt(size/(double)block.x));
            grid.y = grid.x;

            DT_ * f_0_gpu((DT_ *)f_0);
            DT_ * f_1_gpu((DT_ *)f_1);
            DT_ * f_2_gpu((DT_ *)f_2);
            DT_ * f_3_gpu((DT_ *)f_3);
            DT_ * f_4_gpu((DT_ *)f_4);
            DT_ * f_5_gpu((DT_ *)f_5);
            DT_ * f_6_gpu((DT_ *)f_6);
            DT_ * f_7_gpu((DT_ *)f_7);
            DT_ * f_8_gpu((DT_ *)f_8);

            DT_ * h_gpu((DT_ *)h);
            DT_ * u_gpu((DT_ *)u);
            DT_ * v_gpu((DT_ *)v);
            DT_ * distribution_x_gpu((DT_ *)distribution_x);
            DT_ * distribution_y_gpu((DT_ *)distribution_y);

            honei::cuda::extraction_grid_dry_gpu<DT_><<<grid, block, 18 * sizeof(DT_)>>>(
                    f_0_gpu, f_1_gpu, f_2_gpu, f_3_gpu, f_4_gpu,
                    f_5_gpu, f_6_gpu, f_7_gpu, f_8_gpu,
                    h_gpu, u_gpu, v_gpu,
                    distribution_x_gpu, distribution_y_gpu, epsilon,
                    start, size);

            CUDA_ERROR();
        }

    template <typename DT_>
        void cuda_extraction_grid_wet(
                unsigned long start, unsigned long end,
                void * f_0, void * f_1, void * f_2,
                void * f_3, void * f_4, void * f_5,
                void * f_6, void * f_7, void * f_8,
                void * h, void * u, void * v,
                void * distribution_x, void * distribution_y, DT_ epsilon,
                unsigned long blocksize)
        {
            unsigned long size(end);
            dim3 grid;
            dim3 block;
            block.x = blocksize;
            grid.x = (unsigned)ceil(sqrt(size/(double)block.x));
            grid.y = grid.x;

            DT_ * f_0_gpu((DT_ *)f_0);
            DT_ * f_1_gpu((DT_ *)f_1);
            DT_ * f_2_gpu((DT_ *)f_2);
            DT_ * f_3_gpu((DT_ *)f_3);
            DT_ * f_4_gpu((DT_ *)f_4);
            DT_ * f_5_gpu((DT_ *)f_5);
            DT_ * f_6_gpu((DT_ *)f_6);
            DT_ * f_7_gpu((DT_ *)f_7);
            DT_ * f_8_gpu((DT_ *)f_8);

            DT_ * h_gpu((DT_ *)h);
            DT_ * u_gpu((DT_ *)u);
            DT_ * v_gpu((DT_ *)v);
            DT_ * distribution_x_gpu((DT_ *)distribution_x);
            DT_ * distribution_y_gpu((DT_ *)distribution_y);

            honei::cuda::extraction_grid_wet_gpu<DT_><<<grid, block, 18 * sizeof(DT_)>>>(
                    f_0_gpu, f_1_gpu, f_2_gpu, f_3_gpu, f_4_gpu,
                    f_5_gpu, f_6_gpu, f_7_gpu, f_8_gpu,
                    h_gpu, u_gpu, v_gpu,
                    distribution_x_gpu, distribution_y_gpu, epsilon,
                    start, size);

            CUDA_ERROR();
        }
}

extern "C" void cuda_extraction_grid_dry_float(
        unsigned long start, unsigned long end,
        void * f_0, void * f_1, void * f_2,
        void * f_3, void * f_4, void * f_5,
        void * f_6, void * f_7, void * f_8,
        void * h, void * u, void * v,
        void * distribution_x, void * distribution_y, float epsilon,
        unsigned long blocksize)
{
    honei::cuda_extraction_grid_dry<float>(start, end, f_0, f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, h, u, v,
            distribution_x, distribution_y, epsilon, blocksize);
}

#ifdef HONEI_CUDA_DOUBLE
extern "C" void cuda_extraction_grid_dry_double(
        unsigned long start, unsigned long end,
        void * f_0, void * f_1, void * f_2,
        void * f_3, void * f_4, void * f_5,
        void * f_6, void * f_7, void * f_8,
        void * h, void * u, void * v,
        void * distribution_x, void * distribution_y, double epsilon,
        unsigned long blocksize)
{
    honei::cuda_extraction_grid_dry<double>(start, end, f_0, f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, h, u, v,
            distribution_x, distribution_y, epsilon, blocksize);
}
#endif

extern "C" void cuda_extraction_grid_wet_float(
        unsigned long start, unsigned long end,
        void * f_0, void * f_1, void * f_2,
        void * f_3, void * f_4, void * f_5,
        void * f_6, void * f_7, void * f_8,
        void * h, void * u, void * v,
        void * distribution_x, void * distribution_y, float epsilon,
        unsigned long blocksize)
{
    honei::cuda_extraction_grid_wet<float>(start, end, f_0, f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, h, u, v,
            distribution_x, distribution_y, epsilon, blocksize);
}

#ifdef HONEI_CUDA_DOUBLE
extern "C" void cuda_extraction_grid_wet_double(
        unsigned long start, unsigned long end,
        void * f_0, void * f_1, void * f_2,
        void * f_3, void * f_4, void * f_5,
        void * f_6, void * f_7, void * f_8,
        void * h, void * u, void * v,
        void * distribution_x, void * distribution_y, double epsilon,
        unsigned long blocksize)
{
    honei::cuda_extraction_grid_wet<double>(start, end, f_0, f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, h, u, v,
            distribution_x, distribution_y, epsilon, blocksize);
}
#endif
