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
        __global__ void eq_dist_grid_pollutant_gpu(DT_ * tu, DT_ * tv, DT_ * th, DT_ * tc,
                DT_ * distribution_x_data, DT_ * distribution_y_data,
                DT_ * f_eq_0, DT_ * f_eq_1, DT_ * f_eq_2,
                DT_ * f_eq_3, DT_ * f_eq_4, DT_ * f_eq_5,
                DT_ * f_eq_6, DT_ * f_eq_7, DT_ * f_eq_8,
                DT_ e,
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
                DT_ e2_by_3(DT_(3.) * e2);
                DT_ e2_by_12(DT_(12.) * e2);
                DT_ one_over_9(1./9.);
                DT_ five_over_9(5./9.);
                DT_ one_over_36(1./36.);

                DT_ dxu, dyv;
                DT_ t1, t2;

                const DT_ u((tu)[index]);
                const DT_ v((tv)[index]);
                const DT_ h((th)[index]);
                const DT_ c((tc)[index]);

                (f_eq_0)[index] = c * DT_(h - five_over_9);


                dxu = (distribution_x)[1] * u;
                dyv = (distribution_y)[1] * v;
                t1 = (h / e2_by_3) * dxu;
                t2 = (h / e2_by_3) * dyv;
                (f_eq_1)[index] = c * (one_over_9 + t1 + t2);

                dxu = (distribution_x)[3] * u;
                dyv = (distribution_y)[3] * v;
                t1 = (h / e2_by_3) * dxu;
                t2 = (h / e2_by_3) * dyv;
                (f_eq_3)[index] = c * (one_over_9 + t1 + t2);

                dxu = (distribution_x)[5] * u;
                dyv = (distribution_y)[5] * v;
                t1 = (h / e2_by_3) * dxu;
                t2 = (h / e2_by_3) * dyv;
                (f_eq_5)[index] = c * (one_over_9 + t1 + t2);

                dxu = (distribution_x)[7] * u;
                dyv = (distribution_y)[7] * v;
                t1 = (h / e2_by_3) * dxu;
                t2 = (h / e2_by_3) * dyv;
                (f_eq_7)[index] = c * (one_over_9 + t1 + t2);

                dxu = (distribution_x)[2] * u;
                dyv = (distribution_y)[2] * v;
                t1 = (h / e2_by_12) * dxu;
                t2 = (h / e2_by_12) * dyv;
                (f_eq_2)[index] = c * (one_over_36 + t1 + t2);

                dxu = (distribution_x)[4] * u;
                dyv = (distribution_y)[4] * v;
                t1 = (h / e2_by_12) * dxu;
                t2 = (h / e2_by_12) * dyv;
                (f_eq_4)[index] = c * (one_over_36 + t1 + t2);

                dxu = (distribution_x)[6] * u;
                dyv = (distribution_y)[6] * v;
                t1 = (h / e2_by_12) * dxu;
                t2 = (h / e2_by_12) * dyv;
                (f_eq_6)[index] = c * (one_over_36 + t1 + t2);

                dxu = (distribution_x)[8] * u;
                dyv = (distribution_y)[8] * v;
                t1 = (h / e2_by_12) * dxu;
                t2 = (h / e2_by_12) * dyv;
                (f_eq_8)[index] = c * (one_over_36 + t1 + t2);
            }
        }
    }

    template <typename DT_>
    void cuda_eq_dist_grid_pollutant(unsigned long start, unsigned long end, void * u, void * v, void * h, void * c,
            void * distribution_x, void * distribution_y,
            void * f_eq_0, void * f_eq_1, void * f_eq_2,
            void * f_eq_3, void * f_eq_4, void * f_eq_5,
            void * f_eq_6, void * f_eq_7, void * f_eq_8,
            DT_ e,
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
        DT_ * c_gpu((DT_ *)c);
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

        honei::cuda::eq_dist_grid_pollutant_gpu<DT_><<<grid, block, 18 * sizeof(DT_)>>>(u_gpu, v_gpu, h_gpu, c_gpu,
                distribution_x_gpu, distribution_y_gpu,
                f_eq_0_gpu, f_eq_1_gpu, f_eq_2_gpu,
                f_eq_3_gpu, f_eq_4_gpu, f_eq_5_gpu,
                f_eq_6_gpu, f_eq_7_gpu, f_eq_8_gpu,
                e,
                start, size);

        CUDA_ERROR();
    }
}

extern "C" void cuda_eq_dist_grid_pollutant_float(unsigned long start, unsigned long end, void * u, void * v, void * h, void * c,
        void * distribution_x, void * distribution_y,
        void * f_eq_0, void * f_eq_1, void * f_eq_2,
        void * f_eq_3, void * f_eq_4, void * f_eq_5,
        void * f_eq_6, void * f_eq_7, void * f_eq_8,
        float e,
        unsigned long blocksize)
{
    honei::cuda_eq_dist_grid_pollutant<float>(start, end, u, v, h, c, distribution_x, distribution_y, f_eq_0, f_eq_1, f_eq_2, f_eq_3, f_eq_4, f_eq_5,
            f_eq_6, f_eq_7, f_eq_8, e, blocksize);
}

#ifdef HONEI_CUDA_DOUBLE
extern "C" void cuda_eq_dist_grid_pollutant_double(unsigned long start, unsigned long end, void * u, void * v, void * h, void * c,
        void * distribution_x, void * distribution_y,
        void * f_eq_0, void * f_eq_1, void * f_eq_2,
        void * f_eq_3, void * f_eq_4, void * f_eq_5,
        void * f_eq_6, void * f_eq_7, void * f_eq_8,
        double e,
        unsigned long blocksize)
{
    honei::cuda_eq_dist_grid_pollutant<double>(start, end, u, v, h, c, distribution_x, distribution_y, f_eq_0, f_eq_1, f_eq_2, f_eq_3, f_eq_4, f_eq_5,
            f_eq_6, f_eq_7, f_eq_8, e, blocksize);
}
#endif
