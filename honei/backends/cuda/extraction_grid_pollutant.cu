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
        __global__ void extraction_grid_pollutant_dry_gpu(
                DT_ * f_0, DT_ * f_1, DT_ * f_2,
                DT_ * f_3, DT_ * f_4, DT_ * f_5,
                DT_ * f_6, DT_ * f_7, DT_ * f_8,
                DT_ * h, DT_ * c, DT_ epsilon,
                unsigned long offset, unsigned long size)
        {
            unsigned long idx = (blockDim.y * blockIdx.y * gridDim.x * blockDim.x) + (blockDim.x * blockIdx.x) + threadIdx.x;

            if (idx >= offset && idx < size)
            {
                unsigned long i(idx);

                c[i] = (f_0[i] +
                    f_1[i] +
                    f_2[i] +
                    f_3[i] +
                    f_4[i] +
                    f_5[i] +
                    f_6[i] +
                    f_7[i] +
                    f_8[i]) / h[i];
            }
        }
    }

    template <typename DT_>
        void cuda_extraction_grid_pollutant_dry(
                unsigned long start, unsigned long end,
                void * f_0, void * f_1, void * f_2,
                void * f_3, void * f_4, void * f_5,
                void * f_6, void * f_7, void * f_8,
                void * h, void * c, DT_ epsilon,
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
            DT_ * c_gpu((DT_ *)c);

            honei::cuda::extraction_grid_pollutant_dry_gpu<DT_><<<grid, block>>>(
                    f_0_gpu, f_1_gpu, f_2_gpu, f_3_gpu, f_4_gpu,
                    f_5_gpu, f_6_gpu, f_7_gpu, f_8_gpu,
                    h_gpu, c_gpu, epsilon,
                    start, size);

            CUDA_ERROR();
        }
}

extern "C" void cuda_extraction_grid_pollutant_dry_float(
        unsigned long start, unsigned long end,
        void * f_0, void * f_1, void * f_2,
        void * f_3, void * f_4, void * f_5,
        void * f_6, void * f_7, void * f_8,
        void * h, void * c, float epsilon,
        unsigned long blocksize)
{
    honei::cuda_extraction_grid_pollutant_dry<float>(start, end, f_0, f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, h, c,
            epsilon, blocksize);
}

#ifdef HONEI_CUDA_DOUBLE
extern "C" void cuda_extraction_grid_pollutant_dry_double(
        unsigned long start, unsigned long end,
        void * f_0, void * f_1, void * f_2,
        void * f_3, void * f_4, void * f_5,
        void * f_6, void * f_7, void * f_8,
        void * h, void * c, double epsilon,
        unsigned long blocksize)
{
    honei::cuda_extraction_grid_pollutant_dry<double>(start, end, f_0, f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, h, c,
            epsilon, blocksize);
}
#endif
