/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c)  2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/lbm/force_grid.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/util/configuration.hh>


using namespace honei;

void ForceGrid<tags::GPU::CUDA, lbm_applications::LABSWE, lbm_force::CENTRED, lbm_source_schemes::BED_SLOPE>::value(
        PackedGridData<D2Q9, float> & data, PackedGridInfo<D2Q9> & info, float g, float d_x, float d_y, float d_t, float manning)
{
    CONTEXT("When computing LABSWE source term (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::force_grid_float", 128ul));

    void * cuda_dir_1_gpu(info.cuda_dir_1->lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * cuda_dir_2_gpu(info.cuda_dir_2->lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * cuda_dir_3_gpu(info.cuda_dir_3->lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * cuda_dir_4_gpu(info.cuda_dir_4->lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * cuda_dir_5_gpu(info.cuda_dir_5->lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * cuda_dir_6_gpu(info.cuda_dir_6->lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * cuda_dir_7_gpu(info.cuda_dir_7->lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * cuda_dir_8_gpu(info.cuda_dir_8->lock(lm_read_only, tags::GPU::CUDA::memory_value));

    void * h_gpu(data.h->lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * b_gpu(data.b->lock(lm_read_only, tags::GPU::CUDA::memory_value));

    void * f_temp_1_gpu(data.f_temp_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
    void * f_temp_2_gpu(data.f_temp_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
    void * f_temp_3_gpu(data.f_temp_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
    void * f_temp_4_gpu(data.f_temp_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
    void * f_temp_5_gpu(data.f_temp_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
    void * f_temp_6_gpu(data.f_temp_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
    void * f_temp_7_gpu(data.f_temp_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
    void * f_temp_8_gpu(data.f_temp_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

    cuda_force_grid_float(
            cuda_dir_1_gpu, cuda_dir_2_gpu, cuda_dir_3_gpu, cuda_dir_4_gpu,
            cuda_dir_5_gpu, cuda_dir_6_gpu, cuda_dir_7_gpu, cuda_dir_8_gpu,
            h_gpu, b_gpu,
            f_temp_1_gpu, f_temp_2_gpu,
            f_temp_3_gpu, f_temp_4_gpu, f_temp_5_gpu,
            f_temp_6_gpu, f_temp_7_gpu, f_temp_8_gpu,
            (*data.distribution_x)[1], (*data.distribution_y)[1],
            (*data.distribution_x)[2], (*data.distribution_y)[2],
            (*data.distribution_x)[3], (*data.distribution_y)[3],
            (*data.distribution_x)[4], (*data.distribution_y)[4],
            (*data.distribution_x)[5], (*data.distribution_y)[5],
            (*data.distribution_x)[6], (*data.distribution_y)[6],
            (*data.distribution_x)[7], (*data.distribution_y)[7],
            (*data.distribution_x)[8], (*data.distribution_y)[8],
            g, d_x, d_y, d_t,
            data.h->size(),
            blocksize);

    info.cuda_dir_1->unlock(lm_read_only);
    info.cuda_dir_2->unlock(lm_read_only);
    info.cuda_dir_3->unlock(lm_read_only);
    info.cuda_dir_4->unlock(lm_read_only);
    info.cuda_dir_5->unlock(lm_read_only);
    info.cuda_dir_6->unlock(lm_read_only);
    info.cuda_dir_7->unlock(lm_read_only);
    info.cuda_dir_8->unlock(lm_read_only);

    data.h->unlock(lm_read_only);
    data.b->unlock(lm_read_only);

    data.f_temp_1->unlock(lm_read_and_write);
    data.f_temp_2->unlock(lm_read_and_write);
    data.f_temp_3->unlock(lm_read_and_write);
    data.f_temp_4->unlock(lm_read_and_write);
    data.f_temp_5->unlock(lm_read_and_write);
    data.f_temp_6->unlock(lm_read_and_write);
    data.f_temp_7->unlock(lm_read_and_write);
    data.f_temp_8->unlock(lm_read_and_write);
}
