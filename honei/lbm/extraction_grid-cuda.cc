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

#include <honei/lbm/extraction_grid.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/util/configuration.hh>


using namespace honei;

void ExtractionGrid<tags::GPU::CUDA, lbm_applications::LABSWE>::value(
                PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                PackedGridData<lbm_lattice_types::D2Q9, float> & data)
{
    CONTEXT("When extracting h, u and v (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::extraction_grid_float", 128ul));

    info.limits->lock(lm_read_only);

    void * f_0_gpu(data.f_0->lock(lm_write_only, tags::GPU::CUDA::memory_value));
    void * f_1_gpu(data.f_1->lock(lm_write_only, tags::GPU::CUDA::memory_value));
    void * f_2_gpu(data.f_2->lock(lm_write_only, tags::GPU::CUDA::memory_value));
    void * f_3_gpu(data.f_3->lock(lm_write_only, tags::GPU::CUDA::memory_value));
    void * f_4_gpu(data.f_4->lock(lm_write_only, tags::GPU::CUDA::memory_value));
    void * f_5_gpu(data.f_5->lock(lm_write_only, tags::GPU::CUDA::memory_value));
    void * f_6_gpu(data.f_6->lock(lm_write_only, tags::GPU::CUDA::memory_value));
    void * f_7_gpu(data.f_7->lock(lm_write_only, tags::GPU::CUDA::memory_value));
    void * f_8_gpu(data.f_8->lock(lm_write_only, tags::GPU::CUDA::memory_value));

    void * f_temp_0_gpu(data.f_temp_0->lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * f_temp_1_gpu(data.f_temp_1->lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * f_temp_2_gpu(data.f_temp_2->lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * f_temp_3_gpu(data.f_temp_3->lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * f_temp_4_gpu(data.f_temp_4->lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * f_temp_5_gpu(data.f_temp_5->lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * f_temp_6_gpu(data.f_temp_6->lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * f_temp_7_gpu(data.f_temp_7->lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * f_temp_8_gpu(data.f_temp_8->lock(lm_read_only, tags::GPU::CUDA::memory_value));

    void * h_gpu(data.h->lock(lm_write_only, tags::GPU::CUDA::memory_value));
    void * u_gpu(data.u->lock(lm_write_only, tags::GPU::CUDA::memory_value));
    void * v_gpu(data.v->lock(lm_write_only, tags::GPU::CUDA::memory_value));

    void * distribution_x_gpu(data.distribution_x->lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * distribution_y_gpu(data.distribution_y->lock(lm_read_only, tags::GPU::CUDA::memory_value));

    unsigned long start((*info.limits)[0]);
    unsigned long end((*info.limits)[info.limits->size() - 1]);

    cuda_extraction_grid_float(start, end,
            f_0_gpu, f_1_gpu, f_2_gpu,
            f_3_gpu, f_4_gpu, f_5_gpu,
            f_6_gpu, f_7_gpu, f_8_gpu,
            f_temp_0_gpu, f_temp_1_gpu, f_temp_2_gpu,
            f_temp_3_gpu, f_temp_4_gpu, f_temp_5_gpu,
            f_temp_6_gpu, f_temp_7_gpu, f_temp_8_gpu,
            h_gpu, u_gpu, v_gpu,
            distribution_x_gpu, distribution_y_gpu,
            blocksize);

    info.limits->unlock(lm_read_only);

    data.f_temp_0->unlock(lm_read_only);
    data.f_temp_1->unlock(lm_read_only);
    data.f_temp_2->unlock(lm_read_only);
    data.f_temp_3->unlock(lm_read_only);
    data.f_temp_4->unlock(lm_read_only);
    data.f_temp_5->unlock(lm_read_only);
    data.f_temp_6->unlock(lm_read_only);
    data.f_temp_7->unlock(lm_read_only);
    data.f_temp_8->unlock(lm_read_only);

    data.f_0->unlock(lm_write_only);
    data.f_1->unlock(lm_write_only);
    data.f_2->unlock(lm_write_only);
    data.f_3->unlock(lm_write_only);
    data.f_4->unlock(lm_write_only);
    data.f_5->unlock(lm_write_only);
    data.f_6->unlock(lm_write_only);
    data.f_7->unlock(lm_write_only);
    data.f_8->unlock(lm_write_only);

    data.h->unlock(lm_write_only);

    data.distribution_x->unlock(lm_read_only);
    data.distribution_y->unlock(lm_read_only);

    data.u->unlock(lm_write_only);
    data.v->unlock(lm_write_only);
}

