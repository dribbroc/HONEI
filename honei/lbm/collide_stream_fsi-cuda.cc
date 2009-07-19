/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c)  2009 Markus Geveler <apryde@gmx.de>
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

#include <honei/lbm/collide_stream_fsi.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/util/configuration.hh>

using namespace honei;
using namespace lbm;

void CollideStreamFSI<tags::GPU::CUDA, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>::value(
                            PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                            PackedGridData<lbm_lattice_types::D2Q9, float> & data,
                            PackedSolidData<lbm_lattice_types::D2Q9, float> & solids,
                            float d_x,
                            float d_y)
{
    CONTEXT("When performing collision and streaming FSI correction (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::collide_stream_grid_float", 128ul));

    info.limits->lock(lm_read_only);

    void * cuda_dir_1_gpu(info.cuda_dir_1->lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * cuda_dir_2_gpu(info.cuda_dir_2->lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * cuda_dir_3_gpu(info.cuda_dir_3->lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * cuda_dir_4_gpu(info.cuda_dir_4->lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * cuda_dir_5_gpu(info.cuda_dir_5->lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * cuda_dir_6_gpu(info.cuda_dir_6->lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * cuda_dir_7_gpu(info.cuda_dir_7->lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * cuda_dir_8_gpu(info.cuda_dir_8->lock(lm_read_only, tags::GPU::CUDA::memory_value));

    void * f_temp_1_gpu(data.f_temp_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
    void * f_temp_2_gpu(data.f_temp_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
    void * f_temp_3_gpu(data.f_temp_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
    void * f_temp_4_gpu(data.f_temp_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
    void * f_temp_5_gpu(data.f_temp_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
    void * f_temp_6_gpu(data.f_temp_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
    void * f_temp_7_gpu(data.f_temp_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
    void * f_temp_8_gpu(data.f_temp_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

    void * f_mea_1_gpu(solids.f_mea_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
    void * f_mea_2_gpu(solids.f_mea_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
    void * f_mea_3_gpu(solids.f_mea_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
    void * f_mea_4_gpu(solids.f_mea_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
    void * f_mea_5_gpu(solids.f_mea_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
    void * f_mea_6_gpu(solids.f_mea_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
    void * f_mea_7_gpu(solids.f_mea_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
    void * f_mea_8_gpu(solids.f_mea_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

    void * line_flags_gpu(solids.solid_flags->lock(lm_read_only, tags::GPU::CUDA::memory_value));//We need indeed SOLID flags
    void * dist_x_gpu(data.distribution_x->lock(lm_read_only, tags::GPU::CUDA::memory_value));
    void * dist_y_gpu(data.distribution_y->lock(lm_read_only, tags::GPU::CUDA::memory_value));

    unsigned long start((*info.limits)[0]);
    unsigned long end((*info.limits)[info.limits->size() - 1]);

    float u(solids.current_u);
    float v(solids.current_v);

    cuda_collide_stream_fsi_float(start, end,
            cuda_dir_1_gpu, cuda_dir_2_gpu, cuda_dir_3_gpu, cuda_dir_4_gpu,
            cuda_dir_5_gpu, cuda_dir_6_gpu, cuda_dir_7_gpu, cuda_dir_8_gpu,
            f_temp_1_gpu, f_temp_2_gpu,
            f_temp_3_gpu, f_temp_4_gpu, f_temp_5_gpu,
            f_temp_6_gpu, f_temp_7_gpu, f_temp_8_gpu,
            f_mea_1_gpu, f_mea_2_gpu,
            f_mea_3_gpu, f_mea_4_gpu, f_mea_5_gpu,
            f_mea_6_gpu, f_mea_7_gpu, f_mea_8_gpu,
            line_flags_gpu, dist_x_gpu, dist_y_gpu,
            d_x * u, d_y * v,
            data.h->size(),
            blocksize);


    info.limits->unlock(lm_read_only);

    info.cuda_dir_1->unlock(lm_read_only);
    info.cuda_dir_2->unlock(lm_read_only);
    info.cuda_dir_3->unlock(lm_read_only);
    info.cuda_dir_4->unlock(lm_read_only);
    info.cuda_dir_5->unlock(lm_read_only);
    info.cuda_dir_6->unlock(lm_read_only);
    info.cuda_dir_7->unlock(lm_read_only);
    info.cuda_dir_8->unlock(lm_read_only);

    data.f_temp_1->unlock(lm_read_and_write);
    data.f_temp_2->unlock(lm_read_and_write);
    data.f_temp_3->unlock(lm_read_and_write);
    data.f_temp_4->unlock(lm_read_and_write);
    data.f_temp_5->unlock(lm_read_and_write);
    data.f_temp_6->unlock(lm_read_and_write);
    data.f_temp_7->unlock(lm_read_and_write);
    data.f_temp_8->unlock(lm_read_and_write);

    solids.f_mea_1->unlock(lm_read_and_write);
    solids.f_mea_2->unlock(lm_read_and_write);
    solids.f_mea_3->unlock(lm_read_and_write);
    solids.f_mea_4->unlock(lm_read_and_write);
    solids.f_mea_5->unlock(lm_read_and_write);
    solids.f_mea_6->unlock(lm_read_and_write);
    solids.f_mea_7->unlock(lm_read_and_write);
    solids.f_mea_8->unlock(lm_read_and_write);

    solids.solid_flags->unlock(lm_read_only);
    data.distribution_x->unlock(lm_read_only);
    data.distribution_y->unlock(lm_read_only);
}

