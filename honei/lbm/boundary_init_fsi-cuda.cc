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

#include <honei/lbm/boundary_init_fsi.hh>
#include <honei/backends/cuda/operations.hh>
#include <honei/backends/cuda/gpu_pool.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/util/configuration.hh>

using namespace honei;
using namespace lbm;

namespace
{
    class cudaBoundaryInitFSIGridfloatDir1
    {
        private:
            PackedGridInfo<D2Q9> & info;
            PackedGridData<D2Q9, float> & data;
            PackedSolidData<lbm_lattice_types::D2Q9, float> & solids;
            unsigned long blocksize;
        public:
            cudaBoundaryInitFSIGridfloatDir1(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data, PackedSolidData<lbm_lattice_types::D2Q9, float> & solids, unsigned long blocksize) :
                info(info),
                data(data),
                solids(solids),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * cuda_dir_5_gpu(info.cuda_dir_5->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * h_gpu(data.h->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * u_gpu(data.u->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * v_gpu(data.v->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

                void * boundary_flags_gpu(solids.boundary_flags->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * solid_flags_gpu(solids.solid_flags->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * solid_old_flags_gpu(solids.solid_old_flags->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

                void * f_0_gpu(data.f_0->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_1_gpu(data.f_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_2_gpu(data.f_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_3_gpu(data.f_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_4_gpu(data.f_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_5_gpu(data.f_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_6_gpu(data.f_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_7_gpu(data.f_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_8_gpu(data.f_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_0_gpu(data.f_eq_0->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_1_gpu(data.f_eq_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_2_gpu(data.f_eq_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_3_gpu(data.f_eq_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_4_gpu(data.f_eq_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_5_gpu(data.f_eq_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_6_gpu(data.f_eq_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_7_gpu(data.f_eq_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_8_gpu(data.f_eq_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_0_gpu(data.f_temp_0->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_1_gpu(data.f_temp_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_2_gpu(data.f_temp_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_3_gpu(data.f_temp_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_4_gpu(data.f_temp_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_5_gpu(data.f_temp_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_6_gpu(data.f_temp_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_7_gpu(data.f_temp_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_8_gpu(data.f_temp_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

                cuda_boundary_init_fsi_dir_1_float(
                        cuda_dir_5_gpu,
                        boundary_flags_gpu,
                        solid_flags_gpu,
                        solid_old_flags_gpu,
                        h_gpu,
                        u_gpu,
                        v_gpu,
                        f_0_gpu,
                        f_1_gpu,
                        f_2_gpu,
                        f_3_gpu,
                        f_4_gpu,
                        f_5_gpu,
                        f_6_gpu,
                        f_7_gpu,
                        f_8_gpu,
                        f_eq_0_gpu,
                        f_eq_1_gpu,
                        f_eq_2_gpu,
                        f_eq_3_gpu,
                        f_eq_4_gpu,
                        f_eq_5_gpu,
                        f_eq_6_gpu,
                        f_eq_7_gpu,
                        f_eq_8_gpu,
                        f_temp_0_gpu,
                        f_temp_1_gpu,
                        f_temp_2_gpu,
                        f_temp_3_gpu,
                        f_temp_4_gpu,
                        f_temp_5_gpu,
                        f_temp_6_gpu,
                        f_temp_7_gpu,
                        f_temp_8_gpu,
                        solids.current_u,
                        solids.current_v,
                        data.h->size(),
                        blocksize);

                info.cuda_dir_5->unlock(lm_read_only);

                solids.boundary_flags->unlock(lm_read_only);
                solids.solid_flags->unlock(lm_read_only);

                solids.solid_old_flags->unlock(lm_read_and_write);

                data.h->unlock(lm_read_and_write);
                data.u->unlock(lm_read_and_write);
                data.v->unlock(lm_read_and_write);

                data.f_0->unlock(lm_read_and_write);
                data.f_1->unlock(lm_read_and_write);
                data.f_2->unlock(lm_read_and_write);
                data.f_3->unlock(lm_read_and_write);
                data.f_4->unlock(lm_read_and_write);
                data.f_5->unlock(lm_read_and_write);
                data.f_6->unlock(lm_read_and_write);
                data.f_7->unlock(lm_read_and_write);
                data.f_8->unlock(lm_read_and_write);
                data.f_eq_0->unlock(lm_read_and_write);
                data.f_eq_1->unlock(lm_read_and_write);
                data.f_eq_2->unlock(lm_read_and_write);
                data.f_eq_3->unlock(lm_read_and_write);
                data.f_eq_4->unlock(lm_read_and_write);
                data.f_eq_5->unlock(lm_read_and_write);
                data.f_eq_6->unlock(lm_read_and_write);
                data.f_eq_7->unlock(lm_read_and_write);
                data.f_eq_8->unlock(lm_read_and_write);
                data.f_temp_0->unlock(lm_read_and_write);
                data.f_temp_1->unlock(lm_read_and_write);
                data.f_temp_2->unlock(lm_read_and_write);
                data.f_temp_3->unlock(lm_read_and_write);
                data.f_temp_4->unlock(lm_read_and_write);
                data.f_temp_5->unlock(lm_read_and_write);
                data.f_temp_6->unlock(lm_read_and_write);
                data.f_temp_7->unlock(lm_read_and_write);
                data.f_temp_8->unlock(lm_read_and_write);
            }
    };

    class cudaBoundaryInitFSIGridfloatDir2
    {
        private:
            PackedGridInfo<D2Q9> & info;
            PackedGridData<D2Q9, float> & data;
            PackedSolidData<lbm_lattice_types::D2Q9, float> & solids;
            unsigned long blocksize;
        public:
            cudaBoundaryInitFSIGridfloatDir2(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data, PackedSolidData<lbm_lattice_types::D2Q9, float> & solids, unsigned long blocksize) :
                info(info),
                data(data),
                solids(solids),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * cuda_dir_6_gpu(info.cuda_dir_6->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * h_gpu(data.h->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * u_gpu(data.u->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * v_gpu(data.v->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

                void * boundary_flags_gpu(solids.boundary_flags->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * solid_flags_gpu(solids.solid_flags->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * solid_old_flags_gpu(solids.solid_old_flags->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

                void * f_0_gpu(data.f_0->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_1_gpu(data.f_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_2_gpu(data.f_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_3_gpu(data.f_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_4_gpu(data.f_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_5_gpu(data.f_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_6_gpu(data.f_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_7_gpu(data.f_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_8_gpu(data.f_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_0_gpu(data.f_eq_0->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_1_gpu(data.f_eq_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_2_gpu(data.f_eq_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_3_gpu(data.f_eq_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_4_gpu(data.f_eq_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_5_gpu(data.f_eq_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_6_gpu(data.f_eq_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_7_gpu(data.f_eq_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_8_gpu(data.f_eq_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_0_gpu(data.f_temp_0->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_1_gpu(data.f_temp_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_2_gpu(data.f_temp_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_3_gpu(data.f_temp_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_4_gpu(data.f_temp_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_5_gpu(data.f_temp_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_6_gpu(data.f_temp_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_7_gpu(data.f_temp_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_8_gpu(data.f_temp_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

                cuda_boundary_init_fsi_dir_1_float(
                        cuda_dir_6_gpu,
                        boundary_flags_gpu,
                        solid_flags_gpu,
                        solid_old_flags_gpu,
                        h_gpu,
                        u_gpu,
                        v_gpu,
                        f_0_gpu,
                        f_1_gpu,
                        f_2_gpu,
                        f_3_gpu,
                        f_4_gpu,
                        f_5_gpu,
                        f_6_gpu,
                        f_7_gpu,
                        f_8_gpu,
                        f_eq_0_gpu,
                        f_eq_1_gpu,
                        f_eq_2_gpu,
                        f_eq_3_gpu,
                        f_eq_4_gpu,
                        f_eq_5_gpu,
                        f_eq_6_gpu,
                        f_eq_7_gpu,
                        f_eq_8_gpu,
                        f_temp_0_gpu,
                        f_temp_1_gpu,
                        f_temp_2_gpu,
                        f_temp_3_gpu,
                        f_temp_4_gpu,
                        f_temp_5_gpu,
                        f_temp_6_gpu,
                        f_temp_7_gpu,
                        f_temp_8_gpu,
                        solids.current_u,
                        solids.current_v,
                        data.h->size(),
                        blocksize);

                info.cuda_dir_6->unlock(lm_read_only);

                solids.boundary_flags->unlock(lm_read_only);
                solids.solid_flags->unlock(lm_read_only);

                solids.solid_old_flags->unlock(lm_read_and_write);

                data.h->unlock(lm_read_and_write);
                data.u->unlock(lm_read_and_write);
                data.v->unlock(lm_read_and_write);

                data.f_0->unlock(lm_read_and_write);
                data.f_1->unlock(lm_read_and_write);
                data.f_2->unlock(lm_read_and_write);
                data.f_3->unlock(lm_read_and_write);
                data.f_4->unlock(lm_read_and_write);
                data.f_5->unlock(lm_read_and_write);
                data.f_6->unlock(lm_read_and_write);
                data.f_7->unlock(lm_read_and_write);
                data.f_8->unlock(lm_read_and_write);
                data.f_eq_0->unlock(lm_read_and_write);
                data.f_eq_1->unlock(lm_read_and_write);
                data.f_eq_2->unlock(lm_read_and_write);
                data.f_eq_3->unlock(lm_read_and_write);
                data.f_eq_4->unlock(lm_read_and_write);
                data.f_eq_5->unlock(lm_read_and_write);
                data.f_eq_6->unlock(lm_read_and_write);
                data.f_eq_7->unlock(lm_read_and_write);
                data.f_eq_8->unlock(lm_read_and_write);
                data.f_temp_0->unlock(lm_read_and_write);
                data.f_temp_1->unlock(lm_read_and_write);
                data.f_temp_2->unlock(lm_read_and_write);
                data.f_temp_3->unlock(lm_read_and_write);
                data.f_temp_4->unlock(lm_read_and_write);
                data.f_temp_5->unlock(lm_read_and_write);
                data.f_temp_6->unlock(lm_read_and_write);
                data.f_temp_7->unlock(lm_read_and_write);
                data.f_temp_8->unlock(lm_read_and_write);
            }
    };

    class cudaBoundaryInitFSIGridfloatDir3
    {
        private:
            PackedGridInfo<D2Q9> & info;
            PackedGridData<D2Q9, float> & data;
            PackedSolidData<lbm_lattice_types::D2Q9, float> & solids;
            unsigned long blocksize;
        public:
            cudaBoundaryInitFSIGridfloatDir3(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data, PackedSolidData<lbm_lattice_types::D2Q9, float> & solids, unsigned long blocksize) :
                info(info),
                data(data),
                solids(solids),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * cuda_dir_7_gpu(info.cuda_dir_7->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * h_gpu(data.h->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * u_gpu(data.u->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * v_gpu(data.v->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

                void * boundary_flags_gpu(solids.boundary_flags->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * solid_flags_gpu(solids.solid_flags->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * solid_old_flags_gpu(solids.solid_old_flags->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

                void * f_0_gpu(data.f_0->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_1_gpu(data.f_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_2_gpu(data.f_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_3_gpu(data.f_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_4_gpu(data.f_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_5_gpu(data.f_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_6_gpu(data.f_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_7_gpu(data.f_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_8_gpu(data.f_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_0_gpu(data.f_eq_0->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_1_gpu(data.f_eq_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_2_gpu(data.f_eq_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_3_gpu(data.f_eq_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_4_gpu(data.f_eq_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_5_gpu(data.f_eq_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_6_gpu(data.f_eq_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_7_gpu(data.f_eq_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_8_gpu(data.f_eq_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_0_gpu(data.f_temp_0->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_1_gpu(data.f_temp_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_2_gpu(data.f_temp_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_3_gpu(data.f_temp_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_4_gpu(data.f_temp_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_5_gpu(data.f_temp_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_6_gpu(data.f_temp_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_7_gpu(data.f_temp_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_8_gpu(data.f_temp_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

                cuda_boundary_init_fsi_dir_1_float(
                        cuda_dir_7_gpu,
                        boundary_flags_gpu,
                        solid_flags_gpu,
                        solid_old_flags_gpu,
                        h_gpu,
                        u_gpu,
                        v_gpu,
                        f_0_gpu,
                        f_1_gpu,
                        f_2_gpu,
                        f_3_gpu,
                        f_4_gpu,
                        f_5_gpu,
                        f_6_gpu,
                        f_7_gpu,
                        f_8_gpu,
                        f_eq_0_gpu,
                        f_eq_1_gpu,
                        f_eq_2_gpu,
                        f_eq_3_gpu,
                        f_eq_4_gpu,
                        f_eq_5_gpu,
                        f_eq_6_gpu,
                        f_eq_7_gpu,
                        f_eq_8_gpu,
                        f_temp_0_gpu,
                        f_temp_1_gpu,
                        f_temp_2_gpu,
                        f_temp_3_gpu,
                        f_temp_4_gpu,
                        f_temp_5_gpu,
                        f_temp_6_gpu,
                        f_temp_7_gpu,
                        f_temp_8_gpu,
                        solids.current_u,
                        solids.current_v,
                        data.h->size(),
                        blocksize);

                info.cuda_dir_7->unlock(lm_read_only);

                solids.boundary_flags->unlock(lm_read_only);
                solids.solid_flags->unlock(lm_read_only);

                solids.solid_old_flags->unlock(lm_read_and_write);

                data.h->unlock(lm_read_and_write);
                data.u->unlock(lm_read_and_write);
                data.v->unlock(lm_read_and_write);

                data.f_0->unlock(lm_read_and_write);
                data.f_1->unlock(lm_read_and_write);
                data.f_2->unlock(lm_read_and_write);
                data.f_3->unlock(lm_read_and_write);
                data.f_4->unlock(lm_read_and_write);
                data.f_5->unlock(lm_read_and_write);
                data.f_6->unlock(lm_read_and_write);
                data.f_7->unlock(lm_read_and_write);
                data.f_8->unlock(lm_read_and_write);
                data.f_eq_0->unlock(lm_read_and_write);
                data.f_eq_1->unlock(lm_read_and_write);
                data.f_eq_2->unlock(lm_read_and_write);
                data.f_eq_3->unlock(lm_read_and_write);
                data.f_eq_4->unlock(lm_read_and_write);
                data.f_eq_5->unlock(lm_read_and_write);
                data.f_eq_6->unlock(lm_read_and_write);
                data.f_eq_7->unlock(lm_read_and_write);
                data.f_eq_8->unlock(lm_read_and_write);
                data.f_temp_0->unlock(lm_read_and_write);
                data.f_temp_1->unlock(lm_read_and_write);
                data.f_temp_2->unlock(lm_read_and_write);
                data.f_temp_3->unlock(lm_read_and_write);
                data.f_temp_4->unlock(lm_read_and_write);
                data.f_temp_5->unlock(lm_read_and_write);
                data.f_temp_6->unlock(lm_read_and_write);
                data.f_temp_7->unlock(lm_read_and_write);
                data.f_temp_8->unlock(lm_read_and_write);
            }
    };

    class cudaBoundaryInitFSIGridfloatDir4
    {
        private:
            PackedGridInfo<D2Q9> & info;
            PackedGridData<D2Q9, float> & data;
            PackedSolidData<lbm_lattice_types::D2Q9, float> & solids;
            unsigned long blocksize;
        public:
            cudaBoundaryInitFSIGridfloatDir4(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data, PackedSolidData<lbm_lattice_types::D2Q9, float> & solids, unsigned long blocksize) :
                info(info),
                data(data),
                solids(solids),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * cuda_dir_8_gpu(info.cuda_dir_8->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * h_gpu(data.h->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * u_gpu(data.u->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * v_gpu(data.v->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

                void * boundary_flags_gpu(solids.boundary_flags->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * solid_flags_gpu(solids.solid_flags->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * solid_old_flags_gpu(solids.solid_old_flags->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

                void * f_0_gpu(data.f_0->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_1_gpu(data.f_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_2_gpu(data.f_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_3_gpu(data.f_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_4_gpu(data.f_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_5_gpu(data.f_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_6_gpu(data.f_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_7_gpu(data.f_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_8_gpu(data.f_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_0_gpu(data.f_eq_0->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_1_gpu(data.f_eq_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_2_gpu(data.f_eq_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_3_gpu(data.f_eq_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_4_gpu(data.f_eq_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_5_gpu(data.f_eq_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_6_gpu(data.f_eq_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_7_gpu(data.f_eq_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_8_gpu(data.f_eq_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_0_gpu(data.f_temp_0->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_1_gpu(data.f_temp_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_2_gpu(data.f_temp_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_3_gpu(data.f_temp_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_4_gpu(data.f_temp_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_5_gpu(data.f_temp_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_6_gpu(data.f_temp_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_7_gpu(data.f_temp_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_8_gpu(data.f_temp_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

                cuda_boundary_init_fsi_dir_1_float(
                        cuda_dir_8_gpu,
                        boundary_flags_gpu,
                        solid_flags_gpu,
                        solid_old_flags_gpu,
                        h_gpu,
                        u_gpu,
                        v_gpu,
                        f_0_gpu,
                        f_1_gpu,
                        f_2_gpu,
                        f_3_gpu,
                        f_4_gpu,
                        f_5_gpu,
                        f_6_gpu,
                        f_7_gpu,
                        f_8_gpu,
                        f_eq_0_gpu,
                        f_eq_1_gpu,
                        f_eq_2_gpu,
                        f_eq_3_gpu,
                        f_eq_4_gpu,
                        f_eq_5_gpu,
                        f_eq_6_gpu,
                        f_eq_7_gpu,
                        f_eq_8_gpu,
                        f_temp_0_gpu,
                        f_temp_1_gpu,
                        f_temp_2_gpu,
                        f_temp_3_gpu,
                        f_temp_4_gpu,
                        f_temp_5_gpu,
                        f_temp_6_gpu,
                        f_temp_7_gpu,
                        f_temp_8_gpu,
                        solids.current_u,
                        solids.current_v,
                        data.h->size(),
                        blocksize);

                info.cuda_dir_8->unlock(lm_read_only);

                solids.boundary_flags->unlock(lm_read_only);
                solids.solid_flags->unlock(lm_read_only);

                solids.solid_old_flags->unlock(lm_read_and_write);

                data.h->unlock(lm_read_and_write);
                data.u->unlock(lm_read_and_write);
                data.v->unlock(lm_read_and_write);

                data.f_0->unlock(lm_read_and_write);
                data.f_1->unlock(lm_read_and_write);
                data.f_2->unlock(lm_read_and_write);
                data.f_3->unlock(lm_read_and_write);
                data.f_4->unlock(lm_read_and_write);
                data.f_5->unlock(lm_read_and_write);
                data.f_6->unlock(lm_read_and_write);
                data.f_7->unlock(lm_read_and_write);
                data.f_8->unlock(lm_read_and_write);
                data.f_eq_0->unlock(lm_read_and_write);
                data.f_eq_1->unlock(lm_read_and_write);
                data.f_eq_2->unlock(lm_read_and_write);
                data.f_eq_3->unlock(lm_read_and_write);
                data.f_eq_4->unlock(lm_read_and_write);
                data.f_eq_5->unlock(lm_read_and_write);
                data.f_eq_6->unlock(lm_read_and_write);
                data.f_eq_7->unlock(lm_read_and_write);
                data.f_eq_8->unlock(lm_read_and_write);
                data.f_temp_0->unlock(lm_read_and_write);
                data.f_temp_1->unlock(lm_read_and_write);
                data.f_temp_2->unlock(lm_read_and_write);
                data.f_temp_3->unlock(lm_read_and_write);
                data.f_temp_4->unlock(lm_read_and_write);
                data.f_temp_5->unlock(lm_read_and_write);
                data.f_temp_6->unlock(lm_read_and_write);
                data.f_temp_7->unlock(lm_read_and_write);
                data.f_temp_8->unlock(lm_read_and_write);
            }
    };

    class cudaBoundaryInitFSIGridfloatDir5
    {
        private:
            PackedGridInfo<D2Q9> & info;
            PackedGridData<D2Q9, float> & data;
            PackedSolidData<lbm_lattice_types::D2Q9, float> & solids;
            unsigned long blocksize;
        public:
            cudaBoundaryInitFSIGridfloatDir5(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data, PackedSolidData<lbm_lattice_types::D2Q9, float> & solids, unsigned long blocksize) :
                info(info),
                data(data),
                solids(solids),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * cuda_dir_1_gpu(info.cuda_dir_1->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * h_gpu(data.h->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * u_gpu(data.u->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * v_gpu(data.v->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

                void * boundary_flags_gpu(solids.boundary_flags->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * solid_flags_gpu(solids.solid_flags->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * solid_old_flags_gpu(solids.solid_old_flags->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

                void * f_0_gpu(data.f_0->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_1_gpu(data.f_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_2_gpu(data.f_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_3_gpu(data.f_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_4_gpu(data.f_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_5_gpu(data.f_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_6_gpu(data.f_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_7_gpu(data.f_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_8_gpu(data.f_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_0_gpu(data.f_eq_0->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_1_gpu(data.f_eq_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_2_gpu(data.f_eq_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_3_gpu(data.f_eq_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_4_gpu(data.f_eq_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_5_gpu(data.f_eq_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_6_gpu(data.f_eq_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_7_gpu(data.f_eq_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_8_gpu(data.f_eq_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_0_gpu(data.f_temp_0->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_1_gpu(data.f_temp_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_2_gpu(data.f_temp_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_3_gpu(data.f_temp_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_4_gpu(data.f_temp_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_5_gpu(data.f_temp_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_6_gpu(data.f_temp_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_7_gpu(data.f_temp_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_8_gpu(data.f_temp_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

                cuda_boundary_init_fsi_dir_1_float(
                        cuda_dir_1_gpu,
                        boundary_flags_gpu,
                        solid_flags_gpu,
                        solid_old_flags_gpu,
                        h_gpu,
                        u_gpu,
                        v_gpu,
                        f_0_gpu,
                        f_1_gpu,
                        f_2_gpu,
                        f_3_gpu,
                        f_4_gpu,
                        f_5_gpu,
                        f_6_gpu,
                        f_7_gpu,
                        f_8_gpu,
                        f_eq_0_gpu,
                        f_eq_1_gpu,
                        f_eq_2_gpu,
                        f_eq_3_gpu,
                        f_eq_4_gpu,
                        f_eq_5_gpu,
                        f_eq_6_gpu,
                        f_eq_7_gpu,
                        f_eq_8_gpu,
                        f_temp_0_gpu,
                        f_temp_1_gpu,
                        f_temp_2_gpu,
                        f_temp_3_gpu,
                        f_temp_4_gpu,
                        f_temp_5_gpu,
                        f_temp_6_gpu,
                        f_temp_7_gpu,
                        f_temp_8_gpu,
                        solids.current_u,
                        solids.current_v,
                        data.h->size(),
                        blocksize);

                info.cuda_dir_1->unlock(lm_read_only);

                solids.boundary_flags->unlock(lm_read_only);
                solids.solid_flags->unlock(lm_read_only);

                solids.solid_old_flags->unlock(lm_read_and_write);

                data.h->unlock(lm_read_and_write);
                data.u->unlock(lm_read_and_write);
                data.v->unlock(lm_read_and_write);

                data.f_0->unlock(lm_read_and_write);
                data.f_1->unlock(lm_read_and_write);
                data.f_2->unlock(lm_read_and_write);
                data.f_3->unlock(lm_read_and_write);
                data.f_4->unlock(lm_read_and_write);
                data.f_5->unlock(lm_read_and_write);
                data.f_6->unlock(lm_read_and_write);
                data.f_7->unlock(lm_read_and_write);
                data.f_8->unlock(lm_read_and_write);
                data.f_eq_0->unlock(lm_read_and_write);
                data.f_eq_1->unlock(lm_read_and_write);
                data.f_eq_2->unlock(lm_read_and_write);
                data.f_eq_3->unlock(lm_read_and_write);
                data.f_eq_4->unlock(lm_read_and_write);
                data.f_eq_5->unlock(lm_read_and_write);
                data.f_eq_6->unlock(lm_read_and_write);
                data.f_eq_7->unlock(lm_read_and_write);
                data.f_eq_8->unlock(lm_read_and_write);
                data.f_temp_0->unlock(lm_read_and_write);
                data.f_temp_1->unlock(lm_read_and_write);
                data.f_temp_2->unlock(lm_read_and_write);
                data.f_temp_3->unlock(lm_read_and_write);
                data.f_temp_4->unlock(lm_read_and_write);
                data.f_temp_5->unlock(lm_read_and_write);
                data.f_temp_6->unlock(lm_read_and_write);
                data.f_temp_7->unlock(lm_read_and_write);
                data.f_temp_8->unlock(lm_read_and_write);
            }
    };

    class cudaBoundaryInitFSIGridfloatDir6
    {
        private:
            PackedGridInfo<D2Q9> & info;
            PackedGridData<D2Q9, float> & data;
            PackedSolidData<lbm_lattice_types::D2Q9, float> & solids;
            unsigned long blocksize;
        public:
            cudaBoundaryInitFSIGridfloatDir6(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data, PackedSolidData<lbm_lattice_types::D2Q9, float> & solids, unsigned long blocksize) :
                info(info),
                data(data),
                solids(solids),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * cuda_dir_2_gpu(info.cuda_dir_2->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * h_gpu(data.h->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * u_gpu(data.u->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * v_gpu(data.v->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

                void * boundary_flags_gpu(solids.boundary_flags->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * solid_flags_gpu(solids.solid_flags->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * solid_old_flags_gpu(solids.solid_old_flags->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

                void * f_0_gpu(data.f_0->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_1_gpu(data.f_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_2_gpu(data.f_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_3_gpu(data.f_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_4_gpu(data.f_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_5_gpu(data.f_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_6_gpu(data.f_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_7_gpu(data.f_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_8_gpu(data.f_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_0_gpu(data.f_eq_0->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_1_gpu(data.f_eq_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_2_gpu(data.f_eq_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_3_gpu(data.f_eq_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_4_gpu(data.f_eq_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_5_gpu(data.f_eq_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_6_gpu(data.f_eq_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_7_gpu(data.f_eq_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_8_gpu(data.f_eq_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_0_gpu(data.f_temp_0->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_1_gpu(data.f_temp_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_2_gpu(data.f_temp_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_3_gpu(data.f_temp_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_4_gpu(data.f_temp_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_5_gpu(data.f_temp_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_6_gpu(data.f_temp_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_7_gpu(data.f_temp_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_8_gpu(data.f_temp_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

                cuda_boundary_init_fsi_dir_1_float(
                        cuda_dir_2_gpu,
                        boundary_flags_gpu,
                        solid_flags_gpu,
                        solid_old_flags_gpu,
                        h_gpu,
                        u_gpu,
                        v_gpu,
                        f_0_gpu,
                        f_1_gpu,
                        f_2_gpu,
                        f_3_gpu,
                        f_4_gpu,
                        f_5_gpu,
                        f_6_gpu,
                        f_7_gpu,
                        f_8_gpu,
                        f_eq_0_gpu,
                        f_eq_1_gpu,
                        f_eq_2_gpu,
                        f_eq_3_gpu,
                        f_eq_4_gpu,
                        f_eq_5_gpu,
                        f_eq_6_gpu,
                        f_eq_7_gpu,
                        f_eq_8_gpu,
                        f_temp_0_gpu,
                        f_temp_1_gpu,
                        f_temp_2_gpu,
                        f_temp_3_gpu,
                        f_temp_4_gpu,
                        f_temp_5_gpu,
                        f_temp_6_gpu,
                        f_temp_7_gpu,
                        f_temp_8_gpu,
                        solids.current_u,
                        solids.current_v,
                        data.h->size(),
                        blocksize);

                info.cuda_dir_2->unlock(lm_read_only);

                solids.boundary_flags->unlock(lm_read_only);
                solids.solid_flags->unlock(lm_read_only);

                solids.solid_old_flags->unlock(lm_read_and_write);

                data.h->unlock(lm_read_and_write);
                data.u->unlock(lm_read_and_write);
                data.v->unlock(lm_read_and_write);

                data.f_0->unlock(lm_read_and_write);
                data.f_1->unlock(lm_read_and_write);
                data.f_2->unlock(lm_read_and_write);
                data.f_3->unlock(lm_read_and_write);
                data.f_4->unlock(lm_read_and_write);
                data.f_5->unlock(lm_read_and_write);
                data.f_6->unlock(lm_read_and_write);
                data.f_7->unlock(lm_read_and_write);
                data.f_8->unlock(lm_read_and_write);
                data.f_eq_0->unlock(lm_read_and_write);
                data.f_eq_1->unlock(lm_read_and_write);
                data.f_eq_2->unlock(lm_read_and_write);
                data.f_eq_3->unlock(lm_read_and_write);
                data.f_eq_4->unlock(lm_read_and_write);
                data.f_eq_5->unlock(lm_read_and_write);
                data.f_eq_6->unlock(lm_read_and_write);
                data.f_eq_7->unlock(lm_read_and_write);
                data.f_eq_8->unlock(lm_read_and_write);
                data.f_temp_0->unlock(lm_read_and_write);
                data.f_temp_1->unlock(lm_read_and_write);
                data.f_temp_2->unlock(lm_read_and_write);
                data.f_temp_3->unlock(lm_read_and_write);
                data.f_temp_4->unlock(lm_read_and_write);
                data.f_temp_5->unlock(lm_read_and_write);
                data.f_temp_6->unlock(lm_read_and_write);
                data.f_temp_7->unlock(lm_read_and_write);
                data.f_temp_8->unlock(lm_read_and_write);
            }
    };

    class cudaBoundaryInitFSIGridfloatDir7
    {
        private:
            PackedGridInfo<D2Q9> & info;
            PackedGridData<D2Q9, float> & data;
            PackedSolidData<lbm_lattice_types::D2Q9, float> & solids;
            unsigned long blocksize;
        public:
            cudaBoundaryInitFSIGridfloatDir7(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data, PackedSolidData<lbm_lattice_types::D2Q9, float> & solids, unsigned long blocksize) :
                info(info),
                data(data),
                solids(solids),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * cuda_dir_3_gpu(info.cuda_dir_3->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * h_gpu(data.h->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * u_gpu(data.u->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * v_gpu(data.v->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

                void * boundary_flags_gpu(solids.boundary_flags->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * solid_flags_gpu(solids.solid_flags->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * solid_old_flags_gpu(solids.solid_old_flags->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

                void * f_0_gpu(data.f_0->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_1_gpu(data.f_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_2_gpu(data.f_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_3_gpu(data.f_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_4_gpu(data.f_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_5_gpu(data.f_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_6_gpu(data.f_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_7_gpu(data.f_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_8_gpu(data.f_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_0_gpu(data.f_eq_0->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_1_gpu(data.f_eq_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_2_gpu(data.f_eq_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_3_gpu(data.f_eq_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_4_gpu(data.f_eq_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_5_gpu(data.f_eq_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_6_gpu(data.f_eq_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_7_gpu(data.f_eq_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_8_gpu(data.f_eq_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_0_gpu(data.f_temp_0->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_1_gpu(data.f_temp_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_2_gpu(data.f_temp_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_3_gpu(data.f_temp_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_4_gpu(data.f_temp_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_5_gpu(data.f_temp_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_6_gpu(data.f_temp_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_7_gpu(data.f_temp_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_8_gpu(data.f_temp_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

                cuda_boundary_init_fsi_dir_1_float(
                        cuda_dir_3_gpu,
                        boundary_flags_gpu,
                        solid_flags_gpu,
                        solid_old_flags_gpu,
                        h_gpu,
                        u_gpu,
                        v_gpu,
                        f_0_gpu,
                        f_1_gpu,
                        f_2_gpu,
                        f_3_gpu,
                        f_4_gpu,
                        f_5_gpu,
                        f_6_gpu,
                        f_7_gpu,
                        f_8_gpu,
                        f_eq_0_gpu,
                        f_eq_1_gpu,
                        f_eq_2_gpu,
                        f_eq_3_gpu,
                        f_eq_4_gpu,
                        f_eq_5_gpu,
                        f_eq_6_gpu,
                        f_eq_7_gpu,
                        f_eq_8_gpu,
                        f_temp_0_gpu,
                        f_temp_1_gpu,
                        f_temp_2_gpu,
                        f_temp_3_gpu,
                        f_temp_4_gpu,
                        f_temp_5_gpu,
                        f_temp_6_gpu,
                        f_temp_7_gpu,
                        f_temp_8_gpu,
                        solids.current_u,
                        solids.current_v,
                        data.h->size(),
                        blocksize);

                info.cuda_dir_3->unlock(lm_read_only);

                solids.boundary_flags->unlock(lm_read_only);
                solids.solid_flags->unlock(lm_read_only);

                solids.solid_old_flags->unlock(lm_read_and_write);

                data.h->unlock(lm_read_and_write);
                data.u->unlock(lm_read_and_write);
                data.v->unlock(lm_read_and_write);

                data.f_0->unlock(lm_read_and_write);
                data.f_1->unlock(lm_read_and_write);
                data.f_2->unlock(lm_read_and_write);
                data.f_3->unlock(lm_read_and_write);
                data.f_4->unlock(lm_read_and_write);
                data.f_5->unlock(lm_read_and_write);
                data.f_6->unlock(lm_read_and_write);
                data.f_7->unlock(lm_read_and_write);
                data.f_8->unlock(lm_read_and_write);
                data.f_eq_0->unlock(lm_read_and_write);
                data.f_eq_1->unlock(lm_read_and_write);
                data.f_eq_2->unlock(lm_read_and_write);
                data.f_eq_3->unlock(lm_read_and_write);
                data.f_eq_4->unlock(lm_read_and_write);
                data.f_eq_5->unlock(lm_read_and_write);
                data.f_eq_6->unlock(lm_read_and_write);
                data.f_eq_7->unlock(lm_read_and_write);
                data.f_eq_8->unlock(lm_read_and_write);
                data.f_temp_0->unlock(lm_read_and_write);
                data.f_temp_1->unlock(lm_read_and_write);
                data.f_temp_2->unlock(lm_read_and_write);
                data.f_temp_3->unlock(lm_read_and_write);
                data.f_temp_4->unlock(lm_read_and_write);
                data.f_temp_5->unlock(lm_read_and_write);
                data.f_temp_6->unlock(lm_read_and_write);
                data.f_temp_7->unlock(lm_read_and_write);
                data.f_temp_8->unlock(lm_read_and_write);
            }
    };

    class cudaBoundaryInitFSIGridfloatDir8
    {
        private:
            PackedGridInfo<D2Q9> & info;
            PackedGridData<D2Q9, float> & data;
            PackedSolidData<lbm_lattice_types::D2Q9, float> & solids;
            unsigned long blocksize;
        public:
            cudaBoundaryInitFSIGridfloatDir8(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data, PackedSolidData<lbm_lattice_types::D2Q9, float> & solids, unsigned long blocksize) :
                info(info),
                data(data),
                solids(solids),
                blocksize(blocksize)
            {
            }

            void operator() ()
            {
                void * cuda_dir_4_gpu(info.cuda_dir_4->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * h_gpu(data.h->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * u_gpu(data.u->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * v_gpu(data.v->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

                void * boundary_flags_gpu(solids.boundary_flags->lock(lm_read_only, tags::GPU::CUDA::memory_value));
                void * solid_flags_gpu(solids.solid_flags->lock(lm_read_only, tags::GPU::CUDA::memory_value));

                void * solid_old_flags_gpu(solids.solid_old_flags->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

                void * f_0_gpu(data.f_0->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_1_gpu(data.f_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_2_gpu(data.f_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_3_gpu(data.f_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_4_gpu(data.f_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_5_gpu(data.f_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_6_gpu(data.f_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_7_gpu(data.f_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_8_gpu(data.f_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_0_gpu(data.f_eq_0->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_1_gpu(data.f_eq_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_2_gpu(data.f_eq_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_3_gpu(data.f_eq_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_4_gpu(data.f_eq_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_5_gpu(data.f_eq_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_6_gpu(data.f_eq_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_7_gpu(data.f_eq_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_eq_8_gpu(data.f_eq_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_0_gpu(data.f_temp_0->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_1_gpu(data.f_temp_1->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_2_gpu(data.f_temp_2->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_3_gpu(data.f_temp_3->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_4_gpu(data.f_temp_4->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_5_gpu(data.f_temp_5->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_6_gpu(data.f_temp_6->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_7_gpu(data.f_temp_7->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));
                void * f_temp_8_gpu(data.f_temp_8->lock(lm_read_and_write, tags::GPU::CUDA::memory_value));

                cuda_boundary_init_fsi_dir_1_float(
                        cuda_dir_4_gpu,
                        boundary_flags_gpu,
                        solid_flags_gpu,
                        solid_old_flags_gpu,
                        h_gpu,
                        u_gpu,
                        v_gpu,
                        f_0_gpu,
                        f_1_gpu,
                        f_2_gpu,
                        f_3_gpu,
                        f_4_gpu,
                        f_5_gpu,
                        f_6_gpu,
                        f_7_gpu,
                        f_8_gpu,
                        f_eq_0_gpu,
                        f_eq_1_gpu,
                        f_eq_2_gpu,
                        f_eq_3_gpu,
                        f_eq_4_gpu,
                        f_eq_5_gpu,
                        f_eq_6_gpu,
                        f_eq_7_gpu,
                        f_eq_8_gpu,
                        f_temp_0_gpu,
                        f_temp_1_gpu,
                        f_temp_2_gpu,
                        f_temp_3_gpu,
                        f_temp_4_gpu,
                        f_temp_5_gpu,
                        f_temp_6_gpu,
                        f_temp_7_gpu,
                        f_temp_8_gpu,
                        solids.current_u,
                        solids.current_v,
                        data.h->size(),
                        blocksize);

                info.cuda_dir_4->unlock(lm_read_only);

                solids.boundary_flags->unlock(lm_read_only);
                solids.solid_flags->unlock(lm_read_only);

                solids.solid_old_flags->unlock(lm_read_and_write);

                data.h->unlock(lm_read_and_write);
                data.u->unlock(lm_read_and_write);
                data.v->unlock(lm_read_and_write);

                data.f_0->unlock(lm_read_and_write);
                data.f_1->unlock(lm_read_and_write);
                data.f_2->unlock(lm_read_and_write);
                data.f_3->unlock(lm_read_and_write);
                data.f_4->unlock(lm_read_and_write);
                data.f_5->unlock(lm_read_and_write);
                data.f_6->unlock(lm_read_and_write);
                data.f_7->unlock(lm_read_and_write);
                data.f_8->unlock(lm_read_and_write);
                data.f_eq_0->unlock(lm_read_and_write);
                data.f_eq_1->unlock(lm_read_and_write);
                data.f_eq_2->unlock(lm_read_and_write);
                data.f_eq_3->unlock(lm_read_and_write);
                data.f_eq_4->unlock(lm_read_and_write);
                data.f_eq_5->unlock(lm_read_and_write);
                data.f_eq_6->unlock(lm_read_and_write);
                data.f_eq_7->unlock(lm_read_and_write);
                data.f_eq_8->unlock(lm_read_and_write);
                data.f_temp_0->unlock(lm_read_and_write);
                data.f_temp_1->unlock(lm_read_and_write);
                data.f_temp_2->unlock(lm_read_and_write);
                data.f_temp_3->unlock(lm_read_and_write);
                data.f_temp_4->unlock(lm_read_and_write);
                data.f_temp_5->unlock(lm_read_and_write);
                data.f_temp_6->unlock(lm_read_and_write);
                data.f_temp_7->unlock(lm_read_and_write);
                data.f_temp_8->unlock(lm_read_and_write);
            }
    };
}

void BoundaryInitFSI<tags::GPU::CUDA, D2Q9::DIR_1>::value(
        PackedGridInfo<lbm_lattice_types::D2Q9> & info,
        PackedGridData<lbm_lattice_types::D2Q9, float> & data,
        PackedSolidData<lbm_lattice_types::D2Q9, float> & solids)
{
    CONTEXT("When performing boundary initialization FSI correction with object moving along 1 (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::collide_stream_grid_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaBoundaryInitFSIGridfloatDir1 task(info, data, solids, blocksize);
        task();
    }
    else
    {
        cudaBoundaryInitFSIGridfloatDir1 task(info, data, solids, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }
}

void BoundaryInitFSI<tags::GPU::CUDA, D2Q9::DIR_2>::value(
        PackedGridInfo<lbm_lattice_types::D2Q9> & info,
        PackedGridData<lbm_lattice_types::D2Q9, float> & data,
        PackedSolidData<lbm_lattice_types::D2Q9, float> & solids)
{
    CONTEXT("When performing boundary initialization FSI correction with object moving along 2 (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::collide_stream_grid_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaBoundaryInitFSIGridfloatDir2 task(info, data, solids, blocksize);
        task();
    }
    else
    {
        cudaBoundaryInitFSIGridfloatDir2 task(info, data, solids, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }
}

void BoundaryInitFSI<tags::GPU::CUDA, D2Q9::DIR_3>::value(
        PackedGridInfo<lbm_lattice_types::D2Q9> & info,
        PackedGridData<lbm_lattice_types::D2Q9, float> & data,
        PackedSolidData<lbm_lattice_types::D2Q9, float> & solids)
{
    CONTEXT("When performing boundary initialization FSI correction with object moving along 3 (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::collide_stream_grid_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaBoundaryInitFSIGridfloatDir3 task(info, data, solids, blocksize);
        task();
    }
    else
    {
        cudaBoundaryInitFSIGridfloatDir3 task(info, data, solids, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }
}

void BoundaryInitFSI<tags::GPU::CUDA, D2Q9::DIR_4>::value(
        PackedGridInfo<lbm_lattice_types::D2Q9> & info,
        PackedGridData<lbm_lattice_types::D2Q9, float> & data,
        PackedSolidData<lbm_lattice_types::D2Q9, float> & solids)
{
    CONTEXT("When performing boundary initialization FSI correction with object moving along 4 (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::collide_stream_grid_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaBoundaryInitFSIGridfloatDir4 task(info, data, solids, blocksize);
        task();
    }
    else
    {
        cudaBoundaryInitFSIGridfloatDir4 task(info, data, solids, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }
}

void BoundaryInitFSI<tags::GPU::CUDA, D2Q9::DIR_5>::value(
        PackedGridInfo<lbm_lattice_types::D2Q9> & info,
        PackedGridData<lbm_lattice_types::D2Q9, float> & data,
        PackedSolidData<lbm_lattice_types::D2Q9, float> & solids)
{
    CONTEXT("When performing boundary initialization FSI correction with object moving along 5 (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::collide_stream_grid_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaBoundaryInitFSIGridfloatDir5 task(info, data, solids, blocksize);
        task();
    }
    else
    {
        cudaBoundaryInitFSIGridfloatDir5 task(info, data, solids, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }
}

void BoundaryInitFSI<tags::GPU::CUDA, D2Q9::DIR_6>::value(
        PackedGridInfo<lbm_lattice_types::D2Q9> & info,
        PackedGridData<lbm_lattice_types::D2Q9, float> & data,
        PackedSolidData<lbm_lattice_types::D2Q9, float> & solids)
{
    CONTEXT("When performing boundary initialization FSI correction with object moving along 6 (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::collide_stream_grid_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaBoundaryInitFSIGridfloatDir6 task(info, data, solids, blocksize);
        task();
    }
    else
    {
        cudaBoundaryInitFSIGridfloatDir6 task(info, data, solids, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }
}

void BoundaryInitFSI<tags::GPU::CUDA, D2Q9::DIR_7>::value(
        PackedGridInfo<lbm_lattice_types::D2Q9> & info,
        PackedGridData<lbm_lattice_types::D2Q9, float> & data,
        PackedSolidData<lbm_lattice_types::D2Q9, float> & solids)
{
    CONTEXT("When performing boundary initialization FSI correction with object moving along 7 (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::collide_stream_grid_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaBoundaryInitFSIGridfloatDir7 task(info, data, solids, blocksize);
        task();
    }
    else
    {
        cudaBoundaryInitFSIGridfloatDir7 task(info, data, solids, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }
}

void BoundaryInitFSI<tags::GPU::CUDA, D2Q9::DIR_8>::value(
        PackedGridInfo<lbm_lattice_types::D2Q9> & info,
        PackedGridData<lbm_lattice_types::D2Q9, float> & data,
        PackedSolidData<lbm_lattice_types::D2Q9, float> & solids)
{
    CONTEXT("When performing boundary initialization FSI correction with object moving along 8 (CUDA):");

    unsigned long blocksize(Configuration::instance()->get_value("cuda::collide_stream_grid_float", 128ul));

    if (! cuda::GPUPool::instance()->idle())
    {
        cudaBoundaryInitFSIGridfloatDir4 task(info, data, solids, blocksize);
        task();
    }
    else
    {
        cudaBoundaryInitFSIGridfloatDir4 task(info, data, solids, blocksize);
        cuda::GPUPool::instance()->enqueue(task, 0).wait();
    }
}

