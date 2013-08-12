/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2013 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
 *
 * This file is part of the LA C++ library. LibLa is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibLa is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/lbm/collide_stream_grid.hh>
#include <honei/backends/opencl/operations.hh>
#include <honei/util/profiler.hh>
#include <typeinfo>

namespace honei
{
    namespace opencl
    {
        template <typename Tag_, typename DT_>
        void common_collide_stream_grid(PackedGridInfo<lbm_lattice_types::D2Q9> & info,
             PackedGridData<lbm_lattice_types::D2Q9, DT_> & data, DT_ tau)
        {
            info.limits->lock(lm_read_only);

            void * cuda_dir_1_cl(info.cuda_dir_1->lock(lm_read_only, Tag_::memory_value));
            void * cuda_dir_2_cl(info.cuda_dir_2->lock(lm_read_only, Tag_::memory_value));
            void * cuda_dir_3_cl(info.cuda_dir_3->lock(lm_read_only, Tag_::memory_value));
            void * cuda_dir_4_cl(info.cuda_dir_4->lock(lm_read_only, Tag_::memory_value));
            void * cuda_dir_5_cl(info.cuda_dir_5->lock(lm_read_only, Tag_::memory_value));
            void * cuda_dir_6_cl(info.cuda_dir_6->lock(lm_read_only, Tag_::memory_value));
            void * cuda_dir_7_cl(info.cuda_dir_7->lock(lm_read_only, Tag_::memory_value));
            void * cuda_dir_8_cl(info.cuda_dir_8->lock(lm_read_only, Tag_::memory_value));

            void * f_eq_0_cl(data.f_eq_0->lock(lm_read_only, Tag_::memory_value));
            void * f_eq_1_cl(data.f_eq_1->lock(lm_read_only, Tag_::memory_value));
            void * f_eq_2_cl(data.f_eq_2->lock(lm_read_only, Tag_::memory_value));
            void * f_eq_3_cl(data.f_eq_3->lock(lm_read_only, Tag_::memory_value));
            void * f_eq_4_cl(data.f_eq_4->lock(lm_read_only, Tag_::memory_value));
            void * f_eq_5_cl(data.f_eq_5->lock(lm_read_only, Tag_::memory_value));
            void * f_eq_6_cl(data.f_eq_6->lock(lm_read_only, Tag_::memory_value));
            void * f_eq_7_cl(data.f_eq_7->lock(lm_read_only, Tag_::memory_value));
            void * f_eq_8_cl(data.f_eq_8->lock(lm_read_only, Tag_::memory_value));

            void * f_0_cl(data.f_0->lock(lm_read_only, Tag_::memory_value));
            void * f_1_cl(data.f_1->lock(lm_read_only, Tag_::memory_value));
            void * f_2_cl(data.f_2->lock(lm_read_only, Tag_::memory_value));
            void * f_3_cl(data.f_3->lock(lm_read_only, Tag_::memory_value));
            void * f_4_cl(data.f_4->lock(lm_read_only, Tag_::memory_value));
            void * f_5_cl(data.f_5->lock(lm_read_only, Tag_::memory_value));
            void * f_6_cl(data.f_6->lock(lm_read_only, Tag_::memory_value));
            void * f_7_cl(data.f_7->lock(lm_read_only, Tag_::memory_value));
            void * f_8_cl(data.f_8->lock(lm_read_only, Tag_::memory_value));

            void * f_temp_0_cl(data.f_temp_0->lock(lm_write_only, Tag_::memory_value));
            void * f_temp_1_cl(data.f_temp_1->lock(lm_write_only, Tag_::memory_value));
            void * f_temp_2_cl(data.f_temp_2->lock(lm_write_only, Tag_::memory_value));
            void * f_temp_3_cl(data.f_temp_3->lock(lm_write_only, Tag_::memory_value));
            void * f_temp_4_cl(data.f_temp_4->lock(lm_write_only, Tag_::memory_value));
            void * f_temp_5_cl(data.f_temp_5->lock(lm_write_only, Tag_::memory_value));
            void * f_temp_6_cl(data.f_temp_6->lock(lm_write_only, Tag_::memory_value));
            void * f_temp_7_cl(data.f_temp_7->lock(lm_write_only, Tag_::memory_value));
            void * f_temp_8_cl(data.f_temp_8->lock(lm_write_only, Tag_::memory_value));


            unsigned long start((*info.limits)[0]);
            unsigned long end((*info.limits)[info.limits->size() - 1]);

            std::string opname("is_not_used");
            opname += typeid(DT_).name();
            opencl::collide_stream_grid(start, end,
                    cuda_dir_1_cl, cuda_dir_2_cl, cuda_dir_3_cl, cuda_dir_4_cl,
                    cuda_dir_5_cl, cuda_dir_6_cl, cuda_dir_7_cl, cuda_dir_8_cl,
                    f_eq_0_cl, f_eq_1_cl, f_eq_2_cl,
                    f_eq_3_cl, f_eq_4_cl, f_eq_5_cl,
                    f_eq_6_cl, f_eq_7_cl, f_eq_8_cl,
                    f_0_cl, f_1_cl, f_2_cl,
                    f_3_cl, f_4_cl, f_5_cl,
                    f_6_cl, f_7_cl, f_8_cl,
                    f_temp_0_cl, f_temp_1_cl, f_temp_2_cl,
                    f_temp_3_cl, f_temp_4_cl, f_temp_5_cl,
                    f_temp_6_cl, f_temp_7_cl, f_temp_8_cl,
                    tau, data.h->size(),
                    tag_to_device<Tag_>(), opname);

            info.limits->unlock(lm_read_only);

            info.cuda_dir_1->unlock(lm_read_only);
            info.cuda_dir_2->unlock(lm_read_only);
            info.cuda_dir_3->unlock(lm_read_only);
            info.cuda_dir_4->unlock(lm_read_only);
            info.cuda_dir_5->unlock(lm_read_only);
            info.cuda_dir_6->unlock(lm_read_only);
            info.cuda_dir_7->unlock(lm_read_only);
            info.cuda_dir_8->unlock(lm_read_only);

            data.f_eq_0->unlock(lm_read_only);
            data.f_eq_1->unlock(lm_read_only);
            data.f_eq_2->unlock(lm_read_only);
            data.f_eq_3->unlock(lm_read_only);
            data.f_eq_4->unlock(lm_read_only);
            data.f_eq_5->unlock(lm_read_only);
            data.f_eq_6->unlock(lm_read_only);
            data.f_eq_7->unlock(lm_read_only);
            data.f_eq_8->unlock(lm_read_only);

            data.f_0->unlock(lm_read_only);
            data.f_1->unlock(lm_read_only);
            data.f_2->unlock(lm_read_only);
            data.f_3->unlock(lm_read_only);
            data.f_4->unlock(lm_read_only);
            data.f_5->unlock(lm_read_only);
            data.f_6->unlock(lm_read_only);
            data.f_7->unlock(lm_read_only);
            data.f_8->unlock(lm_read_only);

            data.f_temp_0->unlock(lm_write_only);
            data.f_temp_1->unlock(lm_write_only);
            data.f_temp_2->unlock(lm_write_only);
            data.f_temp_3->unlock(lm_write_only);
            data.f_temp_4->unlock(lm_write_only);
            data.f_temp_5->unlock(lm_write_only);
            data.f_temp_6->unlock(lm_write_only);
            data.f_temp_7->unlock(lm_write_only);
            data.f_temp_8->unlock(lm_write_only);
        }
    }
}

using namespace honei;

template <typename DT_>
void CollideStreamGrid<tags::OpenCL::CPU, lbm_boundary_types::NOSLIP,
     lbm_lattice_types::D2Q9>::value(
             PackedGridInfo<lbm_lattice_types::D2Q9> & info,
             PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
             DT_ tau)
{
    CONTEXT("When performing collision and streaming (OpenCL CPU):");

    opencl::common_collide_stream_grid<tags::OpenCL::CPU>(info, data, tau);
}

template void CollideStreamGrid<tags::OpenCL::CPU, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>::value(PackedGridInfo<lbm_lattice_types::D2Q9> &, PackedGridData<lbm_lattice_types::D2Q9, float> &, float);
template void CollideStreamGrid<tags::OpenCL::CPU, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>::value(PackedGridInfo<lbm_lattice_types::D2Q9> &, PackedGridData<lbm_lattice_types::D2Q9, double> &, double);

template <typename DT_>
void CollideStreamGrid<tags::OpenCL::GPU, lbm_boundary_types::NOSLIP,
     lbm_lattice_types::D2Q9>::value(
             PackedGridInfo<lbm_lattice_types::D2Q9> & info,
             PackedGridData<lbm_lattice_types::D2Q9, DT_> & data,
             DT_ tau)
{
    CONTEXT("When performing collision and streaming (OpenCL GPU):");

    opencl::common_collide_stream_grid<tags::OpenCL::GPU>(info, data, tau);
}

template void CollideStreamGrid<tags::OpenCL::GPU, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>::value(PackedGridInfo<lbm_lattice_types::D2Q9> &, PackedGridData<lbm_lattice_types::D2Q9, float> &, float);
template void CollideStreamGrid<tags::OpenCL::GPU, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>::value(PackedGridInfo<lbm_lattice_types::D2Q9> &, PackedGridData<lbm_lattice_types::D2Q9, double> &, double);
