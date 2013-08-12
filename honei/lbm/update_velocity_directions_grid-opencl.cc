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

#include <honei/lbm/update_velocity_directions_grid.hh>
#include <honei/backends/opencl/operations.hh>
#include <honei/util/profiler.hh>
#include <typeinfo>

namespace honei
{
    namespace opencl
    {
        template <typename Tag_, typename DT_>
        void common_update_velocity_directions_grid(PackedGridInfo<lbm_lattice_types::D2Q9> & info,
             PackedGridData<lbm_lattice_types::D2Q9, DT_> & data)
        {
            info.limits->lock(lm_read_only);

            void * cuda_types_cl(info.cuda_types->lock(lm_read_only, Tag_::memory_value));

            void * f_temp_1_cl(data.f_temp_1->lock(lm_read_and_write, Tag_::memory_value));
            void * f_temp_2_cl(data.f_temp_2->lock(lm_read_and_write, Tag_::memory_value));
            void * f_temp_3_cl(data.f_temp_3->lock(lm_read_and_write, Tag_::memory_value));
            void * f_temp_4_cl(data.f_temp_4->lock(lm_read_and_write, Tag_::memory_value));
            void * f_temp_5_cl(data.f_temp_5->lock(lm_read_and_write, Tag_::memory_value));
            void * f_temp_6_cl(data.f_temp_6->lock(lm_read_and_write, Tag_::memory_value));
            void * f_temp_7_cl(data.f_temp_7->lock(lm_read_and_write, Tag_::memory_value));
            void * f_temp_8_cl(data.f_temp_8->lock(lm_read_and_write, Tag_::memory_value));

            unsigned long start((*info.limits)[0]);
            unsigned long end((*info.limits)[info.limits->size() - 1]);

            std::string opname("update_velocity_directions_grid_");
            opname += typeid(DT_).name();
            opencl::update_velocity_directions_grid(start, end, cuda_types_cl,
                    f_temp_1_cl, f_temp_2_cl,
                    f_temp_3_cl, f_temp_4_cl, f_temp_5_cl,
                    f_temp_6_cl, f_temp_7_cl, f_temp_8_cl,
                    tag_to_device<Tag_>(), opname);

            info.limits->unlock(lm_read_only);

            info.cuda_types->unlock(lm_read_only);

            data.f_temp_1->unlock(lm_read_and_write);
            data.f_temp_2->unlock(lm_read_and_write);
            data.f_temp_3->unlock(lm_read_and_write);
            data.f_temp_4->unlock(lm_read_and_write);
            data.f_temp_5->unlock(lm_read_and_write);
            data.f_temp_6->unlock(lm_read_and_write);
            data.f_temp_7->unlock(lm_read_and_write);
            data.f_temp_8->unlock(lm_read_and_write);
        }
    }
}

using namespace honei;

template <typename DT_>
void UpdateVelocityDirectionsGrid<tags::OpenCL::CPU, lbm_boundary_types::NOSLIP>::value(
        PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data)
{
    CONTEXT("When updating velocity directions (OpenCL CPU):");

    opencl::common_update_velocity_directions_grid<tags::OpenCL::CPU>(info, data);
}

template void UpdateVelocityDirectionsGrid<tags::OpenCL::CPU, lbm_boundary_types::NOSLIP>::value(
        PackedGridInfo<lbm_lattice_types::D2Q9> &, PackedGridData<lbm_lattice_types::D2Q9, float> &);
template void UpdateVelocityDirectionsGrid<tags::OpenCL::CPU, lbm_boundary_types::NOSLIP>::value(
        PackedGridInfo<lbm_lattice_types::D2Q9> &, PackedGridData<lbm_lattice_types::D2Q9, double> &);

template <typename DT_>
void UpdateVelocityDirectionsGrid<tags::OpenCL::GPU, lbm_boundary_types::NOSLIP>::value(
        PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data)
{
    CONTEXT("When updating velocity directions (OpenCL GPU):");

    opencl::common_update_velocity_directions_grid<tags::OpenCL::GPU>(info, data);
}

template void UpdateVelocityDirectionsGrid<tags::OpenCL::GPU, lbm_boundary_types::NOSLIP>::value(
        PackedGridInfo<lbm_lattice_types::D2Q9> &, PackedGridData<lbm_lattice_types::D2Q9, float> &);
template void UpdateVelocityDirectionsGrid<tags::OpenCL::GPU, lbm_boundary_types::NOSLIP>::value(
        PackedGridInfo<lbm_lattice_types::D2Q9> &, PackedGridData<lbm_lattice_types::D2Q9, double> &);
