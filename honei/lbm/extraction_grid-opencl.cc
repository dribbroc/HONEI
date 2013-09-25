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

#include <honei/lbm/extraction_grid.hh>
#include <honei/backends/opencl/operations.hh>
#include <honei/util/profiler.hh>
#include <typeinfo>

namespace honei
{
    namespace opencl
    {
        template <typename Tag_, typename DT_>
        void common_extraction_grid_dry(PackedGridInfo<lbm_lattice_types::D2Q9> & info,
             PackedGridData<lbm_lattice_types::D2Q9, DT_> & data, DT_ epsilon)
        {
            info.limits->lock(lm_read_only);

            void * f_0_cl(data.f_0->lock(lm_read_and_write, Tag_::memory_value));
            void * f_1_cl(data.f_1->lock(lm_read_and_write, Tag_::memory_value));
            void * f_2_cl(data.f_2->lock(lm_read_and_write, Tag_::memory_value));
            void * f_3_cl(data.f_3->lock(lm_read_and_write, Tag_::memory_value));
            void * f_4_cl(data.f_4->lock(lm_read_and_write, Tag_::memory_value));
            void * f_5_cl(data.f_5->lock(lm_read_and_write, Tag_::memory_value));
            void * f_6_cl(data.f_6->lock(lm_read_and_write, Tag_::memory_value));
            void * f_7_cl(data.f_7->lock(lm_read_and_write, Tag_::memory_value));
            void * f_8_cl(data.f_8->lock(lm_read_and_write, Tag_::memory_value));

            void * h_cl(data.h->lock(lm_write_only, Tag_::memory_value));
            void * u_cl(data.u->lock(lm_write_only, Tag_::memory_value));
            void * v_cl(data.v->lock(lm_write_only, Tag_::memory_value));

            void * distribution_x_cl(data.distribution_x->lock(lm_read_only, Tag_::memory_value));
            void * distribution_y_cl(data.distribution_y->lock(lm_read_only, Tag_::memory_value));

            unsigned long start((*info.limits)[0]);
            unsigned long end((*info.limits)[info.limits->size() - 1]);

            std::string opname("extraction_grid_dry_");
            opname += typeid(DT_).name();
            opencl::extraction_grid_dry(start, end,
                    f_0_cl, f_1_cl, f_2_cl,
                    f_3_cl, f_4_cl, f_5_cl,
                    f_6_cl, f_7_cl, f_8_cl,
                    h_cl, u_cl, v_cl,
                    distribution_x_cl, distribution_y_cl, epsilon,
                    tag_to_device<Tag_>(), opname);

            info.limits->unlock(lm_read_only);

            data.f_0->unlock(lm_read_and_write);
            data.f_1->unlock(lm_read_and_write);
            data.f_2->unlock(lm_read_and_write);
            data.f_3->unlock(lm_read_and_write);
            data.f_4->unlock(lm_read_and_write);
            data.f_5->unlock(lm_read_and_write);
            data.f_6->unlock(lm_read_and_write);
            data.f_7->unlock(lm_read_and_write);
            data.f_8->unlock(lm_read_and_write);

            data.h->unlock(lm_write_only);

            data.distribution_x->unlock(lm_read_only);
            data.distribution_y->unlock(lm_read_only);

            data.u->unlock(lm_write_only);
            data.v->unlock(lm_write_only);
        }

        template <typename Tag_, typename DT_>
        void common_extraction_grid_wet(PackedGridInfo<lbm_lattice_types::D2Q9> & info,
             PackedGridData<lbm_lattice_types::D2Q9, DT_> & data, DT_ epsilon)
        {
            info.limits->lock(lm_read_only);

            void * f_0_cl(data.f_0->lock(lm_read_and_write, Tag_::memory_value));
            void * f_1_cl(data.f_1->lock(lm_read_and_write, Tag_::memory_value));
            void * f_2_cl(data.f_2->lock(lm_read_and_write, Tag_::memory_value));
            void * f_3_cl(data.f_3->lock(lm_read_and_write, Tag_::memory_value));
            void * f_4_cl(data.f_4->lock(lm_read_and_write, Tag_::memory_value));
            void * f_5_cl(data.f_5->lock(lm_read_and_write, Tag_::memory_value));
            void * f_6_cl(data.f_6->lock(lm_read_and_write, Tag_::memory_value));
            void * f_7_cl(data.f_7->lock(lm_read_and_write, Tag_::memory_value));
            void * f_8_cl(data.f_8->lock(lm_read_and_write, Tag_::memory_value));

            void * h_cl(data.h->lock(lm_write_only, Tag_::memory_value));
            void * u_cl(data.u->lock(lm_write_only, Tag_::memory_value));
            void * v_cl(data.v->lock(lm_write_only, Tag_::memory_value));

            void * distribution_x_cl(data.distribution_x->lock(lm_read_only, Tag_::memory_value));
            void * distribution_y_cl(data.distribution_y->lock(lm_read_only, Tag_::memory_value));

            unsigned long start((*info.limits)[0]);
            unsigned long end((*info.limits)[info.limits->size() - 1]);

            std::string opname("extraction_grid_wet_");
            opname += typeid(DT_).name();
            opencl::extraction_grid_wet(start, end,
                    f_0_cl, f_1_cl, f_2_cl,
                    f_3_cl, f_4_cl, f_5_cl,
                    f_6_cl, f_7_cl, f_8_cl,
                    h_cl, u_cl, v_cl,
                    distribution_x_cl, distribution_y_cl, epsilon,
                    tag_to_device<Tag_>(), opname);

            info.limits->unlock(lm_read_only);

            data.f_0->unlock(lm_read_and_write);
            data.f_1->unlock(lm_read_and_write);
            data.f_2->unlock(lm_read_and_write);
            data.f_3->unlock(lm_read_and_write);
            data.f_4->unlock(lm_read_and_write);
            data.f_5->unlock(lm_read_and_write);
            data.f_6->unlock(lm_read_and_write);
            data.f_7->unlock(lm_read_and_write);
            data.f_8->unlock(lm_read_and_write);

            data.h->unlock(lm_write_only);

            data.distribution_x->unlock(lm_read_only);
            data.distribution_y->unlock(lm_read_only);

            data.u->unlock(lm_write_only);
            data.v->unlock(lm_write_only);
        }
    }
}

using namespace honei;

template <typename Tag_>
template <typename DT_>
void ExtractionGrid<Tag_, lbm_modes::DRY>::value(
        PackedGridInfo<lbm_lattice_types::D2Q9> & info,
        PackedGridData<lbm_lattice_types::D2Q9, DT_> & data, DT_ epsilon)
{
    CONTEXT("When extracting h, u and v (OpenCL):");

    //set f to t_temp
    DenseVector<DT_> * swap;
    swap = data.f_0;
    data.f_0 = data.f_temp_0;
    data.f_temp_0 = swap;
    swap = data.f_1;
    data.f_1 = data.f_temp_1;
    data.f_temp_1 = swap;
    swap = data.f_2;
    data.f_2 = data.f_temp_2;
    data.f_temp_2 = swap;
    swap = data.f_3;
    data.f_3 = data.f_temp_3;
    data.f_temp_3 = swap;
    swap = data.f_4;
    data.f_4 = data.f_temp_4;
    data.f_temp_4 = swap;
    swap = data.f_5;
    data.f_5 = data.f_temp_5;
    data.f_temp_5 = swap;
    swap = data.f_6;
    data.f_6 = data.f_temp_6;
    data.f_temp_6 = swap;
    swap = data.f_7;
    data.f_7 = data.f_temp_7;
    data.f_temp_7 = swap;
    swap = data.f_8;
    data.f_8 = data.f_temp_8;
    data.f_temp_8 = swap;

    opencl::common_extraction_grid_dry<Tag_>(info, data, epsilon);
}

template void ExtractionGrid<tags::OpenCL::CPU, lbm_modes::DRY>::value(PackedGridInfo<lbm_lattice_types::D2Q9> &,
    PackedGridData<lbm_lattice_types::D2Q9, float> &, float);
template void ExtractionGrid<tags::OpenCL::CPU, lbm_modes::DRY>::value(PackedGridInfo<lbm_lattice_types::D2Q9> &,
    PackedGridData<lbm_lattice_types::D2Q9, double> &, double);

#ifdef HONEI_OPENCL_GPU
template void ExtractionGrid<tags::OpenCL::GPU, lbm_modes::DRY>::value(PackedGridInfo<lbm_lattice_types::D2Q9> &,
    PackedGridData<lbm_lattice_types::D2Q9, float> &, float);
template void ExtractionGrid<tags::OpenCL::GPU, lbm_modes::DRY>::value(PackedGridInfo<lbm_lattice_types::D2Q9> &,
    PackedGridData<lbm_lattice_types::D2Q9, double> &, double);
#endif

#ifdef HONEI_OPENCL_ACC
template void ExtractionGrid<tags::OpenCL::Accelerator, lbm_modes::DRY>::value(PackedGridInfo<lbm_lattice_types::D2Q9> &,
    PackedGridData<lbm_lattice_types::D2Q9, float> &, float);
template void ExtractionGrid<tags::OpenCL::Accelerator, lbm_modes::DRY>::value(PackedGridInfo<lbm_lattice_types::D2Q9> &,
    PackedGridData<lbm_lattice_types::D2Q9, double> &, double);
#endif

template <typename Tag_>
template <typename DT_>
void ExtractionGrid<Tag_, lbm_modes::WET>::value(
        PackedGridInfo<lbm_lattice_types::D2Q9> & info,
        PackedGridData<lbm_lattice_types::D2Q9, DT_> & data, DT_ epsilon)
{
    CONTEXT("When extracting h, u and v (OpenCL):");

    //set f to t_temp
    DenseVector<DT_> * swap;
    swap = data.f_0;
    data.f_0 = data.f_temp_0;
    data.f_temp_0 = swap;
    swap = data.f_1;
    data.f_1 = data.f_temp_1;
    data.f_temp_1 = swap;
    swap = data.f_2;
    data.f_2 = data.f_temp_2;
    data.f_temp_2 = swap;
    swap = data.f_3;
    data.f_3 = data.f_temp_3;
    data.f_temp_3 = swap;
    swap = data.f_4;
    data.f_4 = data.f_temp_4;
    data.f_temp_4 = swap;
    swap = data.f_5;
    data.f_5 = data.f_temp_5;
    data.f_temp_5 = swap;
    swap = data.f_6;
    data.f_6 = data.f_temp_6;
    data.f_temp_6 = swap;
    swap = data.f_7;
    data.f_7 = data.f_temp_7;
    data.f_temp_7 = swap;
    swap = data.f_8;
    data.f_8 = data.f_temp_8;
    data.f_temp_8 = swap;

    opencl::common_extraction_grid_wet<tags::OpenCL::CPU>(info, data, epsilon);
}

template void ExtractionGrid<tags::OpenCL::CPU, lbm_modes::WET>::value(PackedGridInfo<lbm_lattice_types::D2Q9> &,
    PackedGridData<lbm_lattice_types::D2Q9, float> &, float);
template void ExtractionGrid<tags::OpenCL::CPU, lbm_modes::WET>::value(PackedGridInfo<lbm_lattice_types::D2Q9> &,
    PackedGridData<lbm_lattice_types::D2Q9, double> &, double);

#ifdef HONEI_OPENCL_GPU
template void ExtractionGrid<tags::OpenCL::GPU, lbm_modes::WET>::value(PackedGridInfo<lbm_lattice_types::D2Q9> &,
    PackedGridData<lbm_lattice_types::D2Q9, float> &, float);
template void ExtractionGrid<tags::OpenCL::GPU, lbm_modes::WET>::value(PackedGridInfo<lbm_lattice_types::D2Q9> &,
    PackedGridData<lbm_lattice_types::D2Q9, double> &, double);
#endif

#ifdef HONEI_OPENCL_ACC
template void ExtractionGrid<tags::OpenCL::Accelerator, lbm_modes::WET>::value(PackedGridInfo<lbm_lattice_types::D2Q9> &,
    PackedGridData<lbm_lattice_types::D2Q9, float> &, float);
template void ExtractionGrid<tags::OpenCL::Accelerator, lbm_modes::WET>::value(PackedGridInfo<lbm_lattice_types::D2Q9> &,
    PackedGridData<lbm_lattice_types::D2Q9, double> &, double);
#endif
