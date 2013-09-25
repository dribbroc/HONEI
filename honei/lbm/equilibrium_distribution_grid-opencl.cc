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

#include <honei/lbm/equilibrium_distribution_grid.hh>
#include <honei/backends/opencl/operations.hh>
#include <honei/util/profiler.hh>
#include <typeinfo>

namespace honei
{
    namespace opencl
    {
        template <typename Tag_, typename DT_>
        void common_eq_dist_grid(PackedGridInfo<lbm_lattice_types::D2Q9> & info,
             PackedGridData<lbm_lattice_types::D2Q9, DT_> & data, DT_ g, DT_ e)
        {
            info.limits->lock(lm_read_only);

            void * f_eq_0_cl(data.f_eq_0->lock(lm_write_only, Tag_::memory_value));
            void * f_eq_1_cl(data.f_eq_1->lock(lm_write_only, Tag_::memory_value));
            void * f_eq_2_cl(data.f_eq_2->lock(lm_write_only, Tag_::memory_value));
            void * f_eq_3_cl(data.f_eq_3->lock(lm_write_only, Tag_::memory_value));
            void * f_eq_4_cl(data.f_eq_4->lock(lm_write_only, Tag_::memory_value));
            void * f_eq_5_cl(data.f_eq_5->lock(lm_write_only, Tag_::memory_value));
            void * f_eq_6_cl(data.f_eq_6->lock(lm_write_only, Tag_::memory_value));
            void * f_eq_7_cl(data.f_eq_7->lock(lm_write_only, Tag_::memory_value));
            void * f_eq_8_cl(data.f_eq_8->lock(lm_write_only, Tag_::memory_value));

            void * h_cl(data.h->lock(lm_read_only, Tag_::memory_value));
            void * u_cl(data.u->lock(lm_read_only, Tag_::memory_value));
            void * v_cl(data.v->lock(lm_read_only, Tag_::memory_value));

            void * distribution_x_cl(data.distribution_x->lock(lm_read_only, Tag_::memory_value));
            void * distribution_y_cl(data.distribution_y->lock(lm_read_only, Tag_::memory_value));

            unsigned long start((*info.limits)[0]);
            unsigned long end((*info.limits)[info.limits->size() - 1]);

            std::string opname("eq_dist_grid_");
            opname += typeid(DT_).name();
            opencl::eq_dist_grid(start, end,
                    f_eq_0_cl, f_eq_1_cl, f_eq_2_cl,
                    f_eq_3_cl, f_eq_4_cl, f_eq_5_cl,
                    f_eq_6_cl, f_eq_7_cl, f_eq_8_cl,
                    h_cl, u_cl, v_cl,
                    distribution_x_cl, distribution_y_cl, g, e,
                    tag_to_device<Tag_>(), opname);

            info.limits->unlock(lm_read_only);

            data.f_eq_0->unlock(lm_write_only);
            data.f_eq_1->unlock(lm_write_only);
            data.f_eq_2->unlock(lm_write_only);
            data.f_eq_3->unlock(lm_write_only);
            data.f_eq_4->unlock(lm_write_only);
            data.f_eq_5->unlock(lm_write_only);
            data.f_eq_6->unlock(lm_write_only);
            data.f_eq_7->unlock(lm_write_only);
            data.f_eq_8->unlock(lm_write_only);

            data.distribution_x->unlock(lm_read_only);
            data.distribution_y->unlock(lm_read_only);

            data.h->unlock(lm_read_only);
            data.u->unlock(lm_read_only);
            data.v->unlock(lm_read_only);
        }
    }
}

using namespace honei;

template <typename Tag_>
template <typename DT_>
void EquilibriumDistributionGrid<Tag_, lbm_applications::LABSWE>::value(
        DT_ g, DT_ e,
        PackedGridInfo<lbm_lattice_types::D2Q9> & info,
        PackedGridData<lbm_lattice_types::D2Q9, DT_> & data)
{
    CONTEXT("When computing LABSWE local equilibrium distribution function (OpenCL):");

    opencl::common_eq_dist_grid<Tag_>(info, data, g, e);
}


template void EquilibriumDistributionGrid<tags::OpenCL::CPU, lbm_applications::LABSWE>::value(float, float,
        PackedGridInfo<lbm_lattice_types::D2Q9> &,
        PackedGridData<lbm_lattice_types::D2Q9, float> &);
template void EquilibriumDistributionGrid<tags::OpenCL::CPU, lbm_applications::LABSWE>::value(double, double,
        PackedGridInfo<lbm_lattice_types::D2Q9> &,
        PackedGridData<lbm_lattice_types::D2Q9, double> &);

#ifdef HONEI_OPENCL_GPU
template void EquilibriumDistributionGrid<tags::OpenCL::GPU, lbm_applications::LABSWE>::value(float, float,
        PackedGridInfo<lbm_lattice_types::D2Q9> &,
        PackedGridData<lbm_lattice_types::D2Q9, float> &);
template void EquilibriumDistributionGrid<tags::OpenCL::GPU, lbm_applications::LABSWE>::value(double, double,
        PackedGridInfo<lbm_lattice_types::D2Q9> &,
        PackedGridData<lbm_lattice_types::D2Q9, double> &);
#endif

#ifdef HONEI_OPENCL_ACC
template void EquilibriumDistributionGrid<tags::OpenCL::Accelerator, lbm_applications::LABSWE>::value(float, float,
        PackedGridInfo<lbm_lattice_types::D2Q9> &,
        PackedGridData<lbm_lattice_types::D2Q9, float> &);
template void EquilibriumDistributionGrid<tags::OpenCL::Accelerator, lbm_applications::LABSWE>::value(double, double,
        PackedGridInfo<lbm_lattice_types::D2Q9> &,
        PackedGridData<lbm_lattice_types::D2Q9, double> &);
#endif
