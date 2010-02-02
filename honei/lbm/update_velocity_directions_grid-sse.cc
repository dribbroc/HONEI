/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2009 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
#include <honei/backends/sse/operations.hh>

using namespace honei;

template <typename DT_>
void UpdateVelocityDirectionsGrid<tags::CPU::SSE, lbm_boundary_types::NOSLIP>::value(
        PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data, DT_ tau)
{
    CONTEXT("When updating velocity directions (SSE):");

    info.limits->lock(lm_read_only);
    info.types->lock(lm_read_only);

    data.f_temp_1->lock(lm_read_and_write);
    data.f_temp_2->lock(lm_read_and_write);
    data.f_temp_3->lock(lm_read_and_write);
    data.f_temp_4->lock(lm_read_and_write);
    data.f_temp_5->lock(lm_read_and_write);
    data.f_temp_6->lock(lm_read_and_write);
    data.f_temp_7->lock(lm_read_and_write);
    data.f_temp_8->lock(lm_read_and_write);

    unsigned long begin(0);
    unsigned long end(info.types->size() - 1);
    sse::up_vel_dir_grid(begin, end, info.limits->elements(), info.types->elements(),
            data.f_temp_1->elements(), data.f_temp_2->elements(), data.f_temp_3->elements(),
            data.f_temp_4->elements(), data.f_temp_5->elements(), data.f_temp_6->elements(),
            data.f_temp_7->elements(), data.f_temp_8->elements(), tau);

    info.limits->unlock(lm_read_only);
    info.types->unlock(lm_read_only);

    data.f_temp_1->unlock(lm_read_and_write);
    data.f_temp_2->unlock(lm_read_and_write);
    data.f_temp_3->unlock(lm_read_and_write);
    data.f_temp_4->unlock(lm_read_and_write);
    data.f_temp_5->unlock(lm_read_and_write);
    data.f_temp_6->unlock(lm_read_and_write);
    data.f_temp_7->unlock(lm_read_and_write);
    data.f_temp_8->unlock(lm_read_and_write);
}
template void UpdateVelocityDirectionsGrid<tags::CPU::SSE, lbm_boundary_types::NOSLIP>::value<float>(
        PackedGridInfo<D2Q9> &, PackedGridData<D2Q9, float> &, float);

template void UpdateVelocityDirectionsGrid<tags::CPU::SSE, lbm_boundary_types::NOSLIP>::value<double>(
        PackedGridInfo<D2Q9> &, PackedGridData<D2Q9, double> &, double);
