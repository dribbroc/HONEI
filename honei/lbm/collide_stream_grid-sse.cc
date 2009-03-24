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

#include <honei/lbm/collide_stream_grid.hh>
#include <honei/backends/sse/operations.hh>

using namespace honei;

void CollideStreamGrid<tags::CPU::SSE, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>::value(
        PackedGridInfo<lbm_lattice_types::D2Q9> & info,
        PackedGridData<lbm_lattice_types::D2Q9, float> & data,
        float tau)
{
    CONTEXT("When performing collision and streaming (SSE):");

    info.limits->lock(lm_read_only);
    info.dir_1->lock(lm_read_only);
    info.dir_2->lock(lm_read_only);
    info.dir_3->lock(lm_read_only);
    info.dir_4->lock(lm_read_only);
    info.dir_5->lock(lm_read_only);
    info.dir_6->lock(lm_read_only);
    info.dir_7->lock(lm_read_only);
    info.dir_8->lock(lm_read_only);
    info.dir_index_1->lock(lm_read_only);
    info.dir_index_2->lock(lm_read_only);
    info.dir_index_3->lock(lm_read_only);
    info.dir_index_4->lock(lm_read_only);
    info.dir_index_5->lock(lm_read_only);
    info.dir_index_6->lock(lm_read_only);
    info.dir_index_7->lock(lm_read_only);
    info.dir_index_8->lock(lm_read_only);

    data.f_eq_0->lock(lm_read_only);
    data.f_eq_1->lock(lm_read_only);
    data.f_eq_2->lock(lm_read_only);
    data.f_eq_3->lock(lm_read_only);
    data.f_eq_4->lock(lm_read_only);
    data.f_eq_5->lock(lm_read_only);
    data.f_eq_6->lock(lm_read_only);
    data.f_eq_7->lock(lm_read_only);
    data.f_eq_8->lock(lm_read_only);
    data.f_0->lock(lm_read_only);
    data.f_1->lock(lm_read_only);
    data.f_2->lock(lm_read_only);
    data.f_3->lock(lm_read_only);
    data.f_4->lock(lm_read_only);
    data.f_5->lock(lm_read_only);
    data.f_6->lock(lm_read_only);
    data.f_7->lock(lm_read_only);
    data.f_8->lock(lm_read_only);

    data.f_temp_0->lock(lm_write_only);
    data.f_temp_1->lock(lm_write_only);
    data.f_temp_2->lock(lm_write_only);
    data.f_temp_3->lock(lm_write_only);
    data.f_temp_4->lock(lm_write_only);
    data.f_temp_5->lock(lm_write_only);
    data.f_temp_6->lock(lm_write_only);
    data.f_temp_7->lock(lm_write_only);
    data.f_temp_8->lock(lm_write_only);


    sse::collide_stream_grid_dir_0((*info.limits)[0], (*info.limits)[info.limits->size() - 1], tau,
            data.f_temp_0->elements(), data.f_0->elements(), data.f_eq_0->elements());

    sse::collide_stream_grid_dir_n(info.dir_index_1->size() - 1, tau,
            info.dir_1->elements(), info.dir_index_1->elements(),
            data.f_temp_1->elements(), data.f_1->elements(), data.f_eq_1->elements());

    sse::collide_stream_grid_dir_n(info.dir_index_2->size() - 1, tau,
            info.dir_2->elements(), info.dir_index_2->elements(),
            data.f_temp_2->elements(), data.f_2->elements(), data.f_eq_2->elements());

    sse::collide_stream_grid_dir_n(info.dir_index_3->size() - 1, tau,
            info.dir_3->elements(), info.dir_index_3->elements(),
            data.f_temp_3->elements(), data.f_3->elements(), data.f_eq_3->elements());

    sse::collide_stream_grid_dir_n(info.dir_index_4->size() - 1, tau,
            info.dir_4->elements(), info.dir_index_4->elements(),
            data.f_temp_4->elements(), data.f_4->elements(), data.f_eq_4->elements());

    sse::collide_stream_grid_dir_n(info.dir_index_5->size() - 1, tau,
            info.dir_5->elements(), info.dir_index_5->elements(),
            data.f_temp_5->elements(), data.f_5->elements(), data.f_eq_5->elements());

    sse::collide_stream_grid_dir_n(info.dir_index_6->size() - 1, tau,
            info.dir_6->elements(), info.dir_index_6->elements(),
            data.f_temp_6->elements(), data.f_6->elements(), data.f_eq_6->elements());

    sse::collide_stream_grid_dir_n(info.dir_index_7->size() - 1, tau,
            info.dir_7->elements(), info.dir_index_7->elements(),
            data.f_temp_7->elements(), data.f_7->elements(), data.f_eq_7->elements());

    sse::collide_stream_grid_dir_n(info.dir_index_8->size() - 1, tau,
            info.dir_8->elements(), info.dir_index_8->elements(),
            data.f_temp_8->elements(), data.f_8->elements(), data.f_eq_8->elements());

    info.limits->unlock(lm_read_only);
    info.dir_1->unlock(lm_read_only);
    info.dir_2->unlock(lm_read_only);
    info.dir_3->unlock(lm_read_only);
    info.dir_4->unlock(lm_read_only);
    info.dir_5->unlock(lm_read_only);
    info.dir_6->unlock(lm_read_only);
    info.dir_7->unlock(lm_read_only);
    info.dir_8->unlock(lm_read_only);
    info.dir_index_1->unlock(lm_read_only);
    info.dir_index_2->unlock(lm_read_only);
    info.dir_index_3->unlock(lm_read_only);
    info.dir_index_4->unlock(lm_read_only);
    info.dir_index_5->unlock(lm_read_only);
    info.dir_index_6->unlock(lm_read_only);
    info.dir_index_7->unlock(lm_read_only);
    info.dir_index_8->unlock(lm_read_only);

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

void CollideStreamGrid<tags::CPU::SSE, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>::value(
        PackedGridInfo<lbm_lattice_types::D2Q9> & info,
        PackedGridData<lbm_lattice_types::D2Q9, double> & data,
        double tau)
{
    CONTEXT("When performing collision and streaming (SSE):");

    info.limits->lock(lm_read_only);
    info.dir_1->lock(lm_read_only);
    info.dir_2->lock(lm_read_only);
    info.dir_3->lock(lm_read_only);
    info.dir_4->lock(lm_read_only);
    info.dir_5->lock(lm_read_only);
    info.dir_6->lock(lm_read_only);
    info.dir_7->lock(lm_read_only);
    info.dir_8->lock(lm_read_only);
    info.dir_index_1->lock(lm_read_only);
    info.dir_index_2->lock(lm_read_only);
    info.dir_index_3->lock(lm_read_only);
    info.dir_index_4->lock(lm_read_only);
    info.dir_index_5->lock(lm_read_only);
    info.dir_index_6->lock(lm_read_only);
    info.dir_index_7->lock(lm_read_only);
    info.dir_index_8->lock(lm_read_only);

    data.f_eq_0->lock(lm_read_only);
    data.f_eq_1->lock(lm_read_only);
    data.f_eq_2->lock(lm_read_only);
    data.f_eq_3->lock(lm_read_only);
    data.f_eq_4->lock(lm_read_only);
    data.f_eq_5->lock(lm_read_only);
    data.f_eq_6->lock(lm_read_only);
    data.f_eq_7->lock(lm_read_only);
    data.f_eq_8->lock(lm_read_only);
    data.f_0->lock(lm_read_only);
    data.f_1->lock(lm_read_only);
    data.f_2->lock(lm_read_only);
    data.f_3->lock(lm_read_only);
    data.f_4->lock(lm_read_only);
    data.f_5->lock(lm_read_only);
    data.f_6->lock(lm_read_only);
    data.f_7->lock(lm_read_only);
    data.f_8->lock(lm_read_only);

    data.f_temp_0->lock(lm_write_only);
    data.f_temp_1->lock(lm_write_only);
    data.f_temp_2->lock(lm_write_only);
    data.f_temp_3->lock(lm_write_only);
    data.f_temp_4->lock(lm_write_only);
    data.f_temp_5->lock(lm_write_only);
    data.f_temp_6->lock(lm_write_only);
    data.f_temp_7->lock(lm_write_only);
    data.f_temp_8->lock(lm_write_only);

    sse::collide_stream_grid_dir_0((*info.limits)[0], (*info.limits)[info.limits->size() - 1], tau,
            data.f_temp_0->elements(), data.f_0->elements(), data.f_eq_0->elements());

    sse::collide_stream_grid_dir_n(info.dir_index_1->size() - 1, tau,
            info.dir_1->elements(), info.dir_index_1->elements(),
            data.f_temp_1->elements(), data.f_1->elements(), data.f_eq_1->elements());

    sse::collide_stream_grid_dir_n(info.dir_index_2->size() - 1, tau,
            info.dir_2->elements(), info.dir_index_2->elements(),
            data.f_temp_2->elements(), data.f_2->elements(), data.f_eq_2->elements());

    sse::collide_stream_grid_dir_n(info.dir_index_3->size() - 1, tau,
            info.dir_3->elements(), info.dir_index_3->elements(),
            data.f_temp_3->elements(), data.f_3->elements(), data.f_eq_3->elements());

    sse::collide_stream_grid_dir_n(info.dir_index_4->size() - 1, tau,
            info.dir_4->elements(), info.dir_index_4->elements(),
            data.f_temp_4->elements(), data.f_4->elements(), data.f_eq_4->elements());

    sse::collide_stream_grid_dir_n(info.dir_index_5->size() - 1, tau,
            info.dir_5->elements(), info.dir_index_5->elements(),
            data.f_temp_5->elements(), data.f_5->elements(), data.f_eq_5->elements());

    sse::collide_stream_grid_dir_n(info.dir_index_6->size() - 1, tau,
            info.dir_6->elements(), info.dir_index_6->elements(),
            data.f_temp_6->elements(), data.f_6->elements(), data.f_eq_6->elements());

    sse::collide_stream_grid_dir_n(info.dir_index_7->size() - 1, tau,
            info.dir_7->elements(), info.dir_index_7->elements(),
            data.f_temp_7->elements(), data.f_7->elements(), data.f_eq_7->elements());

    sse::collide_stream_grid_dir_n(info.dir_index_8->size() - 1, tau,
            info.dir_8->elements(), info.dir_index_8->elements(),
            data.f_temp_8->elements(), data.f_8->elements(), data.f_eq_8->elements());

    info.limits->unlock(lm_read_only);
    info.dir_1->unlock(lm_read_only);
    info.dir_2->unlock(lm_read_only);
    info.dir_3->unlock(lm_read_only);
    info.dir_4->unlock(lm_read_only);
    info.dir_5->unlock(lm_read_only);
    info.dir_6->unlock(lm_read_only);
    info.dir_7->unlock(lm_read_only);
    info.dir_8->unlock(lm_read_only);
    info.dir_index_1->unlock(lm_read_only);
    info.dir_index_2->unlock(lm_read_only);
    info.dir_index_3->unlock(lm_read_only);
    info.dir_index_4->unlock(lm_read_only);
    info.dir_index_5->unlock(lm_read_only);
    info.dir_index_6->unlock(lm_read_only);
    info.dir_index_7->unlock(lm_read_only);
    info.dir_index_8->unlock(lm_read_only);

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
