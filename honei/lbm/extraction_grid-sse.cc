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

#include <honei/lbm/extraction_grid.hh>
#include <honei/backends/sse/operations.hh>
#include <honei/util/type_traits.hh>

using namespace honei;

void ExtractionGrid<tags::CPU::SSE, lbm_modes::DRY>::value(
        PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data, float epsilon)
{
    CONTEXT("When extracting h, u and v(SSE):");

    //set f to t_temp
    DenseVector<float> * swap;
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

    info.limits->lock(lm_read_only);

    data.f_0->lock(lm_read_and_write);
    data.f_1->lock(lm_read_and_write);
    data.f_2->lock(lm_read_and_write);
    data.f_3->lock(lm_read_and_write);
    data.f_4->lock(lm_read_and_write);
    data.f_5->lock(lm_read_and_write);
    data.f_6->lock(lm_read_and_write);
    data.f_7->lock(lm_read_and_write);
    data.f_8->lock(lm_read_and_write);

    data.h->lock(lm_write_only);

    data.distribution_x->lock(lm_read_only);
    data.distribution_y->lock(lm_read_only);

    data.u->lock(lm_write_only);
    data.v->lock(lm_write_only);

    unsigned long begin((*info.limits)[0]);
    unsigned long end((*info.limits)[info.limits->size() - 1]);

    sse::extraction_grid_dry(begin, end,
            data.distribution_x->elements(), data.distribution_y->elements(),
            data.h->elements(), data.u->elements(), data.v->elements(),
            data.f_0->elements(), data.f_1->elements(), data.f_2->elements(),
            data.f_3->elements(), data.f_4->elements(), data.f_5->elements(),
            data.f_6->elements(), data.f_7->elements(), data.f_8->elements(), epsilon);

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

void ExtractionGrid<tags::CPU::SSE, lbm_modes::DRY>::value(
        PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, double> & data, double epsilon)
{
    CONTEXT("When extracting h, u and v(SSE):");

    //set f to t_temp
    DenseVector<double> * swap;
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

    info.limits->lock(lm_read_only);

    data.f_0->lock(lm_read_and_write);
    data.f_1->lock(lm_read_and_write);
    data.f_2->lock(lm_read_and_write);
    data.f_3->lock(lm_read_and_write);
    data.f_4->lock(lm_read_and_write);
    data.f_5->lock(lm_read_and_write);
    data.f_6->lock(lm_read_and_write);
    data.f_7->lock(lm_read_and_write);
    data.f_8->lock(lm_read_and_write);

    data.h->lock(lm_write_only);

    data.distribution_x->lock(lm_read_only);
    data.distribution_y->lock(lm_read_only);

    data.u->lock(lm_write_only);
    data.v->lock(lm_write_only);

    unsigned long begin((*info.limits)[0]);
    unsigned long end((*info.limits)[info.limits->size() - 1]);

    sse::extraction_grid_dry(begin, end,
            data.distribution_x->elements(), data.distribution_y->elements(),
            data.h->elements(), data.u->elements(), data.v->elements(),
            data.f_0->elements(), data.f_1->elements(), data.f_2->elements(),
            data.f_3->elements(), data.f_4->elements(), data.f_5->elements(),
            data.f_6->elements(), data.f_7->elements(), data.f_8->elements(), epsilon);

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

void ExtractionGrid<tags::CPU::SSE, lbm_modes::WET>::value(
        PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data, float epsilon)
{
    CONTEXT("When extracting h, u and v(SSE):");

    //set f to t_temp
    DenseVector<float> * swap;
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

    info.limits->lock(lm_read_only);

    data.f_0->lock(lm_read_and_write);
    data.f_1->lock(lm_read_and_write);
    data.f_2->lock(lm_read_and_write);
    data.f_3->lock(lm_read_and_write);
    data.f_4->lock(lm_read_and_write);
    data.f_5->lock(lm_read_and_write);
    data.f_6->lock(lm_read_and_write);
    data.f_7->lock(lm_read_and_write);
    data.f_8->lock(lm_read_and_write);

    data.h->lock(lm_write_only);

    data.distribution_x->lock(lm_read_only);
    data.distribution_y->lock(lm_read_only);

    data.u->lock(lm_write_only);
    data.v->lock(lm_write_only);

    unsigned long begin((*info.limits)[0]);
    unsigned long end((*info.limits)[info.limits->size() - 1]);

    sse::extraction_grid_wet(begin, end,
            data.distribution_x->elements(), data.distribution_y->elements(),
            data.h->elements(), data.u->elements(), data.v->elements(),
            data.f_0->elements(), data.f_1->elements(), data.f_2->elements(),
            data.f_3->elements(), data.f_4->elements(), data.f_5->elements(),
            data.f_6->elements(), data.f_7->elements(), data.f_8->elements(), epsilon);

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

void ExtractionGrid<tags::CPU::SSE, lbm_modes::WET>::value(
        PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, double> & data, double epsilon)
{
    CONTEXT("When extracting h, u and v(SSE):");

    //set f to t_temp
    DenseVector<double> * swap;
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

    info.limits->lock(lm_read_only);

    data.f_0->lock(lm_read_and_write);
    data.f_1->lock(lm_read_and_write);
    data.f_2->lock(lm_read_and_write);
    data.f_3->lock(lm_read_and_write);
    data.f_4->lock(lm_read_and_write);
    data.f_5->lock(lm_read_and_write);
    data.f_6->lock(lm_read_and_write);
    data.f_7->lock(lm_read_and_write);
    data.f_8->lock(lm_read_and_write);

    data.h->lock(lm_write_only);

    data.distribution_x->lock(lm_read_only);
    data.distribution_y->lock(lm_read_only);

    data.u->lock(lm_write_only);
    data.v->lock(lm_write_only);

    unsigned long begin((*info.limits)[0]);
    unsigned long end((*info.limits)[info.limits->size() - 1]);

    sse::extraction_grid_wet(begin, end,
            data.distribution_x->elements(), data.distribution_y->elements(),
            data.h->elements(), data.u->elements(), data.v->elements(),
            data.f_0->elements(), data.f_1->elements(), data.f_2->elements(),
            data.f_3->elements(), data.f_4->elements(), data.f_5->elements(),
            data.f_6->elements(), data.f_7->elements(), data.f_8->elements(), epsilon);

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
