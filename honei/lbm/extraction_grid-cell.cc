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

#include <honei/la/algorithm.hh>
#include <honei/la/sum.hh>
#include <honei/la/scaled_sum.hh>

using namespace honei;

void ExtractionGrid<tags::Cell, lbm_modes::WET>::value(
        PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data, float epsilon)
{
    CONTEXT("When extracting h, u and v(Cell):");

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

    /*info.limits->lock(lm_read_only);

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
    data.v->unlock(lm_write_only); */
    //honei::fill<tags::Cell>(*data.h, float(0));
    Sum<tags::Cell>::value(*data.h, *data.f_0);
    Sum<tags::Cell>::value(*data.h, *data.f_1);
    Sum<tags::Cell>::value(*data.h, *data.f_2);
    Sum<tags::Cell>::value(*data.h, *data.f_3);
    Sum<tags::Cell>::value(*data.h, *data.f_4);
    Sum<tags::Cell>::value(*data.h, *data.f_5);
    Sum<tags::Cell>::value(*data.h, *data.f_6);
    Sum<tags::Cell>::value(*data.h, *data.f_7);
    Sum<tags::Cell>::value(*data.h, *data.f_8);

    //honei::fill<tags::Cell>(*data.u, float(0));
    ScaledSum<tags::Cell>::value(*data.u, *data.f_0, (*data.distribution_x)[0]);
    ScaledSum<tags::Cell>::value(*data.u, *data.f_1, (*data.distribution_x)[1]);
    ScaledSum<tags::Cell>::value(*data.u, *data.f_2, (*data.distribution_x)[2]);
    ScaledSum<tags::Cell>::value(*data.u, *data.f_3, (*data.distribution_x)[3]);
    ScaledSum<tags::Cell>::value(*data.u, *data.f_4, (*data.distribution_x)[4]);
    ScaledSum<tags::Cell>::value(*data.u, *data.f_5, (*data.distribution_x)[5]);
    ScaledSum<tags::Cell>::value(*data.u, *data.f_6, (*data.distribution_x)[6]);
    ScaledSum<tags::Cell>::value(*data.u, *data.f_7, (*data.distribution_x)[7]);
    ScaledSum<tags::Cell>::value(*data.u, *data.f_8, (*data.distribution_x)[8]);

    //honei::fill<tags::Cell>(*data.v, float(0));
    ScaledSum<tags::Cell>::value(*data.v, *data.f_0, (*data.distribution_y)[0]);
    ScaledSum<tags::Cell>::value(*data.v, *data.f_1, (*data.distribution_y)[1]);
    ScaledSum<tags::Cell>::value(*data.v, *data.f_2, (*data.distribution_y)[2]);
    ScaledSum<tags::Cell>::value(*data.v, *data.f_3, (*data.distribution_y)[3]);
    ScaledSum<tags::Cell>::value(*data.v, *data.f_4, (*data.distribution_y)[4]);
    ScaledSum<tags::Cell>::value(*data.v, *data.f_5, (*data.distribution_y)[5]);
    ScaledSum<tags::Cell>::value(*data.v, *data.f_6, (*data.distribution_y)[6]);
    ScaledSum<tags::Cell>::value(*data.v, *data.f_7, (*data.distribution_y)[7]);
    ScaledSum<tags::Cell>::value(*data.v, *data.f_8, (*data.distribution_y)[8]);
}
