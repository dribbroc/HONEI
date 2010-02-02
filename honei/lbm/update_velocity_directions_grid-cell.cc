/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2009 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/lbm/update_velocity_directions_grid.hh>



using namespace honei;
using namespace cell;

void UpdateVelocityDirectionsGrid<tags::Cell, lbm_boundary_types::NOSLIP>::value(
        PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data, float tau)
{
    CONTEXT("When updating velocity directions (Cell):");
#error "Last patch in this module has not yet been applied. Boundary treatment is wrong!"

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

    unsigned long * types(info.types->elements());
    unsigned long * limits(info.limits->elements());
    float * f_temp_1(data.f_temp_1->elements());
    float * f_temp_2(data.f_temp_2->elements());
    float * f_temp_3(data.f_temp_3->elements());
    float * f_temp_4(data.f_temp_4->elements());
    float * f_temp_5(data.f_temp_5->elements());
    float * f_temp_6(data.f_temp_6->elements());
    float * f_temp_7(data.f_temp_7->elements());
    float * f_temp_8(data.f_temp_8->elements());

    for (unsigned long begin(0) ; begin < info.limits->size() - 1 ; ++begin)
    {
        if((types[begin] & 1<<0) == 1<<0)
            for (unsigned long i(limits[begin]) ; i != limits[begin + 1] ; ++i)
            {
                f_temp_5[i] = f_temp_1[i];
            }
        if((types[begin] & 1<<1) == 1<<1)
            for (unsigned long i(limits[begin]) ; i != limits[begin + 1] ; ++i)
            {
                f_temp_6[i] = f_temp_2[i];
            }
        if((types[begin] & 1<<2) == 1<<2)
            for (unsigned long i(limits[begin]) ; i != limits[begin + 1] ; ++i)
            {
                f_temp_7[i] = f_temp_3[i];
            }
        if((types[begin] & 1<<3) == 1<<3)
            for (unsigned long i(limits[begin]) ; i != limits[begin + 1] ; ++i)
            {
                f_temp_8[i] = f_temp_4[i];
            }
        if((types[begin] & 1<<4) == 1<<4)
            for (unsigned long i(limits[begin]) ; i != limits[begin + 1] ; ++i)
            {
                f_temp_1[i] = f_temp_5[i];
            }
        if((types[begin] & 1<<5) == 1<<5)
            for (unsigned long i(limits[begin]) ; i != limits[begin + 1] ; ++i)
            {
                f_temp_2[i] = f_temp_6[i];
            }
        if((types[begin] & 1<<6) == 1<<6)
            for (unsigned long i(limits[begin]) ; i != limits[begin + 1] ; ++i)
            {
                f_temp_3[i] = f_temp_7[i];
            }
        if((types[begin] & 1<<7) == 1<<7)
            for (unsigned long i(limits[begin]) ; i != limits[begin + 1] ; ++i)
            {
                f_temp_4[i] = f_temp_8[i];
            }

        // Corners
        if((types[begin] & 1<<2) == 1<<2 && (types[begin] & 1<<4) == 1<<4)
            for (unsigned long i(limits[begin]) ; i != limits[begin + 1] ; ++i)
            {
                f_temp_2[i] = f_temp_8[i];
                f_temp_6[i] = f_temp_8[i];
            }
        if((types[begin] & 1<<4) == 1<<4 && (types[begin] & 1<<6) == 1<<6)
            for (unsigned long i(limits[begin]) ; i != limits[begin + 1] ; ++i)
            {
                f_temp_4[i] = f_temp_2[i];
                f_temp_8[i] = f_temp_2[i];
            }
        if((types[begin] & 1<<0) == 1<<0 && (types[begin] & 1<<6) == 1<<6)
            for (unsigned long i(limits[begin]) ; i != limits[begin + 1] ; ++i)
            {
                f_temp_2[i] = f_temp_4[i];
                f_temp_6[i] = f_temp_4[i];
            }
        if((types[begin] & 1<<0) == 1<<0 && (types[begin] & 1<<2) == 1<<2)
            for (unsigned long i(limits[begin]) ; i != limits[begin + 1] ; ++i)
            {
                f_temp_4[i] = f_temp_6[i];
                f_temp_8[i] = f_temp_6[i];
            }
    }

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
