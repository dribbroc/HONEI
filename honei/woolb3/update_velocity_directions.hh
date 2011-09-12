/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2011 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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


#pragma once
#ifndef WOOLB3_GUARD_UPDATE_VELOCITY_DIRECTIONS_HH
#define WOOLB3_GUARD_UPDATE_VELOCITY_DIRECTIONS_HH 1



#include <honei/lbm/tags.hh>
#include <honei/la/dense_vector.hh>
#include <honei/woolb3/grid3.hh>
#include <honei/woolb3/packed_grid3.hh>
#include <cmath>

namespace honei
{
    template <typename Tag_>
    struct UpdateVelocityDirections
    {
    };

    template <>
    struct UpdateVelocityDirections<tags::CPU>
    {
        template <typename DT_, unsigned long directions>
        static void value(Grid3<DT_, directions> & grid, PackedGrid3<DT_, directions> & pgrid)
        {
            for (unsigned long i(0) ; i < pgrid.h->size() ; ++i)
            {
                for (unsigned long direction(1) ; direction < directions ; ++direction)
                if(grid.get_cell(i)->get_neighbours(direction).size() == 0)
                {
                    (*pgrid.f_temp[direction + 4 < 9 ? direction + 4 : direction - 4])[i] = (*pgrid.f_temp[direction])[i];
                }

                //corners
                if(grid.get_cell(i)->get_neighbours(3).size() == 0 && grid.get_cell(i)->get_neighbours(5).size() == 0)
                {
                    (*pgrid.f_temp[2])[i] = (*pgrid.f_temp[8])[i];
                    (*pgrid.f_temp[6])[i] = (*pgrid.f_temp[8])[i];
                }
                if(grid.get_cell(i)->get_neighbours(5).size() == 0 && grid.get_cell(i)->get_neighbours(7).size() == 0)
                {
                    (*pgrid.f_temp[4])[i] = (*pgrid.f_temp[2])[i];
                    (*pgrid.f_temp[8])[i] = (*pgrid.f_temp[2])[i];
                }
                if(grid.get_cell(i)->get_neighbours(1).size() == 0 && grid.get_cell(i)->get_neighbours(7).size() == 0)
                {
                    (*pgrid.f_temp[2])[i] = (*pgrid.f_temp[4])[i];
                    (*pgrid.f_temp[6])[i] = (*pgrid.f_temp[4])[i];
                }
                if(grid.get_cell(i)->get_neighbours(1).size() == 0 && grid.get_cell(i)->get_neighbours(3).size() == 0)
                {
                    (*pgrid.f_temp[4])[i] = (*pgrid.f_temp[6])[i];
                    (*pgrid.f_temp[8])[i] = (*pgrid.f_temp[6])[i];
                }
            }
        }
    };
}
#endif
