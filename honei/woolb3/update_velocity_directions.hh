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
#include <honei/util/profiler.hh>
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
        static void value(Grid3<DT_, directions> & grid, PackedGrid3<DT_, directions> & pgrid, bool inner)
        {
            PROFILER_START("UpVelDir");
            unsigned long start(inner == true ? 0 : pgrid.grid.size() - pgrid.grid.inner_halo_size());
            unsigned long end(inner == true ? pgrid.grid.local_size() : pgrid.grid.size());

            for (typename std::set<Cell<DT_, directions> *>::iterator i(grid.no_neighbour().begin()) ; i != grid.no_neighbour().end() ; ++i)
            {
                if ((*i)->get_id() >= start && (*i)->get_id() < end)
                {
                    for (unsigned long direction(1) ; direction < directions ; ++direction)
                        if((*i)->get_neighbours(direction).size() == 0)
                        {
                            (*pgrid.f_temp[direction + 4 < 9 ? direction + 4 : direction - 4])[(*i)->get_id()] = (*pgrid.f_temp[direction])[(*i)->get_id()];
                        }

                    //corners
                    if((*i)->get_neighbours(3).size() == 0 && (*i)->get_neighbours(5).size() == 0)
                    {
                        (*pgrid.f_temp[2])[(*i)->get_id()] = (*pgrid.f_temp[8])[(*i)->get_id()];
                        (*pgrid.f_temp[6])[(*i)->get_id()] = (*pgrid.f_temp[8])[(*i)->get_id()];
                    }
                    if((*i)->get_neighbours(5).size() == 0 && (*i)->get_neighbours(7).size() == 0)
                    {
                        (*pgrid.f_temp[4])[(*i)->get_id()] = (*pgrid.f_temp[2])[(*i)->get_id()];
                        (*pgrid.f_temp[8])[(*i)->get_id()] = (*pgrid.f_temp[2])[(*i)->get_id()];
                    }
                    if((*i)->get_neighbours(1).size() == 0 && (*i)->get_neighbours(7).size() == 0)
                    {
                        (*pgrid.f_temp[2])[(*i)->get_id()] = (*pgrid.f_temp[4])[(*i)->get_id()];
                        (*pgrid.f_temp[6])[(*i)->get_id()] = (*pgrid.f_temp[4])[(*i)->get_id()];
                    }
                    if((*i)->get_neighbours(1).size() == 0 && (*i)->get_neighbours(3).size() == 0)
                    {
                        (*pgrid.f_temp[4])[(*i)->get_id()] = (*pgrid.f_temp[6])[(*i)->get_id()];
                        (*pgrid.f_temp[8])[(*i)->get_id()] = (*pgrid.f_temp[6])[(*i)->get_id()];
                    }
                }
            }
            PROFILER_STOP("UpVelDir");
        }
    };
}
#endif
