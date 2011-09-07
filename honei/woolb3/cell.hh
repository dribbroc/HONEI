/* vim: set sw=4 sts=4 et foldmethod=syntax : */

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
#ifndef WOOLB3_GUARD_CELL_HH
#define WOOLB3_GUARD_CELL_HH 1

#include <honei/util/assertion.hh>
#include <honei/la/dense_vector.hh>

#include <list>

namespace honei
{
    template <typename DT_, unsigned long directions>
        class Cell
        {
            private:
                unsigned long _id;
                //center coordinates and half of height/width
                double _x, _y, _d_x, _d_y;
                //cell values
                DT_ _h, _u, _v, _b;
                // list of neighbours, ccw in each direction, empty list means boundary
                // for different boundary types one may use dummy cells in the apropriate direction.
                std::list<Cell *> _neighbours[directions];

            public:
                Cell(double x, double y, double d_x, double d_y, DT_ h, DT_ b, DT_ u, DT_ v) :
                    _id(4711),
                    _x(x),
                    _y(y),
                    _d_x(d_x / DT_(2)),
                    _d_y(d_y / DT_(2)),
                    _h(h),
                    _u(u),
                    _v(v),
                    _b(b)
                {
                }

                void set_id(unsigned long id)
                {
                    _id = id;
                }

                unsigned long get_id()
                {
                    return _id;
                }

                double get_x()
                {
                    return _x;
                }

                double get_y()
                {
                    return _y;
                }

                double get_d_x()
                {
                    return _d_x;
                }

                double get_d_y()
                {
                    return _d_y;
                }

                DT_ get_h()
                {
                    return _h;
                }

                DT_ get_b()
                {
                    return _b;
                }

                DT_ get_u()
                {
                    return _u;
                }

                DT_ get_v()
                {
                    return _v;
                }

                void add_neighbour(Cell* neighbour, unsigned long direction)
                {
                    ASSERT(direction < directions, "Cell::add_neighbour called with illegal direction");
                    _neighbours[direction].push_back(neighbour);
                }

                std::list<Cell*> get_neighbours(unsigned long direction)
                {
                    ASSERT(direction < directions, "Cell::get_neighbours called with illegal direction");
                    return _neighbours[direction];
                }
        };
}

#endif
