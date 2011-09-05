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
#ifndef WOOLB3_GUARD_GRID_HH
#define WOOLB3_GUARD_GRID_HH 1

#include <honei/woolb3/cell.hh>
#include <honei/la/dense_matrix.hh>

#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>


namespace honei
{
    template <typename DT_, unsigned long directions>
        class Grid
        {
            private:
                std::vector<Cell<DT_, directions> *> _cells;

            public:
                Grid()
                {
                }

                ~Grid()
                {
                    for (unsigned long i(0) ; i < _cells.size() ; ++i)
                    {
                        delete _cells.at(i);
                    }
                }

                unsigned long size()
                {
                    return _cells.size();
                }

                void add_cell(Cell<DT_, directions> * cell)
                {
                    cell->set_id(_cells.size());
                    _cells.push_back(cell);
                }

                Cell<DT_, directions> * get_cell(unsigned long i)
                {
                    return _cells.at(i);
                }

                Grid(DenseMatrix<bool> & geometry, DenseMatrix<DT_> & h, DenseMatrix<DT_> & b, DenseMatrix<DT_> & u, DenseMatrix<DT_> & v)
                {
                    // read in cells
                    for (unsigned long row(0) ; row < geometry.rows() ; ++row)
                    {
                        for (unsigned long col(0) ; col < geometry.columns() ; ++col)
                        {
                            Cell<DT_, directions>* cell = new Cell<DT_, directions>(row + 0.5, col + 0.5, 1, 1,
                                    h(row, col), b(row, col), u(row, col), v(row, col));
                            this->add_cell(cell);
                        }
                    }

                    // set neighbourhood
                    unsigned long idx(0);
                    for (unsigned long row(0) ; row < geometry.rows() ; ++row)
                    {
                        for (unsigned long col(0) ; col < geometry.columns() ; ++col)
                        {
                            if (geometry(row, col) == false)
                            {
                                // DIR 1
                                if(col < geometry.columns() - 1 && geometry(row, col + 1) == false)
                                    _cells.at(idx)->add_neighbour(_cells.at(idx + 1), 1);

                                // DIR 2
                                if(row > 0 && col < geometry.columns() - 1 && geometry(row - 1, col + 1) == false)
                                    _cells.at(idx)->add_neighbour(_cells.at(geometry.columns() * (row-1) + (col+1)), 2);

                                // DIR 3
                                if(row > 0 && geometry(row - 1, col) == false)
                                    _cells.at(idx)->add_neighbour(_cells.at(geometry.columns() * (row-1) + (col)), 3);

                                // DIR 4
                                if(row > 0 && col > 0 && geometry(row - 1, col - 1) == false)
                                    _cells.at(idx)->add_neighbour(_cells.at(geometry.columns() * (row-1) + (col-1)), 4);

                                // DIR 5
                                if(col > 0 && geometry(row, col - 1) == false)
                                    _cells.at(idx)->add_neighbour(_cells.at(idx - 1), 5);

                                // DIR 6
                                if(row < geometry.rows()  - 1&& col > 0 && geometry(row + 1, col - 1) == false)
                                    _cells.at(idx)->add_neighbour(_cells.at(geometry.columns() * (row+1) + (col-1)), 6);

                                // DIR 7
                                if(row < geometry.rows() - 1 && geometry(row + 1, col) == false)
                                    _cells.at(idx)->add_neighbour(_cells.at(geometry.columns() * (row+1) + (col)), 7);

                                // DIR 8
                                if(row < geometry.rows() - 1 && col < geometry.columns() - 1 && geometry(row + 1, col + 1) == false)
                                    _cells.at(idx)->add_neighbour(_cells.at(geometry.columns() * (row+1) + (col+1)), 8);
                            }
                            ++idx;
                        }
                    }

                    // remove obstacle cells
                    idx = 0;
                    typename std::vector<Cell<DT_, directions> *>::iterator i = _cells.begin();
                    for (unsigned long row(0) ; row < geometry.rows() ; ++row)
                    {
                        for (unsigned long col(0) ; col < geometry.columns() ; ++col)
                        {
                            if (geometry(row, col) == true)
                            {
                                i = _cells.erase(i);
                            }
                            else
                            {
                                (*i)->set_id(idx);
                                ++i;
                                ++idx;
                            }
                        }
                    }
                }
        };
}

#endif
