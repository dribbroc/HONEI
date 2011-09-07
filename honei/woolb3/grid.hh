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

                static unsigned long _coord2idx(unsigned long x, unsigned long y, unsigned long max_col)
                {
                    std::string method("z-curve");
                    if (method.compare("z-curve") == 0)
                    {
                        unsigned long r(0);
                        for (unsigned long i(0) ; i < sizeof(unsigned long) * 8 ; ++i)
                        {
                            r |= (x & 1ul << i) << i | (y & 1ul << i) << (i + 1);
                        }
                        return r;
                    }
                    else if (method.compare("row") == 0)
                    {
                        unsigned long r = max_col * y + x;
                        return r;
                    }
                    else
                    {
                        throw InternalError(method + " is not a valid ordering scheme!");
                        return 0;
                    }
                }

                static void _idx2coord(unsigned long & x, unsigned long & y, unsigned long r, unsigned long max_col)
                {
                    std::string method("z-curve");
                    if (method.compare("z-curve") == 0)
                    {
                        unsigned long s(r);

                        for (unsigned long i(0) ; i < sizeof(unsigned long) * 8 ; i+=2)
                        {
                            x |= (s & 1ul) << (i/2);
                            s >>= 1ul;
                            y |= (s & 1ul) << (i/2);
                            s >>= 1ul;
                        }
                    }
                    else if (method.compare("row") == 0)
                    {
                        x = r % max_col;
                        y = r / max_col;
                    }
                    else
                    {
                        throw InternalError(method + " is not a valid ordering scheme!");
                    }
                }


            public:
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

                static void print_numbering(DenseMatrix<bool> & geometry)
                {
                    DenseMatrix<long> result(geometry.rows(), geometry.columns(), -1);
                    unsigned long i(0);
                    for (unsigned long idx(0) ; idx < result.size() ; ++idx)
                    {
                        unsigned long row(0);
                        unsigned long col(0);
                        _idx2coord(col, row, idx, geometry.columns());
                        if (geometry(row, col) == false)
                        {
                            result(row, col) = i;
                            ++i;
                        }
                    }
                    std::cout<<result;
                }

                Grid(DenseMatrix<bool> & geometry, DenseMatrix<DT_> & h, DenseMatrix<DT_> & b, DenseMatrix<DT_> & u, DenseMatrix<DT_> & v)
                {
                    // read in cells
                    for (unsigned long idx(0) ; idx < geometry.size() ; ++idx)
                    {
                        unsigned long row(0);
                        unsigned long col(0);
                        _idx2coord(col, row, idx, geometry.columns());
                        Cell<DT_, directions>* cell = new Cell<DT_, directions>(row + 0.5, col + 0.5, 1, 1,
                                h(row, col), b(row, col), u(row, col), v(row, col));
                        this->add_cell(cell);
                    }

                    // set neighbourhood
                    for (unsigned long idx(0) ; idx < geometry.size() ; ++idx)
                    {
                        unsigned long row(0);
                        unsigned long col(0);
                        _idx2coord(col, row, idx, geometry.columns());
                        if (geometry(row, col) == false)
                        {
                            // DIR 1
                            if(col < geometry.columns() - 1 && geometry(row, col + 1) == false)
                                _cells.at(idx)->add_neighbour(_cells.at(_coord2idx(col+1, row, geometry.columns())), 1);

                            // DIR 2
                            if(row > 0 && col < geometry.columns() - 1 && geometry(row - 1, col + 1) == false)
                                _cells.at(idx)->add_neighbour(_cells.at(_coord2idx(col+1, row-1, geometry.columns())), 2);

                            // DIR 3
                            if(row > 0 && geometry(row - 1, col) == false)
                                _cells.at(idx)->add_neighbour(_cells.at(_coord2idx(col, row-1, geometry.columns())), 3);

                            // DIR 4
                            if(row > 0 && col > 0 && geometry(row - 1, col - 1) == false)
                                _cells.at(idx)->add_neighbour(_cells.at(_coord2idx(col-1, row-1, geometry.columns())), 4);

                            // DIR 5
                            if(col > 0 && geometry(row, col - 1) == false)
                                _cells.at(idx)->add_neighbour(_cells.at(_coord2idx(col-1, row, geometry.columns())), 5);

                            // DIR 6
                            if(row < geometry.rows()  - 1&& col > 0 && geometry(row + 1, col - 1) == false)
                                _cells.at(idx)->add_neighbour(_cells.at(_coord2idx(col-1, row+1, geometry.columns())), 6);

                            // DIR 7
                            if(row < geometry.rows() - 1 && geometry(row + 1, col) == false)
                                _cells.at(idx)->add_neighbour(_cells.at(_coord2idx(col, row+1, geometry.columns())), 7);

                            // DIR 8
                            if(row < geometry.rows() - 1 && col < geometry.columns() - 1 && geometry(row + 1, col + 1) == false)
                                _cells.at(idx)->add_neighbour(_cells.at(_coord2idx(col+1, row+1, geometry.columns())), 8);
                        }
                    }

                    // remove obstacle cells
                    typename std::vector<Cell<DT_, directions> *>::iterator i = _cells.begin();
                    unsigned long id(0);
                    for (unsigned long idx(0) ; idx < geometry.size() ; ++idx)
                    {
                        unsigned long row(0);
                        unsigned long col(0);
                        _idx2coord(col, row, idx, geometry.columns());
                        if (geometry(row, col) == true)
                        {
                            i = _cells.erase(i);
                        }
                        else
                        {
                            (*i)->set_id(id);
                            ++i;
                            ++id;
                        }
                    }
                }
        };
}

#endif
