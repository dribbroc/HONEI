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
#include <map>
#include <algorithm>
#include <iostream>


namespace honei
{
    template <typename DT_, unsigned long directions>
        class Grid
        {
            private:
                std::vector<Cell<DT_, directions> *> _cells;
                unsigned long _idx_start;
                unsigned long _idx_end;
                unsigned long _local_size;
                unsigned long _inner_halo_size;

                static unsigned long _coord2idx(unsigned long x, unsigned long y, unsigned long max_col)
                {
                    std::string method("z-curve");
                    //std::string method("row");
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
                    //std::string method("row");
                    if (method.compare("z-curve") == 0)
                    {
                        unsigned long s(r);
                        x = 0;
                        y = 0;

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

                static unsigned long _idx2process(unsigned long idx, unsigned long process_count, unsigned long * /*starts*/, unsigned long * ends)
                {
                    for (unsigned long i(0) ; i < process_count ; ++i)
                        if (idx < ends[i])
                            return i;
                    return process_count - 1;
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

                unsigned long local_size()
                {
                    return _local_size;
                }

                unsigned long inner_halo_size()
                {
                    return _inner_halo_size;
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

                Grid(DenseMatrix<bool> & geometry, DenseMatrix<DT_> & h, DenseMatrix<DT_> & b, DenseMatrix<DT_> & u,
                        DenseMatrix<DT_> & v, unsigned long process_id = 1, unsigned long process_count = 3)
                {
                    //TODO iteration over every global cell and every single process is lame
                    // will be solved, when the partition vector is served from outside :)

                    // calc global fluid cell count
                    unsigned long fluid_cells(0);
                    for (unsigned long idx(0) ; idx < geometry.size() ; ++idx)
                    {
                        unsigned long row(0);
                        unsigned long col(0);
                        _idx2coord(col, row, idx, geometry.columns());
                        if (geometry(row, col) == false)
                            ++fluid_cells;
                    }

                    // find start and end of every patch
                    unsigned long idx_starts[process_count];
                    unsigned long idx_ends[process_count];
                    for (unsigned long process(0) ; process < process_count ; ++process)
                    {
                        unsigned long fluid_start(process * (fluid_cells / process_count));
                        unsigned long fluid_end(process == process_count - 1 ? fluid_cells : (process+1) * (fluid_cells / process_count));
                        unsigned long nfluid_cells = 0;
                        idx_starts[process] = 0;
                        idx_ends[process] = 0;
                        for (unsigned long idx(0) ; idx < geometry.size() ; ++idx)
                        {
                            unsigned long row(0);
                            unsigned long col(0);
                            _idx2coord(col, row, idx, geometry.columns());
                            if (geometry(row, col) == false)
                            {
                                if (nfluid_cells == fluid_start)
                                {
                                    idx_starts[process] = idx;
                                    if(fluid_start == 0)
                                        idx_starts[process] = 0;
                                }

                                if (nfluid_cells == fluid_end)
                                    idx_ends[process] = idx;
                                ++nfluid_cells;

                            }
                        }
                        if (idx_ends[process] == 0)
                            idx_ends[process] = geometry.size();
                    }
                    _idx_start = idx_starts[process_id];
                    _idx_end = idx_ends[process_id];


                    // read in cells
                    _local_size = 0;
                    for (unsigned long idx(_idx_start) ; idx < _idx_end ; ++idx)
                    {
                        unsigned long row(0);
                        unsigned long col(0);
                        _idx2coord(col, row, idx, geometry.columns());
                        Cell<DT_, directions>* cell = new Cell<DT_, directions>(col, row, 1, 1,
                                h(row, col), b(row, col), u(row, col), v(row, col));
                        this->add_cell(cell);
                        if (geometry(row, col) == false)
                            ++_local_size;
                    }

                    std::map<unsigned long, Cell<DT_, directions> *> halo;
                    std::map<unsigned long, Cell<DT_, directions> *> inner_halo;
                    typename std::map<unsigned long, Cell<DT_, directions> *>::iterator halo_it;
                    // set neighbourhood
                    for (unsigned long idx(_idx_start) ; idx < _idx_end ; ++idx)
                    {
                        unsigned long row(0);
                        unsigned long col(0);
                        _idx2coord(col, row, idx, geometry.columns());
                        if (geometry(row, col) == false)
                        {
                            unsigned long target_id(0);

                            // DIR 1
                            if(col < geometry.columns() - 1 && geometry(row, col + 1) == false)
                            {
                                target_id = _coord2idx(col+1, row, geometry.columns());
                                if(target_id >= _idx_start && target_id < _idx_end)
                                    _cells.at(idx - _idx_start)->add_neighbour(_cells.at(target_id - _idx_start), 1);
                                else
                                {
                                    inner_halo.insert(std::pair<unsigned long, Cell<DT_, directions> *>
                                            (idx, _cells.at(idx - _idx_start)));
                                    halo_it = halo.find(target_id);
                                    if(halo_it != halo.end())
                                    {
                                        _cells.at(idx - _idx_start)->add_neighbour(halo_it->second, 1);
                                    }
                                    else
                                    {
                                        unsigned long ncol(0);
                                        unsigned long nrow(0);
                                        _idx2coord(ncol, nrow, target_id, geometry.columns());
                                        Cell<DT_, directions>* cell = new Cell<DT_, directions>(ncol, nrow, 1, 1,
                                                h(nrow, ncol), b(nrow, ncol), u(nrow, ncol), v(nrow, ncol));
                                        halo[target_id] = cell;
                                        _cells.at(idx - _idx_start)->add_neighbour(cell, 1);
                                    }
                                }
                            }

                            // DIR 2
                            if(row > 0 && col < geometry.columns() - 1 && geometry(row - 1, col + 1) == false)
                            {
                                target_id = _coord2idx(col+1, row-1, geometry.columns());
                                if(target_id >= _idx_start && target_id < _idx_end)
                                    _cells.at(idx - _idx_start)->add_neighbour(_cells.at(target_id - _idx_start), 2);
                                else
                                {
                                    inner_halo.insert(std::pair<unsigned long, Cell<DT_, directions> *>
                                            (idx, _cells.at(idx - _idx_start)));
                                    halo_it = halo.find(target_id);
                                    if(halo_it != halo.end())
                                    {
                                        _cells.at(idx - _idx_start)->add_neighbour(halo_it->second, 2);
                                    }
                                    else
                                    {
                                        unsigned long ncol(0);
                                        unsigned long nrow(0);
                                        _idx2coord(ncol, nrow, target_id, geometry.columns());
                                        Cell<DT_, directions>* cell = new Cell<DT_, directions>(ncol, nrow, 1, 1,
                                                h(nrow, ncol), b(nrow, ncol), u(nrow, ncol), v(nrow, ncol));
                                        halo[target_id] = cell;
                                        _cells.at(idx - _idx_start)->add_neighbour(cell, 2);
                                    }
                                }
                            }

                            // DIR 3
                            if(row > 0 && geometry(row - 1, col) == false)
                            {
                                target_id = _coord2idx(col, row-1, geometry.columns());
                                if(target_id >= _idx_start && target_id < _idx_end)
                                    _cells.at(idx - _idx_start)->add_neighbour(_cells.at(target_id - _idx_start), 3);
                                else
                                {
                                    inner_halo.insert(std::pair<unsigned long, Cell<DT_, directions> *>
                                            (idx, _cells.at(idx - _idx_start)));
                                    halo_it = halo.find(target_id);
                                    if(halo_it != halo.end())
                                    {
                                        _cells.at(idx - _idx_start)->add_neighbour(halo_it->second, 3);
                                    }
                                    else
                                    {
                                        unsigned long ncol(0);
                                        unsigned long nrow(0);
                                        _idx2coord(ncol, nrow, target_id, geometry.columns());
                                        Cell<DT_, directions>* cell = new Cell<DT_, directions>(ncol, nrow, 1, 1,
                                                h(nrow, ncol), b(nrow, ncol), u(nrow, ncol), v(nrow, ncol));
                                        halo[target_id] = cell;
                                        _cells.at(idx - _idx_start)->add_neighbour(cell, 3);
                                    }
                                }
                            }

                            // DIR 4
                            if(row > 0 && col > 0 && geometry(row - 1, col - 1) == false)
                            {
                                target_id = _coord2idx(col-1, row-1, geometry.columns());
                                if(target_id >= _idx_start && target_id < _idx_end)
                                    _cells.at(idx - _idx_start)->add_neighbour(_cells.at(target_id - _idx_start), 4);
                                else
                                {
                                    inner_halo.insert(std::pair<unsigned long, Cell<DT_, directions> *>
                                            (idx, _cells.at(idx - _idx_start)));
                                    halo_it = halo.find(target_id);
                                    if(halo_it != halo.end())
                                    {
                                        _cells.at(idx - _idx_start)->add_neighbour(halo_it->second, 4);
                                    }
                                    else
                                    {
                                        unsigned long ncol(0);
                                        unsigned long nrow(0);
                                        _idx2coord(ncol, nrow, target_id, geometry.columns());
                                        Cell<DT_, directions>* cell = new Cell<DT_, directions>(ncol, nrow, 1, 1,
                                                h(nrow, ncol), b(nrow, ncol), u(nrow, ncol), v(nrow, ncol));
                                        halo[target_id] = cell;
                                        _cells.at(idx - _idx_start)->add_neighbour(cell, 4);
                                    }
                                }
                            }

                            // DIR 5
                            if(col > 0 && geometry(row, col - 1) == false)
                            {
                                target_id = _coord2idx(col-1, row, geometry.columns());
                                if(target_id >= _idx_start && target_id < _idx_end)
                                    _cells.at(idx - _idx_start)->add_neighbour(_cells.at(target_id - _idx_start), 5);
                                else
                                {
                                    inner_halo.insert(std::pair<unsigned long, Cell<DT_, directions> *>
                                            (idx, _cells.at(idx - _idx_start)));
                                    halo_it = halo.find(target_id);
                                    if(halo_it != halo.end())
                                    {
                                        _cells.at(idx - _idx_start)->add_neighbour(halo_it->second, 5);
                                    }
                                    else
                                    {
                                        unsigned long ncol(0);
                                        unsigned long nrow(0);
                                        _idx2coord(ncol, nrow, target_id, geometry.columns());
                                        Cell<DT_, directions>* cell = new Cell<DT_, directions>(ncol, nrow, 1, 1,
                                                h(nrow, ncol), b(nrow, ncol), u(nrow, ncol), v(nrow, ncol));
                                        halo[target_id] = cell;
                                        _cells.at(idx - _idx_start)->add_neighbour(cell, 5);
                                    }
                                }
                            }

                            // DIR 6
                            if(row < geometry.rows()  - 1&& col > 0 && geometry(row + 1, col - 1) == false)
                            {
                                target_id = _coord2idx(col-1, row+1, geometry.columns());
                                if(target_id >= _idx_start && target_id < _idx_end)
                                    _cells.at(idx - _idx_start)->add_neighbour(_cells.at(target_id - _idx_start), 6);
                                else
                                {
                                    inner_halo.insert(std::pair<unsigned long, Cell<DT_, directions> *>
                                            (idx, _cells.at(idx - _idx_start)));
                                    halo_it = halo.find(target_id);
                                    if(halo_it != halo.end())
                                    {
                                        _cells.at(idx - _idx_start)->add_neighbour(halo_it->second, 6);
                                    }
                                    else
                                    {
                                        unsigned long ncol(0);
                                        unsigned long nrow(0);
                                        _idx2coord(ncol, nrow, target_id, geometry.columns());
                                        Cell<DT_, directions>* cell = new Cell<DT_, directions>(ncol, nrow, 1, 1,
                                                h(nrow, ncol), b(nrow, ncol), u(nrow, ncol), v(nrow, ncol));
                                        halo[target_id] = cell;
                                        _cells.at(idx - _idx_start)->add_neighbour(cell, 6);
                                    }
                                }
                            }

                            // DIR 7
                            if(row < geometry.rows() - 1 && geometry(row + 1, col) == false)
                            {
                                target_id = _coord2idx(col, row+1, geometry.columns());
                                if(target_id >= _idx_start && target_id < _idx_end)
                                    _cells.at(idx - _idx_start)->add_neighbour(_cells.at(target_id - _idx_start), 7);
                                else
                                {
                                    inner_halo.insert(std::pair<unsigned long, Cell<DT_, directions> *>
                                            (idx, _cells.at(idx - _idx_start)));
                                    halo_it = halo.find(target_id);
                                    if(halo_it != halo.end())
                                    {
                                        _cells.at(idx - _idx_start)->add_neighbour(halo_it->second, 7);
                                    }
                                    else
                                    {
                                        unsigned long ncol(0);
                                        unsigned long nrow(0);
                                        _idx2coord(ncol, nrow, target_id, geometry.columns());
                                        Cell<DT_, directions>* cell = new Cell<DT_, directions>(ncol, nrow, 1, 1,
                                                h(nrow, ncol), b(nrow, ncol), u(nrow, ncol), v(nrow, ncol));
                                        halo[target_id] = cell;
                                        _cells.at(idx - _idx_start)->add_neighbour(cell, 7);
                                    }
                                }
                            }

                            // DIR 8
                            if(row < geometry.rows() - 1 && col < geometry.columns() - 1 && geometry(row + 1, col + 1) == false)
                            {
                                target_id = _coord2idx(col+1, row+1, geometry.columns());
                                if(target_id >= _idx_start && target_id < _idx_end)
                                    _cells.at(idx - _idx_start)->add_neighbour(_cells.at(target_id - _idx_start), 8);
                                else
                                {
                                    inner_halo.insert(std::pair<unsigned long, Cell<DT_, directions> *>
                                            (idx, _cells.at(idx - _idx_start)));
                                    halo_it = halo.find(target_id);
                                    if(halo_it != halo.end())
                                    {
                                        _cells.at(idx - _idx_start)->add_neighbour(halo_it->second, 8);
                                    }
                                    else
                                    {
                                        unsigned long ncol(0);
                                        unsigned long nrow(0);
                                        _idx2coord(ncol, nrow, target_id, geometry.columns());
                                        Cell<DT_, directions>* cell = new Cell<DT_, directions>(ncol, nrow, 1, 1,
                                                h(nrow, ncol), b(nrow, ncol), u(nrow, ncol), v(nrow, ncol));
                                        halo[target_id] = cell;
                                        _cells.at(idx - _idx_start)->add_neighbour(cell, 8);
                                    }
                                }
                            }
                        }
                    }

                    std::cout<<"Halo:"<<std::endl;
                    // insert halo cells in cell vector
                    for (halo_it = halo.begin() ; halo_it != halo.end() ; ++halo_it)
                    {
                        std::cout<<halo_it->first<<" "<<halo_it->second->get_y()<<" "<<halo_it->second->get_x()<<std::endl;
                        this->add_cell(halo_it->second);
                    }
                    std::cout<<"Innner Halo:"<<std::endl;
                    // insert halo cells in cell vector
                    for (halo_it = inner_halo.begin() ; halo_it != inner_halo.end() ; ++halo_it)
                    {
                        std::cout<<halo_it->first<<" "<<halo_it->second->get_y()<<" "<<halo_it->second->get_x()<<std::endl;
                    }
                    std::cout<<"offset: "<<_idx_start<<" end: "<<_idx_end<<std::endl;

                    _inner_halo_size = inner_halo.size();
                    _local_size -= _inner_halo_size;

                    std::cout<<"local size: "<<_local_size<<" inner halo size: "<<_inner_halo_size<<std::endl;

                    // rearrange elements: inner cells in the front - inner halo at the end
                    for (halo_it = inner_halo.begin() ; halo_it != inner_halo.end() ; ++halo_it)
                    {
                        Cell<DT_, directions> * temp = halo_it->second;
                        typename std::vector<Cell<DT_, directions> * >::iterator t = std::find(_cells.begin(), _cells.end(), temp);
                        _cells.erase(t);
                        _cells.push_back(temp);
                    }

                    // remove obstacle cells and enumerate fluid cells according to their array position
                    unsigned long id(0);
                    for (typename std::vector<Cell<DT_, directions> *>::iterator i = _cells.begin() ; i != _cells.end() ; )
                    {
                        unsigned long row((*i)->get_y());
                        unsigned long col((*i)->get_x());
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
