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
#include <set>
#include <algorithm>
#include <iostream>


namespace honei
{
    unsigned long coord2idx(unsigned long x, unsigned long y, unsigned long max_col, std::string method)
    {
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


    void idx2coord(unsigned long & x, unsigned long & y, unsigned long r, unsigned long max_col, std::string method)
    {
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

    template <typename DT_, unsigned long directions>
        class CellComparator
        {
            private:
                unsigned long _max_col;
                std::string _method;
            public:
                CellComparator(unsigned long max_col, std::string method) :
                    _max_col(max_col),
                    _method(method)
                {
                }

                bool operator() (Cell<DT_, directions> * i, Cell<DT_, directions> * j)
                {
                    return coord2idx(i->get_x(), i->get_y(), _max_col, _method) < coord2idx(j->get_x(), j->get_y(), _max_col, _method);
                }
        };


    template <typename DT_, unsigned long directions>
        class Grid
        {
            private:
                std::vector<Cell<DT_, directions> *> _cells;
                std::multimap<unsigned long, Cell<DT_, directions> *> _h_targets;
                std::multimap<unsigned long, Cell<DT_, directions> *> _dir_targets[directions];
                unsigned long _idx_start;
                unsigned long _idx_end;
                unsigned long _local_size;
                unsigned long _inner_halo_size;

                unsigned long _idx2process(unsigned long idx, unsigned long process_count, unsigned long * /*starts*/, unsigned long * ends)
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

                Cell<DT_, directions> * get_cell(unsigned long i)
                {
                    return _cells.at(i);
                }

                static void print_numbering(DenseMatrix<bool> & geometry, std::string method)
                {
                    DenseMatrix<long> result(geometry.rows(), geometry.columns(), -1);
                    unsigned long i(0);
                    for (unsigned long idx(0) ; idx < result.size() ; ++idx)
                    {
                        unsigned long row(0);
                        unsigned long col(0);
                        idx2coord(col, row, idx, geometry.columns(), method);
                        if (geometry(row, col) == false)
                        {
                            result(row, col) = i;
                            ++i;
                        }
                    }
                    std::cout<<result;
                }

                Grid(DenseMatrix<bool> & geometry, DenseMatrix<DT_> & h, DenseMatrix<DT_> & b, DenseMatrix<DT_> & u,
                        DenseMatrix<DT_> & v, unsigned long process_id = 0, unsigned long process_count = 1)
                {
                    std::string outer_numbering("z-curve");
                    //std::string outer_numbering("row");
                    std::string inner_numbering("row");
                    //std::string inner_numbering("z-curve");

                    //TODO iteration twice over every global cell is lame
                    // will be solved, when the partition vector is served from outside :)

                    // calc global fluid cell count
                    unsigned long fluid_cells(0);
                    for (unsigned long idx(0) ; idx < geometry.size() ; ++idx)
                    {
                        unsigned long row(0);
                        unsigned long col(0);
                        idx2coord(col, row, idx, geometry.columns(), outer_numbering);
                        if (geometry(row, col) == false)
                            ++fluid_cells;
                    }

                    // find start and end of every patch in idx coords, depending on fluid count
                    unsigned long idx_starts[process_count];
                    unsigned long idx_ends[process_count];
                    unsigned long nfluid_cells = 0;
                    for (unsigned long process(0) ; process < process_count ; ++process)
                    {
                        unsigned long fluid_start(process * (fluid_cells / process_count));
                        unsigned long fluid_end(process == process_count - 1 ? fluid_cells : (process+1) * (fluid_cells / process_count));
                        idx_starts[process] = 0;
                        idx_ends[process] = 0;
                        for (unsigned long idx(idx_ends[process == 0 ? 0 : process - 1]) ; idx < geometry.size() ; ++idx)
                        {
                            unsigned long row(0);
                            unsigned long col(0);
                            idx2coord(col, row, idx, geometry.columns(), outer_numbering);
                            if (geometry(row, col) == false)
                            {
                                if (nfluid_cells == fluid_start)
                                {
                                    idx_starts[process] = idx;
                                    if(fluid_start == 0)
                                        idx_starts[process] = 0;
                                }

                                if (nfluid_cells == fluid_end)
                                {
                                    idx_ends[process] = idx;
                                    break;
                                }
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
                        idx2coord(col, row, idx, geometry.columns(), outer_numbering);
                        Cell<DT_, directions>* cell = new Cell<DT_, directions>(col, row, 1, 1,
                                h(row, col), b(row, col), u(row, col), v(row, col));
                        _cells.push_back(cell);
                        if (geometry(row, col) == false)
                            ++_local_size;
                    }

                    //resort cells by inner numbering
                    CellComparator<DT_, directions> cell_comp(geometry.columns(), inner_numbering);
                    std::sort(_cells.begin(), _cells.end(), cell_comp);

                    // create temp cell map with global id keys
                    std::map<unsigned long, Cell<DT_, directions> *> temp_cells;
                    for (unsigned long i(0) ; i < _cells.size() ; ++i)
                    {
                        temp_cells.insert(std::pair<unsigned long, Cell<DT_, directions> *>
                                (coord2idx(_cells.at(i)->get_x(), _cells.at(i)->get_y(), geometry.columns(), inner_numbering), _cells.at(i)));
                    }

                    // inner and outer halo
                    std::map<unsigned long, Cell<DT_, directions> *> halo;
                    std::set<Cell<DT_, directions> *> inner_halo;
                    typename std::map<unsigned long, Cell<DT_, directions> *>::iterator halo_it;

                    // set neighbourhood
                    for (typename std::vector<Cell<DT_, directions> *>::iterator i = _cells.begin() ; i != _cells.end() ; ++i)
                    {
                        unsigned long row((*i)->get_y());
                        unsigned long col((*i)->get_x());
                        if (geometry(row, col) == false)
                        {
                            std::set<unsigned long> h_targets;
                            unsigned long target_id(0);
                            unsigned long direction(0);

                            // TODO umschreiben: get row col auf directino methode und dann nur ein block in for 1..9 schleife
                            // und nachbar id auch nur einmal berechnen pro direction
                            // DIR 1
                            if(col < geometry.columns() - 1 && geometry(row, col + 1) == false)
                            {
                                direction = 1;
                                target_id = coord2idx(col+1, row, geometry.columns(), inner_numbering);
                                unsigned long outer_target_id(coord2idx(col+1, row, geometry.columns(), outer_numbering));
                                // if our neighbour is another "normal" cell
                                if(_idx2process(outer_target_id, process_count, idx_starts, idx_ends) == process_id)
                                    (*i)->add_neighbour(temp_cells[target_id], direction);
                                // our neighbour lies outside
                                else
                                {
                                    // mark cell itself as having outer neighbours
                                    inner_halo.insert(*i);
                                    h_targets.insert(_idx2process(outer_target_id, process_count, idx_starts, idx_ends));
                                    halo_it = halo.find(target_id);
                                    if(halo_it != halo.end())
                                    {
                                        (*i)->add_neighbour(halo_it->second, direction);
                                        _dir_targets[direction].insert(std::pair<unsigned long, Cell<DT_, directions> *>(_idx2process(outer_target_id, process_count, idx_starts, idx_ends), halo_it->second));
                                    }
                                    else
                                    {
                                        unsigned long ncol(0);
                                        unsigned long nrow(0);
                                        idx2coord(ncol, nrow, target_id, geometry.columns(), inner_numbering);
                                        Cell<DT_, directions>* cell = new Cell<DT_, directions>(ncol, nrow, 1, 1,
                                                h(nrow, ncol), b(nrow, ncol), u(nrow, ncol), v(nrow, ncol));
                                        halo[target_id] = cell;
                                        (*i)->add_neighbour(cell, direction);
                                        _dir_targets[direction].insert(std::pair<unsigned long, Cell<DT_, directions> *>(_idx2process(outer_target_id, process_count, idx_starts, idx_ends), cell));
                                    }
                                }
                            }

                            // DIR 2
                            if(row > 0 && col < geometry.columns() - 1 && geometry(row - 1, col + 1) == false)
                            {
                                direction = 2;
                                target_id = coord2idx(col+1, row-1, geometry.columns(), inner_numbering);
                                unsigned long outer_target_id(coord2idx(col+1, row-1, geometry.columns(), outer_numbering));
                                // if our neighbour is another "normal" cell
                                if(_idx2process(outer_target_id, process_count, idx_starts, idx_ends) == process_id)
                                    (*i)->add_neighbour(temp_cells[target_id], direction);
                                // our neighbour lies outside
                                else
                                {
                                    // mark cell itself as having outer neighbours
                                    inner_halo.insert(*i);
                                    h_targets.insert(_idx2process(outer_target_id, process_count, idx_starts, idx_ends));
                                    halo_it = halo.find(target_id);
                                    if(halo_it != halo.end())
                                    {
                                        (*i)->add_neighbour(halo_it->second, direction);
                                        _dir_targets[direction].insert(std::pair<unsigned long, Cell<DT_, directions> *>(_idx2process(outer_target_id, process_count, idx_starts, idx_ends), halo_it->second));
                                    }
                                    else
                                    {
                                        unsigned long ncol(0);
                                        unsigned long nrow(0);
                                        idx2coord(ncol, nrow, target_id, geometry.columns(), inner_numbering);
                                        Cell<DT_, directions>* cell = new Cell<DT_, directions>(ncol, nrow, 1, 1,
                                                h(nrow, ncol), b(nrow, ncol), u(nrow, ncol), v(nrow, ncol));
                                        halo[target_id] = cell;
                                        (*i)->add_neighbour(cell, direction);
                                        _dir_targets[direction].insert(std::pair<unsigned long, Cell<DT_, directions> *>(_idx2process(outer_target_id, process_count, idx_starts, idx_ends), cell));
                                    }
                                }
                            }

                            // DIR 3
                            if(row > 0 && geometry(row - 1, col) == false)
                            {
                                direction = 3;
                                target_id = coord2idx(col, row-1, geometry.columns(), inner_numbering);
                                unsigned long outer_target_id(coord2idx(col, row-1, geometry.columns(), outer_numbering));
                                // if our neighbour is another "normal" cell
                                if(_idx2process(outer_target_id, process_count, idx_starts, idx_ends) == process_id)
                                    (*i)->add_neighbour(temp_cells[target_id], direction);
                                // our neighbour lies outside
                                else
                                {
                                    // mark cell itself as having outer neighbours
                                    inner_halo.insert(*i);
                                    h_targets.insert(_idx2process(outer_target_id, process_count, idx_starts, idx_ends));
                                    halo_it = halo.find(target_id);
                                    if(halo_it != halo.end())
                                    {
                                        (*i)->add_neighbour(halo_it->second, direction);
                                        _dir_targets[direction].insert(std::pair<unsigned long, Cell<DT_, directions> *>(_idx2process(outer_target_id, process_count, idx_starts, idx_ends), halo_it->second));
                                    }
                                    else
                                    {
                                        unsigned long ncol(0);
                                        unsigned long nrow(0);
                                        idx2coord(ncol, nrow, target_id, geometry.columns(), inner_numbering);
                                        Cell<DT_, directions>* cell = new Cell<DT_, directions>(ncol, nrow, 1, 1,
                                                h(nrow, ncol), b(nrow, ncol), u(nrow, ncol), v(nrow, ncol));
                                        halo[target_id] = cell;
                                        (*i)->add_neighbour(cell, direction);
                                        _dir_targets[direction].insert(std::pair<unsigned long, Cell<DT_, directions> *>(_idx2process(outer_target_id, process_count, idx_starts, idx_ends), cell));
                                    }
                                }
                            }

                            // DIR 4
                            if(row > 0 && col > 0 && geometry(row - 1, col - 1) == false)
                            {
                                direction = 4;
                                target_id = coord2idx(col-1, row-1, geometry.columns(), inner_numbering);
                                unsigned long outer_target_id(coord2idx(col-1, row-1, geometry.columns(), outer_numbering));
                                // if our neighbour is another "normal" cell
                                if(_idx2process(outer_target_id, process_count, idx_starts, idx_ends) == process_id)
                                    (*i)->add_neighbour(temp_cells[target_id], direction);
                                // our neighbour lies outside
                                else
                                {
                                    // mark cell itself as having outer neighbours
                                    inner_halo.insert(*i);
                                    h_targets.insert(_idx2process(outer_target_id, process_count, idx_starts, idx_ends));
                                    halo_it = halo.find(target_id);
                                    if(halo_it != halo.end())
                                    {
                                        (*i)->add_neighbour(halo_it->second, direction);
                                        _dir_targets[direction].insert(std::pair<unsigned long, Cell<DT_, directions> *>(_idx2process(outer_target_id, process_count, idx_starts, idx_ends), halo_it->second));
                                    }
                                    else
                                    {
                                        unsigned long ncol(0);
                                        unsigned long nrow(0);
                                        idx2coord(ncol, nrow, target_id, geometry.columns(), inner_numbering);
                                        Cell<DT_, directions>* cell = new Cell<DT_, directions>(ncol, nrow, 1, 1,
                                                h(nrow, ncol), b(nrow, ncol), u(nrow, ncol), v(nrow, ncol));
                                        halo[target_id] = cell;
                                        (*i)->add_neighbour(cell, direction);
                                        _dir_targets[direction].insert(std::pair<unsigned long, Cell<DT_, directions> *>(_idx2process(outer_target_id, process_count, idx_starts, idx_ends), cell));
                                    }
                                }
                            }

                            // DIR 5
                            if(col > 0 && geometry(row, col - 1) == false)
                            {
                                direction = 5;
                                target_id = coord2idx(col-1, row, geometry.columns(), inner_numbering);
                                unsigned long outer_target_id(coord2idx(col-1, row, geometry.columns(), outer_numbering));
                                // if our neighbour is another "normal" cell
                                if(_idx2process(outer_target_id, process_count, idx_starts, idx_ends) == process_id)
                                    (*i)->add_neighbour(temp_cells[target_id], direction);
                                // our neighbour lies outside
                                else
                                {
                                    // mark cell itself as having outer neighbours
                                    inner_halo.insert(*i);
                                    h_targets.insert(_idx2process(outer_target_id, process_count, idx_starts, idx_ends));
                                    halo_it = halo.find(target_id);
                                    if(halo_it != halo.end())
                                    {
                                        (*i)->add_neighbour(halo_it->second, direction);
                                        _dir_targets[direction].insert(std::pair<unsigned long, Cell<DT_, directions> *>(_idx2process(outer_target_id, process_count, idx_starts, idx_ends), halo_it->second));
                                    }
                                    else
                                    {
                                        unsigned long ncol(0);
                                        unsigned long nrow(0);
                                        idx2coord(ncol, nrow, target_id, geometry.columns(), inner_numbering);
                                        Cell<DT_, directions>* cell = new Cell<DT_, directions>(ncol, nrow, 1, 1,
                                                h(nrow, ncol), b(nrow, ncol), u(nrow, ncol), v(nrow, ncol));
                                        halo[target_id] = cell;
                                        (*i)->add_neighbour(cell, direction);
                                        _dir_targets[direction].insert(std::pair<unsigned long, Cell<DT_, directions> *>(_idx2process(outer_target_id, process_count, idx_starts, idx_ends), cell));
                                    }
                                }
                            }

                            // DIR 6
                            if(row < geometry.rows()  - 1&& col > 0 && geometry(row + 1, col - 1) == false)
                            {
                                direction = 6;
                                target_id = coord2idx(col-1, row+1, geometry.columns(), inner_numbering);
                                unsigned long outer_target_id(coord2idx(col-1, row+1, geometry.columns(), outer_numbering));
                                // if our neighbour is another "normal" cell
                                if(_idx2process(outer_target_id, process_count, idx_starts, idx_ends) == process_id)
                                    (*i)->add_neighbour(temp_cells[target_id], direction);
                                // our neighbour lies outside
                                else
                                {
                                    // mark cell itself as having outer neighbours
                                    inner_halo.insert(*i);
                                    h_targets.insert(_idx2process(outer_target_id, process_count, idx_starts, idx_ends));
                                    halo_it = halo.find(target_id);
                                    if(halo_it != halo.end())
                                    {
                                        (*i)->add_neighbour(halo_it->second, direction);
                                        _dir_targets[direction].insert(std::pair<unsigned long, Cell<DT_, directions> *>(_idx2process(outer_target_id, process_count, idx_starts, idx_ends), halo_it->second));
                                    }
                                    else
                                    {
                                        unsigned long ncol(0);
                                        unsigned long nrow(0);
                                        idx2coord(ncol, nrow, target_id, geometry.columns(), inner_numbering);
                                        Cell<DT_, directions>* cell = new Cell<DT_, directions>(ncol, nrow, 1, 1,
                                                h(nrow, ncol), b(nrow, ncol), u(nrow, ncol), v(nrow, ncol));
                                        halo[target_id] = cell;
                                        (*i)->add_neighbour(cell, direction);
                                        _dir_targets[direction].insert(std::pair<unsigned long, Cell<DT_, directions> *>(_idx2process(outer_target_id, process_count, idx_starts, idx_ends), cell));
                                    }
                                }
                            }

                            // DIR 7
                            if(row < geometry.rows() - 1 && geometry(row + 1, col) == false)
                            {
                                direction = 7;
                                target_id = coord2idx(col, row+1, geometry.columns(), inner_numbering);
                                unsigned long outer_target_id(coord2idx(col, row+1, geometry.columns(), outer_numbering));
                                // if our neighbour is another "normal" cell
                                if(_idx2process(outer_target_id, process_count, idx_starts, idx_ends) == process_id)
                                    (*i)->add_neighbour(temp_cells[target_id], direction);
                                // our neighbour lies outside
                                else
                                {
                                    // mark cell itself as having outer neighbours
                                    inner_halo.insert(*i);
                                    h_targets.insert(_idx2process(outer_target_id, process_count, idx_starts, idx_ends));
                                    halo_it = halo.find(target_id);
                                    if(halo_it != halo.end())
                                    {
                                        (*i)->add_neighbour(halo_it->second, direction);
                                        _dir_targets[direction].insert(std::pair<unsigned long, Cell<DT_, directions> *>(_idx2process(outer_target_id, process_count, idx_starts, idx_ends), halo_it->second));
                                    }
                                    else
                                    {
                                        unsigned long ncol(0);
                                        unsigned long nrow(0);
                                        idx2coord(ncol, nrow, target_id, geometry.columns(), inner_numbering);
                                        Cell<DT_, directions>* cell = new Cell<DT_, directions>(ncol, nrow, 1, 1,
                                                h(nrow, ncol), b(nrow, ncol), u(nrow, ncol), v(nrow, ncol));
                                        halo[target_id] = cell;
                                        (*i)->add_neighbour(cell, direction);
                                        _dir_targets[direction].insert(std::pair<unsigned long, Cell<DT_, directions> *>(_idx2process(outer_target_id, process_count, idx_starts, idx_ends), cell));
                                    }
                                }
                            }

                            // DIR 8
                            if(row < geometry.rows() - 1 && col < geometry.columns() - 1 && geometry(row + 1, col + 1) == false)
                            {
                                direction = 8;
                                target_id = coord2idx(col+1, row+1, geometry.columns(), inner_numbering);
                                unsigned long outer_target_id(coord2idx(col+1, row+1, geometry.columns(), outer_numbering));
                                // if our neighbour is another "normal" cell
                                if(_idx2process(outer_target_id, process_count, idx_starts, idx_ends) == process_id)
                                    (*i)->add_neighbour(temp_cells[target_id], direction);
                                // our neighbour lies outside
                                else
                                {
                                    // mark cell itself as having outer neighbours
                                    inner_halo.insert(*i);
                                    h_targets.insert(_idx2process(outer_target_id, process_count, idx_starts, idx_ends));
                                    halo_it = halo.find(target_id);
                                    if(halo_it != halo.end())
                                    {
                                        (*i)->add_neighbour(halo_it->second, direction);
                                        _dir_targets[direction].insert(std::pair<unsigned long, Cell<DT_, directions> *>(_idx2process(outer_target_id, process_count, idx_starts, idx_ends), halo_it->second));
                                    }
                                    else
                                    {
                                        unsigned long ncol(0);
                                        unsigned long nrow(0);
                                        idx2coord(ncol, nrow, target_id, geometry.columns(), inner_numbering);
                                        Cell<DT_, directions>* cell = new Cell<DT_, directions>(ncol, nrow, 1, 1,
                                                h(nrow, ncol), b(nrow, ncol), u(nrow, ncol), v(nrow, ncol));
                                        halo[target_id] = cell;
                                        (*i)->add_neighbour(cell, direction);
                                        _dir_targets[direction].insert(std::pair<unsigned long, Cell<DT_, directions> *>(_idx2process(outer_target_id, process_count, idx_starts, idx_ends), cell));
                                    }
                                }
                            }
                            // copy inner halo to _h_targets
                            for (typename std::set<unsigned long>::iterator t = h_targets.begin() ; t != h_targets.end() ; ++t)
                            {
                                _h_targets.insert(std::pair<unsigned long, Cell<DT_, directions> *>(*t, *i));
                            }
                        }
                    }


                    // insert outer halo cells in cell vector
                    for (halo_it = halo.begin() ; halo_it != halo.end() ; ++halo_it)
                    {
                        _cells.push_back(halo_it->second);
                    }

                    _inner_halo_size = inner_halo.size();
                    _local_size -= _inner_halo_size;

                    std::cout<<"local size: "<<_local_size<<" inner halo size: "<<_inner_halo_size<< " outer halo size: "<<halo.size()<<std::endl;


                    // rearrange elements: inner cells in the front - inner halo at the end
                    for (typename std::set<Cell<DT_, directions> *>::iterator i = inner_halo.begin() ; i != inner_halo.end() ; ++i)
                    {
                        typename std::vector<Cell<DT_, directions> * >::iterator t = std::find(_cells.begin(), _cells.end(), *i);
                        _cells.erase(t);
                        _cells.push_back(*i);
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
