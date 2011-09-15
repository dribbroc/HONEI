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
#include <honei/util/private_implementation_pattern.hh>
#include <honei/util/private_implementation_pattern-impl.hh>

#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <iostream>


namespace honei
{
    template <typename DT_, unsigned long directions> class Grid3;

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
                    return Grid3<DT_, directions>::coord2idx(i->get_x(), i->get_y(), _max_col, _method) < Grid3<DT_, directions>::coord2idx(j->get_x(), j->get_y(), _max_col, _method);
                }
        };

    template <typename DT_, unsigned long directions>
    struct SyncData
    {
        DT_ data;
        unsigned long idx;
        signed long target_vector; // 1..9 f_n ; -1 h ; -2 u ; -3 v
        unsigned long process;
        Cell<DT_, directions> * cell;

        SyncData(DT_ _data, unsigned long _idx, signed long _target_vector, unsigned long _process, Cell<DT_, directions> * _cell):
            data(_data),
            idx(_idx),
            target_vector(_target_vector),
            process(_process),
            cell(_cell)
        {
        }
    };


    template <typename DT_, unsigned long directions> struct Implementation<Grid3<DT_, directions> >
    {
        std::vector<Cell<DT_, directions> *> cells; // all cells with ordering (inner|outer_halo|inner_halo)
        std::map<unsigned long, unsigned long> halo_map; // mapping of global idx -> local index into _cells
        std::vector<SyncData<DT_, directions> > send_targets; // list of all cells to send

        Implementation()
        {
        }

        ~Implementation()
        {
            for (unsigned long i(0) ; i < cells.size() ; ++i)
            {
                delete cells.at(i);
            }
        }
    };

    template <typename DT_, unsigned long directions>
        class Grid3 :
            public PrivateImplementationPattern<Grid3<DT_, directions>, Shared>
        {
            private:
                unsigned long _idx_start;
                unsigned long _idx_end;
                unsigned long _local_size;
                unsigned long _inner_halo_size;
                std::string _inner_numbering;
                std::string _outer_numbering;

                static bool _sync_data_comp(SyncData<DT_, directions> i, SyncData<DT_, directions> j)
                {
                    return i.process < j.process;
                }

                static unsigned long _idx2process(unsigned long idx, unsigned long process_count, unsigned long * /*starts*/, unsigned long * ends)
                {
                    for (unsigned long i(0) ; i < process_count ; ++i)
                        if (idx < ends[i])
                            return i;
                    return process_count - 1;
                }

                bool _is_valid_direction(unsigned long dir, unsigned long row, unsigned long col, unsigned long & nrow, unsigned long & ncol, DenseMatrix<bool> & geometry)
                {
                    switch(dir)
                    {
                        case 0:
                            throw InternalError("You should not check against direction zero!");
                            break;

                        case 1:
                            if(col < geometry.columns() - 1 && geometry(row, col + 1) == false)
                            {
                                ncol = col + 1;
                                nrow = row;
                                return true;
                            }
                            else
                                return false;
                            break;

                        case 2:
                            if(row > 0 && col < geometry.columns() - 1 && geometry(row - 1, col + 1) == false)
                            {
                                ncol = col + 1;
                                nrow = row - 1;
                                return true;
                            }
                            else
                                return false;
                            break;

                        case 3:
                            if(row > 0 && geometry(row - 1, col) == false)
                            {
                                ncol = col;
                                nrow = row - 1;
                                return true;
                            }
                            else
                                return false;
                            break;

                        case 4:
                            if(row > 0 && col > 0 && geometry(row - 1, col - 1) == false)
                            {
                                ncol = col - 1;
                                nrow = row - 1;
                                return true;
                            }
                            else
                                return false;
                            break;

                        case 5:
                            if(col > 0 && geometry(row, col - 1) == false)
                            {
                                ncol = col - 1;
                                nrow = row;
                                return true;
                            }
                            else
                                return false;
                            break;

                        case 6:
                            if(row < geometry.rows() - 1 && col > 0 && geometry(row + 1, col - 1) == false)
                            {
                                ncol = col - 1;
                                nrow = row + 1;
                                return true;
                            }
                            else
                                return false;
                            break;

                        case 7:
                            if(row < geometry.rows() - 1 && geometry(row + 1, col) == false)
                            {
                                ncol = col;
                                nrow = row + 1;
                                return true;
                            }
                            else
                                return false;
                            break;

                        case 8:
                            if(row < geometry.rows() - 1 && col < geometry.columns() - 1 && geometry(row + 1, col + 1) == false)
                            {
                                ncol = col + 1;
                                nrow = row + 1;
                                return true;
                            }
                            else
                                return false;
                            break;

                        default:
                            throw InternalError("You should not check against unknow direction!");
                            return false;
                            break;
                    }
                }

            public:
                // copy constructor
                Grid3(const Grid3<DT_, directions> & other) :
                    PrivateImplementationPattern<Grid3<DT_, directions>, Shared>(other._imp),
                    _idx_start(other._idx_start),
                    _idx_end(other._idx_end),
                    _local_size(other._local_size),
                    _inner_halo_size(other._inner_halo_size),
                    _inner_numbering(other._inner_numbering),
                    _outer_numbering(other._outer_numbering)
                {
                }

                ~Grid3()
                {
                }

                unsigned long size()
                {
                    return this->_imp->cells.size();
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
                    return this->_imp->cells.at(i);
                }

                std::vector<SyncData<DT_, directions> > & send_targets()
                {
                    return this->_imp->send_targets;
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

                void fill_h(DenseMatrix<DT_> & h, DenseMatrix<bool> & geometry, DenseVector<DT_> & h_v)
                {
                    unsigned long i(0);
                    for (unsigned long idx(0) ; idx < h.size() ; ++idx)
                    {
                        unsigned long row(0);
                        unsigned long col(0);
                        idx2coord(col, row, idx, geometry.columns(), _inner_numbering);
                        if (geometry(row, col) == false)
                        {
                            h(row, col) = h_v[i];
                            ++i;
                        }
                    }
                }

                static unsigned long coord2idx(unsigned long x, unsigned long y, unsigned long max_col, std::string method)
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


                static void idx2coord(unsigned long & x, unsigned long & y, unsigned long r, unsigned long max_col, std::string method)
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

                Grid3(DenseMatrix<bool> & geometry, DenseMatrix<DT_> & h, DenseMatrix<DT_> & b, DenseMatrix<DT_> & u,
                        DenseMatrix<DT_> & v, unsigned long process_id = 1, unsigned long process_count = 2) :
                    PrivateImplementationPattern<Grid3<DT_, directions>, Shared>(new Implementation<Grid3<DT_, directions> >())
                {
                    _outer_numbering = "z-curve";
                    //_outer_numbering = "row";
                    //_inner_numbering = "row";
                    _inner_numbering = "z-curve";

                    //TODO iteration twice over every global cell is lame
                    // will be solved, when the partition vector is served from outside :)

                    // calc global fluid cell count
                    unsigned long fluid_cells(0);
                    for (unsigned long idx(0) ; idx < geometry.size() ; ++idx)
                    {
                        unsigned long row(0);
                        unsigned long col(0);
                        idx2coord(col, row, idx, geometry.columns(), _outer_numbering);
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
                            idx2coord(col, row, idx, geometry.columns(), _outer_numbering);
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
                        idx2coord(col, row, idx, geometry.columns(), _outer_numbering);
                        Cell<DT_, directions>* cell = new Cell<DT_, directions>(col, row, 1, 1,
                                h(row, col), b(row, col), u(row, col), v(row, col));
                        this->_imp->cells.push_back(cell);
                        if (geometry(row, col) == false)
                            ++_local_size;
                    }

                    //resort cells by inner numbering
                    CellComparator<DT_, directions> cell_comp(geometry.columns(), _inner_numbering);
                    std::sort(this->_imp->cells.begin(), this->_imp->cells.end(), cell_comp);

                    // create temp cell map with global id keys
                    std::map<unsigned long, Cell<DT_, directions> *> temp_cells;
                    for (unsigned long i(0) ; i < this->_imp->cells.size() ; ++i)
                    {
                        temp_cells.insert(std::pair<unsigned long, Cell<DT_, directions> *>
                                (coord2idx(this->_imp->cells.at(i)->get_x(), this->_imp->cells.at(i)->get_y(), geometry.columns(), _inner_numbering), this->_imp->cells.at(i)));
                    }

                    // inner and outer halo
                    std::map<unsigned long, Cell<DT_, directions> *> halo;
                    std::set<Cell<DT_, directions> *> inner_halo;
                    typename std::map<unsigned long, Cell<DT_, directions> *>::iterator halo_it;

                    // set neighbourhood
                    for (typename std::vector<Cell<DT_, directions> *>::iterator i = this->_imp->cells.begin() ; i != this->_imp->cells.end() ; ++i)
                    {
                        unsigned long row((*i)->get_y());
                        unsigned long col((*i)->get_x());
                        unsigned long new_col(0);
                        unsigned long new_row(0);
                        if (geometry(row, col) == false)
                        {
                            std::set<unsigned long> h_targets;
                            unsigned long target_id(0);

                            for (unsigned long direction(1) ; direction < directions ; ++direction)
                            {
                                if(_is_valid_direction(direction, row, col, new_row, new_col, geometry))
                                {
                                    target_id = coord2idx(new_col, new_row, geometry.columns(), _inner_numbering);
                                    unsigned long outer_target_id(coord2idx(new_col, new_row, geometry.columns(), _outer_numbering));
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
                                            SyncData<DT_, directions> sync_data(0, outer_target_id, direction, _idx2process(outer_target_id, process_count, idx_starts, idx_ends), halo_it->second);
                                            this->_imp->send_targets.push_back(sync_data);
                                        }
                                        else
                                        {
                                            Cell<DT_, directions>* cell = new Cell<DT_, directions>(new_col, new_row, 1, 1,
                                                    h(new_row, new_col), b(new_row, new_col), u(new_row, new_col), v(new_row, new_col));
                                            halo[target_id] = cell;
                                            (*i)->add_neighbour(cell, direction);
                                            SyncData<DT_, directions> sync_data(0, outer_target_id, direction, _idx2process(outer_target_id, process_count, idx_starts, idx_ends), cell);
                                            this->_imp->send_targets.push_back(sync_data);
                                        }
                                    }
                                }
                            }

                            // copy inner halo to this->_imp->h_targets
                            for (typename std::set<unsigned long>::iterator t = h_targets.begin() ; t != h_targets.end() ; ++t)
                            {
                                SyncData<DT_, directions> sync_data(0, coord2idx(col, row, geometry.columns(), _outer_numbering), -1, *t, *i);
                                this->_imp->send_targets.push_back(sync_data);
                            }
                        }
                    }

                    // sort sync_data by target process
                    std::sort(this->_imp->send_targets.begin(), this->_imp->send_targets.end(), _sync_data_comp);

                    // insert outer halo cells in cell vector
                    for (halo_it = halo.begin() ; halo_it != halo.end() ; ++halo_it)
                    {
                        this->_imp->cells.push_back(halo_it->second);
                    }

                    _inner_halo_size = inner_halo.size();
                    _local_size -= _inner_halo_size;

                    std::cout<<"local size: "<<_local_size<<" inner halo size: "<<_inner_halo_size<< " outer halo size: "<<halo.size()<<std::endl;


                    // rearrange elements: inner cells in the front - inner halo at the end
                    for (typename std::set<Cell<DT_, directions> *>::iterator i = inner_halo.begin() ; i != inner_halo.end() ; ++i)
                    {
                        typename std::vector<Cell<DT_, directions> * >::iterator t = std::find(this->_imp->cells.begin(), this->_imp->cells.end(), *i);
                        this->_imp->cells.erase(t);
                        this->_imp->cells.push_back(*i);
                    }

                    // remove obstacle cells and enumerate fluid cells according to their array position
                    unsigned long id(0);
                    for (typename std::vector<Cell<DT_, directions> *>::iterator i = this->_imp->cells.begin() ; i != this->_imp->cells.end() ; )
                    {
                        unsigned long row((*i)->get_y());
                        unsigned long col((*i)->get_x());
                        if (geometry(row, col) == true)
                        {
                            i = this->_imp->cells.erase(i);
                        }
                        else
                        {
                            (*i)->set_id(id);
                            ++i;
                            ++id;
                        }
                    }

                    // fill this->_imp->halo_map
                    for (unsigned long i(_local_size) ; i < this->_imp->cells.size() ; ++i)
                    {
                        this->_imp->halo_map.insert(std::pair<unsigned long, unsigned long>(
                                    coord2idx(this->_imp->cells.at(i)->get_x(), this->_imp->cells.at(i)->get_y(), geometry.columns(), _outer_numbering), i));
                    }
                }
        };
}

#endif
