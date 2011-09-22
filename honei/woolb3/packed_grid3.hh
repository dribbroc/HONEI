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
#ifndef WOOLB3_GUARD_PACKED_GRID_HH
#define WOOLB3_GUARD_PACKED_GRID_HH 1

#include <honei/woolb3/cell.hh>
#include <honei/woolb3/grid3.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/util/shared_array-impl.hh>

#include <vector>
#include <iostream>


namespace honei
{
    template <typename DT_>
    struct SyncTupple
    {
        DT_ data;
        unsigned long idx;
        long target_vector;

        SyncTupple(DT_ _data, unsigned long _idx, long _target_vector) :
            data(_data),
            idx(_idx),
            target_vector(_target_vector)
        {
        }
    };

    template <typename DT_>
    struct SyncData
    {
        std::vector<SyncTupple<DT_> > data;
        unsigned long process;
    };

    template <typename DT_, unsigned long directions>
    class PackedGrid3
    {
        public:
            void import_synch_data(std::list<SyncData<DT_> > & data)
            {
                for (typename std::list<SyncData<DT_> >::iterator i(data.begin()) ; i != data.end() ; ++i)
                {
                    for (typename std::vector<SyncTupple<DT_> >::iterator j(i->data.begin()) ; j != i->data.end() ; ++j)
                    {
                        if (j->target_vector == -1)
                            (*h)[grid.halo_map()[j->idx]] = j->data;
                        else if (j->target_vector > 0 && j->target_vector < (long) directions)
                            (*f_temp[j->target_vector])[grid.halo_map()[j->idx]] = j->data;
                        else
                            throw InternalError("Wrong target_vector in SyncTupple found!");
                    }
                }
            }

            std::list<SyncData<DT_> > export_synch_data()
            {
                unsigned long process(grid.send_targets().front().process);
                std::list<SyncData<DT_> > sync_list;
                SyncData<DT_> sync_data;
                sync_data.process = process;
                sync_list.push_back(sync_data);

                for (typename std::vector<SyncInfo<DT_, directions> >::iterator i(grid.send_targets().begin()) ; i != grid.send_targets().end() ; ++i)
                {
                    if (i->process != process)
                    {
                        process = i->process;
                        SyncData<DT_> sync_data;
                        sync_data.process = process;
                        sync_list.push_back(sync_data);
                    }

                    if (i->target_vector == -1)
                    {
                        SyncTupple<DT_> sync_tupple((*h2)[i->cell->get_id()], i->idx, i->target_vector);
                        sync_list.back().data.push_back(sync_tupple);
                    }
                    else if (i->target_vector >= 0 && i->target_vector < (long)directions)
                    {
                        SyncTupple<DT_> sync_tupple((*f_temp2[i->target_vector])[i->cell->get_id()], i->idx, i->target_vector);
                        sync_list.back().data.push_back(sync_tupple);
                    }
                    else
                        throw InternalError("Wrong target_vector in SyncInfo found!");
                }

                return sync_list;
            }

            Grid3<DT_, directions> grid;

            SharedArray<shared_ptr<DenseVector<unsigned long> > > neighbours;
            SharedArray<shared_ptr<DenseVector<unsigned long> > > dir;
            SharedArray<shared_ptr<DenseVector<unsigned long> > > dir_index;

            SharedArray<shared_ptr<DenseVector<DT_> > > f;
            SharedArray<shared_ptr<DenseVector<DT_> > > f_eq;
            SharedArray<shared_ptr<DenseVector<DT_> > > f_temp;
            SharedArray<shared_ptr<DenseVector<DT_> > > f_temp2;

            shared_ptr<DenseVector<DT_> > h;
            shared_ptr<DenseVector<DT_> > h2;
            shared_ptr<DenseVector<DT_> > b;
            shared_ptr<DenseVector<DT_> > u;
            shared_ptr<DenseVector<DT_> > v;

            shared_ptr<DenseVector<DT_> > distribution_x;
            shared_ptr<DenseVector<DT_> > distribution_y;


            PackedGrid3(Grid3<DT_, directions> & grid3) :
                grid(grid3),
                neighbours(directions),
                dir(directions),
                dir_index(directions),
                f(directions),
                f_eq(directions),
                f_temp(directions),
                f_temp2(directions)
            {
                // fill raw neighbour vectors
                for (unsigned long i(0) ; i < directions ; ++i)
                {
                    neighbours[i].reset(new DenseVector<unsigned long>(grid.size(), -(1ul)));
                    for (unsigned long idx(0) ; idx < grid.local_size() ; ++idx)
                    {
                        if (i == 0)
                            (*neighbours[i])[idx] = grid.get_cell(idx)->get_id();
                        else if (grid.get_cell(idx)->get_neighbours(i).size() != 0)
                            (*neighbours[i])[idx] = grid.get_cell(idx)->get_neighbours(i).front()->get_id();
                    }
                    for (unsigned long idx(grid.size() - grid.inner_halo_size()) ; idx < grid.size() ; ++idx)
                    {
                        if (i == 0)
                            (*neighbours[i])[idx] = grid.get_cell(idx)->get_id();
                        else if (grid.get_cell(idx)->get_neighbours(i).size() != 0)
                            (*neighbours[i])[idx] = grid.get_cell(idx)->get_neighbours(i).front()->get_id();
                    }
                }

                //fill packed direction vectors
                for (unsigned long i(0) ; i < directions ; ++i)
                {
                    unsigned long idx(0);
                    std::vector<unsigned long> temp_dir;
                    std::vector<unsigned long> temp_dir_index;
                    unsigned long old(-(1ul));
                    while (idx < grid.size())
                    {
                        if (temp_dir_index.size() % 2 == 0 && (*neighbours[i])[idx] != -(1ul))
                        {
                            temp_dir_index.push_back(idx);
                            temp_dir.push_back((*neighbours[i])[idx]);
                            old = (*neighbours[i])[idx];

                            if(idx == grid.size() - 1 ||  (*neighbours[i])[idx+1] != old+1)
                            {
                                ++idx;
                                temp_dir_index.push_back(idx);
                                continue;
                            }
                        }
                        else if (temp_dir_index.size() % 2 == 1 && (idx == grid.size() - 1 || (*neighbours[i])[idx+1] != old+1))
                        {
                            temp_dir_index.push_back(idx+1);
                        }
                        ++idx;
                        if (old != -(1ul))
                            ++old;
                    }

                    dir[i].reset(new DenseVector<unsigned long>(temp_dir.size()));
                    dir_index[i].reset(new DenseVector<unsigned long>(temp_dir_index.size()));
                    for (unsigned long j(0) ; j < temp_dir.size() ; ++j)
                    {
                        (*dir[i])[j] = temp_dir.at(j);
                        (*dir_index[i])[j*2] = temp_dir_index.at(j*2);
                        (*dir_index[i])[j*2 + 1] = temp_dir_index.at(j*2 + 1);
                    }
                }

                // fill h, b, u, v
                h.reset(new DenseVector<DT_>(grid.size()));
                h2.reset(new DenseVector<DT_>(grid.size()));
                b.reset(new DenseVector<DT_>(grid.size()));
                u.reset(new DenseVector<DT_>(grid.size()));
                v.reset(new DenseVector<DT_>(grid.size()));
                for (unsigned long idx(0) ; idx < grid.size() ; ++idx)
                {
                    (*h)[idx] = grid.get_cell(idx)->get_h();
                    (*h2)[idx] = grid.get_cell(idx)->get_h();
                    (*b)[idx] = grid.get_cell(idx)->get_b();
                    (*u)[idx] = grid.get_cell(idx)->get_u();
                    (*v)[idx] = grid.get_cell(idx)->get_v();
                }

                // fill f's
                for (unsigned long i(0) ; i < directions ; ++i)
                {
                    f[i].reset(new DenseVector<DT_>(grid.size(), 0));
                    f_eq[i].reset(new DenseVector<DT_>(grid.size(), 0));
                    f_temp[i].reset(new DenseVector<DT_>(grid.size(), 0));
                    f_temp2[i].reset(new DenseVector<DT_>(grid.size(), 0));
                }

                distribution_x.reset(new DenseVector<DT_>(directions));
                distribution_y.reset(new DenseVector<DT_>(directions));
            }
    };
}

#endif
