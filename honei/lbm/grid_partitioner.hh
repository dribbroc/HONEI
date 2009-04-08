/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the LBM C++ library. LBM is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LBM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */


#ifndef LBM_GUARD_GRID_PARTITIONER_HH
#define LBM_GUARD_GRID_PARTITIONER_HH 1

#include <honei/lbm/grid.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/dense_matrix.hh>

#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>

/**
 * \file
 * Definition of LBM grid partitioner.
 *
 * \ingroup grpliblbm
 **/

using namespace honei;
using namespace lbm;
using namespace lbm_lattice_types;

namespace honei
{
    template <typename LatticeType_, typename DT_> struct GridPartitioner
    {
    };

    template <typename DT_> struct GridPartitioner<D2Q9, DT_>
    {
        private:
            /// Scatter f_temp elements to other patches
            static void _synch_temp_1(unsigned long patch,
                    DenseVector<unsigned long> & index_vector, DenseVector<unsigned long> & targets,
                    std::vector<PackedGridInfo<D2Q9> > & info_list, std::vector<PackedGridData<D2Q9, DT_> > & data_list)
            {
                for (unsigned long i(0) ; i < index_vector.size() - 1 ; i += 2)
                {
                    unsigned long target(targets[i / 2]);
                    for (unsigned long j(index_vector[i]) ; j < index_vector[i + 1] ; ++j)
                    {
                        (*data_list[target].f_temp_1)[j - info_list[target].offset] = (*data_list[patch].f_temp_1)[j - info_list[patch].offset];
                    }
                }
            }

            static void _synch_temp_2(unsigned long patch,
                    DenseVector<unsigned long> & index_vector, DenseVector<unsigned long> & targets,
                    std::vector<PackedGridInfo<D2Q9> > & info_list, std::vector<PackedGridData<D2Q9, DT_> > & data_list)
            {
                for (unsigned long i(0) ; i < index_vector.size() - 1 ; i += 2)
                {
                    unsigned long target(targets[i / 2]);
                    for (unsigned long j(index_vector[i]) ; j < index_vector[i + 1] ; ++j)
                    {
                        (*data_list[target].f_temp_2)[j - info_list[target].offset] = (*data_list[patch].f_temp_2)[j - info_list[patch].offset];
                    }
                }
            }

            static void _synch_temp_3(unsigned long patch,
                    DenseVector<unsigned long> & index_vector, DenseVector<unsigned long> & targets,
                    std::vector<PackedGridInfo<D2Q9> > & info_list, std::vector<PackedGridData<D2Q9, DT_> > & data_list)
            {
                for (unsigned long i(0) ; i < index_vector.size() - 1 ; i += 2)
                {
                    unsigned long target(targets[i / 2]);
                    for (unsigned long j(index_vector[i]) ; j < index_vector[i + 1] ; ++j)
                    {
                        (*data_list[target].f_temp_3)[j - info_list[target].offset] = (*data_list[patch].f_temp_3)[j - info_list[patch].offset];
                    }
                }
            }

            static void _synch_temp_4(unsigned long patch,
                    DenseVector<unsigned long> & index_vector, DenseVector<unsigned long> & targets,
                    std::vector<PackedGridInfo<D2Q9> > & info_list, std::vector<PackedGridData<D2Q9, DT_> > & data_list)
            {
                for (unsigned long i(0) ; i < index_vector.size() - 1 ; i += 2)
                {
                    unsigned long target(targets[i / 2]);
                    for (unsigned long j(index_vector[i]) ; j < index_vector[i + 1] ; ++j)
                    {
                        (*data_list[target].f_temp_4)[j - info_list[target].offset] = (*data_list[patch].f_temp_4)[j - info_list[patch].offset];
                    }
                }
            }

            static void _synch_temp_5(unsigned long patch,
                    DenseVector<unsigned long> & index_vector, DenseVector<unsigned long> & targets,
                    std::vector<PackedGridInfo<D2Q9> > & info_list, std::vector<PackedGridData<D2Q9, DT_> > & data_list)
            {
                for (unsigned long i(0) ; i < index_vector.size() - 1 ; i += 2)
                {
                    unsigned long target(targets[i / 2]);
                    for (unsigned long j(index_vector[i]) ; j < index_vector[i + 1] ; ++j)
                    {
                        (*data_list[target].f_temp_5)[j - info_list[target].offset] = (*data_list[patch].f_temp_5)[j - info_list[patch].offset];
                    }
                }
            }

            static void _synch_temp_6(unsigned long patch,
                    DenseVector<unsigned long> & index_vector, DenseVector<unsigned long> & targets,
                    std::vector<PackedGridInfo<D2Q9> > & info_list, std::vector<PackedGridData<D2Q9, DT_> > & data_list)
            {
                for (unsigned long i(0) ; i < index_vector.size() - 1 ; i += 2)
                {
                    unsigned long target(targets[i / 2]);
                    for (unsigned long j(index_vector[i]) ; j < index_vector[i + 1] ; ++j)
                    {
                        (*data_list[target].f_temp_6)[j - info_list[target].offset] = (*data_list[patch].f_temp_6)[j - info_list[patch].offset];
                    }
                }
            }

            static void _synch_temp_7(unsigned long patch,
                    DenseVector<unsigned long> & index_vector, DenseVector<unsigned long> & targets,
                    std::vector<PackedGridInfo<D2Q9> > & info_list, std::vector<PackedGridData<D2Q9, DT_> > & data_list)
            {
                for (unsigned long i(0) ; i < index_vector.size() - 1 ; i += 2)
                {
                    unsigned long target(targets[i / 2]);
                    for (unsigned long j(index_vector[i]) ; j < index_vector[i + 1] ; ++j)
                    {
                        (*data_list[target].f_temp_7)[j - info_list[target].offset] = (*data_list[patch].f_temp_7)[j - info_list[patch].offset];
                    }
                }
            }

            static void _synch_temp_8(unsigned long patch,
                    DenseVector<unsigned long> & index_vector, DenseVector<unsigned long> & targets,
                    std::vector<PackedGridInfo<D2Q9> > & info_list, std::vector<PackedGridData<D2Q9, DT_> > & data_list)
            {
                for (unsigned long i(0) ; i < index_vector.size() - 1 ; i += 2)
                {
                    unsigned long target(targets[i / 2]);
                    for (unsigned long j(index_vector[i]) ; j < index_vector[i + 1] ; ++j)
                    {
                        (*data_list[target].f_temp_8)[j - info_list[target].offset] = (*data_list[patch].f_temp_8)[j - info_list[patch].offset];
                    }
                }
            }

            /// Gather h,u,v elements from other patches
            static void _synch_h(unsigned long patch,
                    DenseVector<unsigned long> & index_vector, DenseVector<unsigned long> & targets,
                    std::vector<PackedGridInfo<D2Q9> > & info_list, std::vector<PackedGridData<D2Q9, DT_> > & data_list)
            {
                for (unsigned long i(0) ; i < index_vector.size() - 1 ; i += 2)
                {
                    for (unsigned long j(index_vector[i]) ; j < index_vector[i + 1] ; ++j)
                    {
                        unsigned long target(targets[i / 2]);
                        (*data_list[patch].h)[j - info_list[patch].offset] = (*data_list[target].h)[j - info_list[target].offset];
                        //(*data_list[patch].u)[j - info_list[patch].offset] = (*data_list[target].u)[j - info_list[target].offset];
                        //(*data_list[patch].v)[j - info_list[patch].offset] = (*data_list[target].v)[j - info_list[target].offset];
                    }
                }
            }

            /// Store the element indices, that other patches need from us
            static void _create_dir_fringe(unsigned long patch, DenseVector<unsigned long> & dir_index,
                    DenseVector<unsigned long> & dir,
                    std::vector<PackedGridInfo<D2Q9> > & info_list,
                    DenseVector<unsigned long> * &new_dir_index, DenseVector<unsigned long> * &new_dir_targets)
            {
                std::vector<unsigned long> temp_dir_index;
                std::vector<unsigned long> temp_dir_targets;

                for (unsigned long index(0) ; index < dir.size() ; ++index)
                {
                    unsigned long start(dir_index[2 * index]);
                    unsigned long end(dir_index[2 * index + 1]);
                    for (unsigned long offset(0) ; offset < end - start ; ++offset)
                    {
                        // elements before our own data
                        if (dir[index] + offset < (*info_list[patch].limits)[0])
                        {
                            temp_dir_targets.push_back(patch - 1);
                            temp_dir_index.push_back(dir[index] + info_list[patch].offset + offset);
                            temp_dir_index.push_back(dir[index] + info_list[patch].offset + 1 + offset);
                        }
                        // elements behind our own data
                        if (dir[index] + offset >= (*info_list[patch].limits)[info_list[patch].limits->size() - 1])
                        {
                            temp_dir_targets.push_back(patch + 1);
                            temp_dir_index.push_back(dir[index] + info_list[patch].offset + offset);
                            temp_dir_index.push_back(dir[index] + info_list[patch].offset + 1 + offset);
                        }
                    }
                }

                if (temp_dir_targets.size() == 0)
                {
                    temp_dir_targets.push_back(0);
                    temp_dir_index.push_back(0);
                    temp_dir_index.push_back(0);
                }

                // postpack
                std::vector<unsigned long> packed_index;
                std::vector<unsigned long> packed_targets;
                unsigned long start(0);
                while (start < temp_dir_index.size())
                {
                    packed_index.push_back(temp_dir_index[start]);
                    packed_targets.push_back(temp_dir_targets[start / 2]);
                    while (start + 2 < temp_dir_index.size() && temp_dir_index[start + 2] == temp_dir_index[start] + 1)
                    {
                        start += 2;
                    }
                    packed_index.push_back(temp_dir_index[start + 1]);
                    start += 2;
                }

                // Fill the real vectors
                new_dir_index = new DenseVector<unsigned long>(packed_index.size());
                for (unsigned long i(0) ; i < packed_index.size() ; ++i)
                {
                    (*new_dir_index)[i] = packed_index[i];
                }
                new_dir_targets = new DenseVector<unsigned long>(packed_targets.size());
                for (unsigned long i(0) ; i < packed_targets.size() ; ++i)
                {
                    (*new_dir_targets)[i] = packed_targets[i];
                }
            }

            /// Store the element indices, that we need from other patches
            static void _create_h_fringe(unsigned long patch, DenseVector<unsigned long> & limits,
                    DenseVector<DT_> & h,
                    std::vector<PackedGridInfo<D2Q9> > & info_list,
                    DenseVector<unsigned long> * &new_h_index, DenseVector<unsigned long> * &new_h_targets)
            {
                std::vector<unsigned long> temp_h_index;
                std::vector<unsigned long> temp_h_targets;

                // elements before our own data
                for (unsigned long index(0) ; index < limits[0] ; ++index)
                {
                    temp_h_targets.push_back(patch - 1);
                    temp_h_index.push_back(index + info_list[patch].offset);
                    temp_h_index.push_back(index + info_list[patch].offset + 1);
                }
                // elements behind our own data
                for (unsigned long index(limits[limits.size() - 1]) ; index < h.size() ; ++index)
                {
                    temp_h_targets.push_back(patch + 1);
                    temp_h_index.push_back(index + info_list[patch].offset);
                    temp_h_index.push_back(index + info_list[patch].offset + 1);
                }
                if (temp_h_targets.size() == 0)
                {
                    temp_h_targets.push_back(0);
                    temp_h_index.push_back(0);
                    temp_h_index.push_back(0);
                }

                // postpack
                std::vector<unsigned long> packed_index;
                std::vector<unsigned long> packed_targets;
                unsigned long start(0);
                while (start < temp_h_index.size())
                {
                    packed_index.push_back(temp_h_index[start]);
                    packed_targets.push_back(temp_h_targets[start / 2]);
                    while (start + 2 < temp_h_index.size() && temp_h_index[start + 2] == temp_h_index[start] + 1)
                    {
                        start += 2;
                    }
                    packed_index.push_back(temp_h_index[start + 1]);
                    start += 2;
                }

                // Fill the real vectors
                new_h_index = new DenseVector<unsigned long>(packed_index.size());
                for (unsigned long i(0) ; i < packed_index.size() ; ++i)
                {
                    (*new_h_index)[i] = packed_index[i];
                }
                new_h_targets = new DenseVector<unsigned long>(packed_targets.size());
                for (unsigned long i(0) ; i < packed_targets.size() ; ++i)
                {
                    (*new_h_targets)[i] = packed_targets[i];
                }
            }

            // Gather all external modified parts of patch self
            static void _add_external_fringe(unsigned long self, std::vector<PackedGridFringe<D2Q9> > & fringe_list)
            {
                std::vector<unsigned long> temp_external_1;
                for (unsigned long index(0) ; index < fringe_list.size() ; ++index)
                {
                    for (unsigned long i(0) ; i < fringe_list[index].dir_targets_1->size() ; ++i)
                    {
                        if ((*fringe_list[index].dir_targets_1)[i] == self
                                && (
                                    //prevent target==0 index==[0,0] dummy entries
                                    (*fringe_list[index].dir_index_1)[i*2] != 0||
                                    (*fringe_list[index].dir_index_1)[i*2 + 1] != 0
                                    ))
                        {
                            temp_external_1.push_back((*fringe_list[index].dir_index_1)[i*2]);
                            temp_external_1.push_back((*fringe_list[index].dir_index_1)[i*2 + 1]);
                        }
                    }
                }

                if (temp_external_1.size() == 0)
                {
                    temp_external_1.push_back(0);
                    temp_external_1.push_back(0);
                }

                fringe_list[self].external_dir_index_1 = new DenseVector<unsigned long>(temp_external_1.size());
                for (unsigned long i(0) ; i < temp_external_1.size() ; ++i)
                {
                    (*fringe_list[self].external_dir_index_1)[i] = temp_external_1.at(i);
                }

                std::vector<unsigned long> temp_external_2;
                for (unsigned long index(0) ; index < fringe_list.size() ; ++index)
                {
                    for (unsigned long i(0) ; i < fringe_list[index].dir_targets_2->size() ; ++i)
                    {
                        if ((*fringe_list[index].dir_targets_2)[i] == self)
                        {
                            temp_external_2.push_back((*fringe_list[index].dir_index_2)[i*2]);
                            temp_external_2.push_back((*fringe_list[index].dir_index_2)[i*2 + 1]);
                        }
                    }
                }

                if (temp_external_2.size() == 0)
                {
                    temp_external_2.push_back(0);
                    temp_external_2.push_back(0);
                }

                fringe_list[self].external_dir_index_2 = new DenseVector<unsigned long>(temp_external_2.size());
                for (unsigned long i(0) ; i < temp_external_2.size() ; ++i)
                {
                    (*fringe_list[self].external_dir_index_2)[i] = temp_external_2.at(i);
                }

                std::vector<unsigned long> temp_external_3;
                for (unsigned long index(0) ; index < fringe_list.size() ; ++index)
                {
                    for (unsigned long i(0) ; i < fringe_list[index].dir_targets_3->size() ; ++i)
                    {
                        if ((*fringe_list[index].dir_targets_3)[i] == self)
                        {
                            temp_external_3.push_back((*fringe_list[index].dir_index_3)[i*2]);
                            temp_external_3.push_back((*fringe_list[index].dir_index_3)[i*2 + 1]);
                        }
                    }
                }

                if (temp_external_3.size() == 0)
                {
                    temp_external_3.push_back(0);
                    temp_external_3.push_back(0);
                }

                fringe_list[self].external_dir_index_3 = new DenseVector<unsigned long>(temp_external_3.size());
                for (unsigned long i(0) ; i < temp_external_3.size() ; ++i)
                {
                    (*fringe_list[self].external_dir_index_3)[i] = temp_external_3.at(i);
                }

                std::vector<unsigned long> temp_external_4;
                for (unsigned long index(0) ; index < fringe_list.size() ; ++index)
                {
                    for (unsigned long i(0) ; i < fringe_list[index].dir_targets_4->size() ; ++i)
                    {
                        if ((*fringe_list[index].dir_targets_4)[i] == self)
                        {
                            temp_external_4.push_back((*fringe_list[index].dir_index_4)[i*2]);
                            temp_external_4.push_back((*fringe_list[index].dir_index_4)[i*2 + 1]);
                        }
                    }
                }

                if (temp_external_4.size() == 0)
                {
                    temp_external_4.push_back(0);
                    temp_external_4.push_back(0);
                }

                fringe_list[self].external_dir_index_4 = new DenseVector<unsigned long>(temp_external_4.size());
                for (unsigned long i(0) ; i < temp_external_4.size() ; ++i)
                {
                    (*fringe_list[self].external_dir_index_4)[i] = temp_external_4.at(i);
                }

                std::vector<unsigned long> temp_external_5;
                for (unsigned long index(0) ; index < fringe_list.size() ; ++index)
                {
                    for (unsigned long i(0) ; i < fringe_list[index].dir_targets_5->size() ; ++i)
                    {
                        if ((*fringe_list[index].dir_targets_5)[i] == self)
                        {
                            temp_external_5.push_back((*fringe_list[index].dir_index_5)[i*2]);
                            temp_external_5.push_back((*fringe_list[index].dir_index_5)[i*2 + 1]);
                        }
                    }
                }

                if (temp_external_5.size() == 0)
                {
                    temp_external_5.push_back(0);
                    temp_external_5.push_back(0);
                }

                fringe_list[self].external_dir_index_5 = new DenseVector<unsigned long>(temp_external_5.size());
                for (unsigned long i(0) ; i < temp_external_5.size() ; ++i)
                {
                    (*fringe_list[self].external_dir_index_5)[i] = temp_external_5.at(i);
                }

                std::vector<unsigned long> temp_external_6;
                for (unsigned long index(0) ; index < fringe_list.size() ; ++index)
                {
                    for (unsigned long i(0) ; i < fringe_list[index].dir_targets_6->size() ; ++i)
                    {
                        if ((*fringe_list[index].dir_targets_6)[i] == self)
                        {
                            temp_external_6.push_back((*fringe_list[index].dir_index_6)[i*2]);
                            temp_external_6.push_back((*fringe_list[index].dir_index_6)[i*2 + 1]);
                        }
                    }
                }

                if (temp_external_6.size() == 0)
                {
                    temp_external_6.push_back(0);
                    temp_external_6.push_back(0);
                }

                fringe_list[self].external_dir_index_6 = new DenseVector<unsigned long>(temp_external_6.size());
                for (unsigned long i(0) ; i < temp_external_6.size() ; ++i)
                {
                    (*fringe_list[self].external_dir_index_6)[i] = temp_external_6.at(i);
                }

                std::vector<unsigned long> temp_external_7;
                for (unsigned long index(0) ; index < fringe_list.size() ; ++index)
                {
                    for (unsigned long i(0) ; i < fringe_list[index].dir_targets_7->size() ; ++i)
                    {
                        if ((*fringe_list[index].dir_targets_7)[i] == self)
                        {
                            temp_external_7.push_back((*fringe_list[index].dir_index_7)[i*2]);
                            temp_external_7.push_back((*fringe_list[index].dir_index_7)[i*2 + 1]);
                        }
                    }
                }

                if (temp_external_7.size() == 0)
                {
                    temp_external_7.push_back(0);
                    temp_external_7.push_back(0);
                }

                fringe_list[self].external_dir_index_7 = new DenseVector<unsigned long>(temp_external_7.size());
                for (unsigned long i(0) ; i < temp_external_7.size() ; ++i)
                {
                    (*fringe_list[self].external_dir_index_7)[i] = temp_external_7.at(i);
                }

                std::vector<unsigned long> temp_external_8;
                for (unsigned long index(0) ; index < fringe_list.size() ; ++index)
                {
                    for (unsigned long i(0) ; i < fringe_list[index].dir_targets_8->size() ; ++i)
                    {
                        if ((*fringe_list[index].dir_targets_8)[i] == self)
                        {
                            temp_external_8.push_back((*fringe_list[index].dir_index_8)[i*2]);
                            temp_external_8.push_back((*fringe_list[index].dir_index_8)[i*2 + 1]);
                        }
                    }
                }

                if (temp_external_8.size() == 0)
                {
                    temp_external_8.push_back(0);
                    temp_external_8.push_back(0);
                }

                fringe_list[self].external_dir_index_8 = new DenseVector<unsigned long>(temp_external_8.size());
                for (unsigned long i(0) ; i < temp_external_8.size() ; ++i)
                {
                    (*fringe_list[self].external_dir_index_8)[i] = temp_external_8.at(i);
                }
            }

        public:

            static void compose(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data,
                    std::vector<PackedGridInfo<D2Q9> > & info_list, std::vector<PackedGridData<D2Q9, DT_> > & data_list)
            {
                for(unsigned long i(0); i < data.h->size(); ++i)
                {
                    unsigned long index(0);
                    while(info_list[index].offset + (*info_list[index].limits)[info_list[index].limits->size() - 1] <= i)
                    {
                        ++index;
                    }
                    (*data.h)[i] = (*data_list[index].h)[i - (info_list[index].offset)];
                }
            }

            static void synch(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data,
                    std::vector<PackedGridInfo<D2Q9> > & info_list, std::vector<PackedGridData<D2Q9, DT_> > & data_list,
                    std::vector<PackedGridFringe<D2Q9> > & fringe_list)
            {
                for (unsigned long index(0) ; index < fringe_list.size() ; ++index)
                {
                    _synch_temp_1(index, *(fringe_list[index].dir_index_1), *(fringe_list[index].dir_targets_1), info_list, data_list);
                    _synch_temp_2(index, *(fringe_list[index].dir_index_2), *(fringe_list[index].dir_targets_2), info_list, data_list);
                    _synch_temp_3(index, *(fringe_list[index].dir_index_3), *(fringe_list[index].dir_targets_3), info_list, data_list);
                    _synch_temp_4(index, *(fringe_list[index].dir_index_4), *(fringe_list[index].dir_targets_4), info_list, data_list);
                    _synch_temp_5(index, *(fringe_list[index].dir_index_5), *(fringe_list[index].dir_targets_5), info_list, data_list);
                    _synch_temp_6(index, *(fringe_list[index].dir_index_6), *(fringe_list[index].dir_targets_6), info_list, data_list);
                    _synch_temp_7(index, *(fringe_list[index].dir_index_7), *(fringe_list[index].dir_targets_7), info_list, data_list);
                    _synch_temp_8(index, *(fringe_list[index].dir_index_8), *(fringe_list[index].dir_targets_8), info_list, data_list);
                    _synch_h(index, *(fringe_list[index].h_index), *(fringe_list[index].h_targets), info_list, data_list);
                }
            }

            static void partition_directions(std::vector<unsigned long> & dir_index, std::vector<unsigned long> & dir,
                    std::vector<unsigned long> & barriers, unsigned long start, unsigned long end)
            {
                CONTEXT("When partitioning one single direction:");
                if (dir_index.size() == 0 || dir_index[dir_index.size() - 1] < start || dir_index[0] > end)
                {
                    barriers.push_back(0);
                    barriers.push_back(0);
                    return;
                }

                std::vector<unsigned long>::iterator start_it(dir_index.begin());
                unsigned long start_index(0);
                while (start_index < dir_index.size() - 2 && dir_index[start_index] < start && dir_index[start_index + 1] <= start)
                {
                    start_index += 2;
                    start_it += 2;
                }

                if (dir_index[start_index + 1] < start)
                {
                    start_index += 2;
                    start_it += 2;
                }

                std::vector<unsigned long>::iterator end_it(start_it);
                ++end_it;
                unsigned long end_index(start_index);
                ++end_index;
                while (end_index < dir_index.size() - 1 && dir_index[end_index] < end)
                {
                    end_index += 2;
                    end_it += 2;
                }
                if (dir_index[end_index - 1] > end)
                {
                    end_index -= 2;
                    end_it -= 2;
                }
                //std::cout<<"start_index: "<<start_index<<"("<<*start_it<<") end_index: "<<end_index<<"("<<*end_it<<")"<<std::endl;
                if (dir_index[start_index] < start)
                {
                    std::vector<unsigned long>::iterator temp_it(start_it);
                    ++temp_it;
                    temp_it = dir_index.insert(temp_it, start);
                    dir_index.insert(temp_it, start);

                    std::vector<unsigned long>::iterator temp_dir(dir.begin());
                    for (unsigned long i(0) ; i < start_index / 2 ; ++i)
                    {
                        ++temp_dir;
                    }
                    temp_dir = dir.insert(temp_dir, *temp_dir);
                    ++temp_dir;
                    *temp_dir = *temp_dir + (dir_index[start_index + 2] - dir_index[start_index]);
                    barriers.push_back(start_index + 2);

                    end_index += 2;
                }
                else
                {
                    //nur start-barrier einfuegen
                    barriers.push_back(start_index);
                }

                end_it = dir_index.begin();
                for (unsigned long i(0) ; i < end_index ; ++i)
                {
                    ++end_it;
                }
                if (dir_index[end_index] >= end)
                {
                    std::vector<unsigned long>::iterator temp_it(end_it);
                    temp_it = dir_index.insert(temp_it, end);
                    dir_index.insert(temp_it, end);

                    std::vector<unsigned long>::iterator temp_dir(dir.begin());
                    for (unsigned long i(0) ; i < end_index / 2 ; ++i)
                    {
                        ++temp_dir;
                    }
                    temp_dir = dir.insert(temp_dir, *temp_dir);
                    ++temp_dir;
                    *temp_dir = *temp_dir + (dir_index[end_index + 1] - dir_index[end_index - 1]);
                    barriers.push_back(end_index);
                }
                else
                {
                    //nur end-barrier einfuegen
                    barriers.push_back(end_index);
                }
            }

            static void decompose(unsigned long parts, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data,
                    std::vector<PackedGridInfo<D2Q9> > & info_list, std::vector<PackedGridData<D2Q9, DT_> > & data_list,
                    std::vector<PackedGridFringe<D2Q9> > & fringe_list)
            {
                CONTEXT("When creating grid partitions:");
                std::vector<unsigned long> temp_limits;
                std::vector<unsigned long> temp_types;
                std::vector<unsigned long> temp_dir_1;
                std::vector<unsigned long> temp_dir_2;
                std::vector<unsigned long> temp_dir_3;
                std::vector<unsigned long> temp_dir_4;
                std::vector<unsigned long> temp_dir_5;
                std::vector<unsigned long> temp_dir_6;
                std::vector<unsigned long> temp_dir_7;
                std::vector<unsigned long> temp_dir_8;
                std::vector<unsigned long> temp_dir_index_1;
                std::vector<unsigned long> temp_dir_index_2;
                std::vector<unsigned long> temp_dir_index_3;
                std::vector<unsigned long> temp_dir_index_4;
                std::vector<unsigned long> temp_dir_index_5;
                std::vector<unsigned long> temp_dir_index_6;
                std::vector<unsigned long> temp_dir_index_7;
                std::vector<unsigned long> temp_dir_index_8;
                std::vector<unsigned long> barriers;
                std::vector<unsigned long> barriers_1;
                std::vector<unsigned long> barriers_2;
                std::vector<unsigned long> barriers_3;
                std::vector<unsigned long> barriers_4;
                std::vector<unsigned long> barriers_5;
                std::vector<unsigned long> barriers_6;
                std::vector<unsigned long> barriers_7;
                std::vector<unsigned long> barriers_8;

                for (unsigned long i(0) ; i < info.limits->size() ; ++i)
                {
                    temp_limits.push_back((*info.limits)[i]);
                    temp_types.push_back((*info.types)[i]);
                }
                for (unsigned long i(0) ; i < info.dir_1->size() ; ++i)
                {
                    temp_dir_1.push_back((*info.dir_1)[i]);
                    temp_dir_index_1.push_back((*info.dir_index_1)[i * 2]);
                    temp_dir_index_1.push_back((*info.dir_index_1)[i * 2 + 1]);
                }
                for (unsigned long i(0) ; i < info.dir_2->size() ; ++i)
                {
                    temp_dir_2.push_back((*info.dir_2)[i]);
                    temp_dir_index_2.push_back((*info.dir_index_2)[i * 2]);
                    temp_dir_index_2.push_back((*info.dir_index_2)[i * 2 + 1]);
                }
                for (unsigned long i(0) ; i < info.dir_3->size() ; ++i)
                {
                    temp_dir_3.push_back((*info.dir_3)[i]);
                    temp_dir_index_3.push_back((*info.dir_index_3)[i * 2]);
                    temp_dir_index_3.push_back((*info.dir_index_3)[i * 2 + 1]);
                }
                for (unsigned long i(0) ; i < info.dir_4->size() ; ++i)
                {
                    temp_dir_4.push_back((*info.dir_4)[i]);
                    temp_dir_index_4.push_back((*info.dir_index_4)[i * 2]);
                    temp_dir_index_4.push_back((*info.dir_index_4)[i * 2 + 1]);
                }
                for (unsigned long i(0) ; i < info.dir_5->size() ; ++i)
                {
                    temp_dir_5.push_back((*info.dir_5)[i]);
                    temp_dir_index_5.push_back((*info.dir_index_5)[i * 2]);
                    temp_dir_index_5.push_back((*info.dir_index_5)[i * 2 + 1]);
                }
                for (unsigned long i(0) ; i < info.dir_6->size() ; ++i)
                {
                    temp_dir_6.push_back((*info.dir_6)[i]);
                    temp_dir_index_6.push_back((*info.dir_index_6)[i * 2]);
                    temp_dir_index_6.push_back((*info.dir_index_6)[i * 2 + 1]);
                }
                for (unsigned long i(0) ; i < info.dir_7->size() ; ++i)
                {
                    temp_dir_7.push_back((*info.dir_7)[i]);
                    temp_dir_index_7.push_back((*info.dir_index_7)[i * 2]);
                    temp_dir_index_7.push_back((*info.dir_index_7)[i * 2 + 1]);
                }
                for (unsigned long i(0) ; i < info.dir_8->size() ; ++i)
                {
                    temp_dir_8.push_back((*info.dir_8)[i]);
                    temp_dir_index_8.push_back((*info.dir_index_8)[i * 2]);
                    temp_dir_index_8.push_back((*info.dir_index_8)[i * 2 + 1]);
                }

                if (parts == 0)
                    throw InternalError("GridPartitioner: Cannot decompose into 0 parts!");

                if (data.u->size() < parts)
                    parts = data.u->size();

                unsigned long normal_size((data.u->size() / parts));
                unsigned long first_size((data.u->size() / parts) + (data.u->size() % parts));

                unsigned long end(0);
                unsigned long start(0);
                barriers.push_back(0);
                for (unsigned long i(0) ; i < parts ; ++i)
                {
                    end += (i==0 ? first_size : normal_size);
                    std::vector<unsigned long>::iterator start_index(std::lower_bound(temp_limits.begin(), temp_limits.end(), start));
                    std::vector<unsigned long>::iterator end_index(std::lower_bound(temp_limits.begin(), temp_limits.end(), end));
                    //std::cout<<"start: "<<start<<" "<<*start_index<<std::endl;
                    //std::cout<<"end: "<<end<<" "<<*end_index<<std::endl;
                    if (*end_index != end)
                    {
                        // Insert dummy nodes
                        temp_limits.insert(end_index, end);

                        unsigned long i(0);
                        std::vector<unsigned long>::iterator it_types(temp_types.begin());
                        while (i != temp_limits.size() - 1 && temp_limits[i] != end)
                        {
                            ++i;
                            ++it_types;
                        }
                        barriers.push_back(i);
                        temp_types.insert(it_types, *(--it_types));
                    }
                    else
                    {
                        unsigned long i(0);
                        while (i != temp_limits.size() - 1 && temp_limits[i] != end)
                        {
                            ++i;
                        }
                        barriers.push_back(i);
                    }
                    // prepare direction vectors
                    partition_directions(temp_dir_index_1, temp_dir_1, barriers_1, start, end);
                    partition_directions(temp_dir_index_2, temp_dir_2, barriers_2, start, end);
                    partition_directions(temp_dir_index_3, temp_dir_3, barriers_3, start, end);
                    partition_directions(temp_dir_index_4, temp_dir_4, barriers_4, start, end);
                    partition_directions(temp_dir_index_5, temp_dir_5, barriers_5, start, end);
                    partition_directions(temp_dir_index_6, temp_dir_6, barriers_6, start, end);
                    partition_directions(temp_dir_index_7, temp_dir_7, barriers_7, start, end);
                    partition_directions(temp_dir_index_8, temp_dir_8, barriers_8, start, end);

                    start = end;
                }

                // Create partitions
                std::vector<std::vector<unsigned long> > temp_limits_list;
                std::vector<std::vector<unsigned long> > temp_types_list;
                std::vector<std::vector<unsigned long> > temp_dir_1_list;
                std::vector<std::vector<unsigned long> > temp_dir_2_list;
                std::vector<std::vector<unsigned long> > temp_dir_3_list;
                std::vector<std::vector<unsigned long> > temp_dir_4_list;
                std::vector<std::vector<unsigned long> > temp_dir_5_list;
                std::vector<std::vector<unsigned long> > temp_dir_6_list;
                std::vector<std::vector<unsigned long> > temp_dir_7_list;
                std::vector<std::vector<unsigned long> > temp_dir_8_list;
                std::vector<std::vector<unsigned long> > temp_dir_index_1_list;
                std::vector<std::vector<unsigned long> > temp_dir_index_2_list;
                std::vector<std::vector<unsigned long> > temp_dir_index_3_list;
                std::vector<std::vector<unsigned long> > temp_dir_index_4_list;
                std::vector<std::vector<unsigned long> > temp_dir_index_5_list;
                std::vector<std::vector<unsigned long> > temp_dir_index_6_list;
                std::vector<std::vector<unsigned long> > temp_dir_index_7_list;
                std::vector<std::vector<unsigned long> > temp_dir_index_8_list;
                std::vector<unsigned long> temp_min_list;
                std::vector<unsigned long> temp_max_list;

                for (unsigned long i(0) ; i < barriers.size() - 1; ++i)
                {
                    std::vector<unsigned long> new_limits;
                    std::vector<unsigned long> new_types;
                    std::vector<unsigned long> new_dir_1;
                    std::vector<unsigned long> new_dir_2;
                    std::vector<unsigned long> new_dir_3;
                    std::vector<unsigned long> new_dir_4;
                    std::vector<unsigned long> new_dir_5;
                    std::vector<unsigned long> new_dir_6;
                    std::vector<unsigned long> new_dir_7;
                    std::vector<unsigned long> new_dir_8;
                    std::vector<unsigned long> new_dir_index_1;
                    std::vector<unsigned long> new_dir_index_2;
                    std::vector<unsigned long> new_dir_index_3;
                    std::vector<unsigned long> new_dir_index_4;
                    std::vector<unsigned long> new_dir_index_5;
                    std::vector<unsigned long> new_dir_index_6;
                    std::vector<unsigned long> new_dir_index_7;
                    std::vector<unsigned long> new_dir_index_8;

                    std::vector<unsigned long> max_collection;
                    std::vector<unsigned long> min_collection;

                    for(unsigned long index(barriers[i]) ; index <= barriers[i + 1] ; ++index)
                    {
                        new_limits.push_back(temp_limits[index]);
                        new_types.push_back(temp_types[index]);
                    }
                    if (barriers_1[i * 2] != barriers_1[i * 2 + 1])
                    {
                        for (unsigned long index(barriers_1[i * 2]) ; index <= barriers_1[i * 2 + 1] ; ++index)
                        {
                            new_dir_index_1.push_back(temp_dir_index_1[index]);
                            if (index % 2 == 0)
                                new_dir_1.push_back(temp_dir_1[index / 2]);
                        }
                        max_collection.push_back(*max_element(new_dir_1.begin(), new_dir_1.end()) + new_dir_index_1[new_dir_index_1.size() - 1] - new_dir_index_1[new_dir_index_1.size() - 2] - 1);
                        min_collection.push_back(*min_element(new_dir_1.begin(), new_dir_1.end()));
                    }
                    else
                    {
                        new_dir_index_1.push_back(0);
                        new_dir_index_1.push_back(0);
                        new_dir_1.push_back(0);
                    }

                    if (barriers_2[i * 2] != barriers_2[i * 2 + 1])
                    {
                        for (unsigned long index(barriers_2[i * 2]) ; index <= barriers_2[i * 2 + 1] ; ++index)
                        {
                            new_dir_index_2.push_back(temp_dir_index_2[index]);
                            if (index % 2 == 0)
                                new_dir_2.push_back(temp_dir_2[index / 2]);
                        }
                        max_collection.push_back(*max_element(new_dir_2.begin(), new_dir_2.end()));
                        min_collection.push_back(*min_element(new_dir_2.begin(), new_dir_2.end()));
                    }
                    else
                    {
                        new_dir_index_2.push_back(0);
                        new_dir_index_2.push_back(0);
                        new_dir_2.push_back(0);
                    }

                    if (barriers_3[i * 2] != barriers_3[i * 2 + 1])
                    {
                        for (unsigned long index(barriers_3[i * 2]) ; index <= barriers_3[i * 2 + 1] ; ++index)
                        {
                            new_dir_index_3.push_back(temp_dir_index_3[index]);
                            if (index % 2 == 0)
                                new_dir_3.push_back(temp_dir_3[index / 2]);
                        }
                        max_collection.push_back(*max_element(new_dir_3.begin(), new_dir_3.end()));
                        min_collection.push_back(*min_element(new_dir_3.begin(), new_dir_3.end()));
                    }
                    else
                    {
                        new_dir_index_3.push_back(0);
                        new_dir_index_3.push_back(0);
                        new_dir_3.push_back(0);
                    }

                    if (barriers_4[i * 2] != barriers_4[i * 2 + 1])
                    {
                        for (unsigned long index(barriers_4[i * 2]) ; index <= barriers_4[i * 2 + 1] ; ++index)
                        {
                            new_dir_index_4.push_back(temp_dir_index_4[index]);
                            if (index % 2 == 0)
                                new_dir_4.push_back(temp_dir_4[index / 2]);
                        }
                        max_collection.push_back(*max_element(new_dir_4.begin(), new_dir_4.end()));
                        min_collection.push_back(*min_element(new_dir_4.begin(), new_dir_4.end()));
                    }
                    else
                    {
                        new_dir_index_4.push_back(0);
                        new_dir_index_4.push_back(0);
                        new_dir_4.push_back(0);
                    }

                    if (barriers_5[i * 2] != barriers_5[i * 2 + 1])
                    {
                        for (unsigned long index(barriers_5[i * 2]) ; index <= barriers_5[i * 2 + 1] ; ++index)
                        {
                            new_dir_index_5.push_back(temp_dir_index_5[index]);
                            if (index % 2 == 0)
                                new_dir_5.push_back(temp_dir_5[index / 2]);
                        }
                        max_collection.push_back(*max_element(new_dir_5.begin(), new_dir_5.end()));
                        min_collection.push_back(*min_element(new_dir_5.begin(), new_dir_5.end()));
                    }
                    else
                    {
                        new_dir_index_5.push_back(0);
                        new_dir_index_5.push_back(0);
                        new_dir_5.push_back(0);
                    }

                    if (barriers_6[i * 2] != barriers_6[i * 2 + 1])
                    {
                        for (unsigned long index(barriers_6[i * 2]) ; index <= barriers_6[i * 2 + 1] ; ++index)
                        {
                            new_dir_index_6.push_back(temp_dir_index_6[index]);
                            if (index % 2 == 0)
                                new_dir_6.push_back(temp_dir_6[index / 2]);
                        }
                        max_collection.push_back(*max_element(new_dir_6.begin(), new_dir_6.end()) + new_dir_index_6[new_dir_index_6.size() - 1] - new_dir_index_6[new_dir_index_6.size() - 2] - 1);
                        min_collection.push_back(*min_element(new_dir_6.begin(), new_dir_6.end()));
                    }
                    else
                    {
                        new_dir_index_6.push_back(0);
                        new_dir_index_6.push_back(0);
                        new_dir_6.push_back(0);
                    }

                    if (barriers_7[i * 2] != barriers_7[i * 2 + 1])
                    {
                        for (unsigned long index(barriers_7[i * 2]) ; index <= barriers_7[i * 2 + 1] ; ++index)
                        {
                            new_dir_index_7.push_back(temp_dir_index_7[index]);
                            if (index % 2 == 0)
                                new_dir_7.push_back(temp_dir_7[index / 2]);
                        }
                        max_collection.push_back(*max_element(new_dir_7.begin(), new_dir_7.end()) + new_dir_index_7[new_dir_index_7.size() - 1] - new_dir_index_7[new_dir_index_7.size() - 2] - 1);
                        min_collection.push_back(*min_element(new_dir_7.begin(), new_dir_7.end()));
                    }
                    else
                    {
                        new_dir_index_7.push_back(0);
                        new_dir_index_7.push_back(0);
                        new_dir_7.push_back(0);
                    }

                    if (barriers_8[i * 2] != barriers_8[i * 2 + 1])
                    {
                        for (unsigned long index(barriers_8[i * 2]) ; index <= barriers_8[i * 2 + 1] ; ++index)
                        {
                            new_dir_index_8.push_back(temp_dir_index_8[index]);
                            if (index % 2 == 0)
                                new_dir_8.push_back(temp_dir_8[index / 2]);
                        }
                        max_collection.push_back(*max_element(new_dir_8.begin(), new_dir_8.end()) + new_dir_index_8[new_dir_index_8.size() - 1] - new_dir_index_8[new_dir_index_8.size() - 2] - 1);
                        min_collection.push_back(*min_element(new_dir_7.begin(), new_dir_7.end()));
                    }
                    else
                    {
                        new_dir_index_8.push_back(0);
                        new_dir_index_8.push_back(0);
                        new_dir_8.push_back(0);
                    }

                    temp_limits_list.push_back(new_limits);
                    temp_types_list.push_back(new_types);
                    temp_dir_1_list.push_back(new_dir_1);
                    temp_dir_2_list.push_back(new_dir_2);
                    temp_dir_3_list.push_back(new_dir_3);
                    temp_dir_4_list.push_back(new_dir_4);
                    temp_dir_5_list.push_back(new_dir_5);
                    temp_dir_6_list.push_back(new_dir_6);
                    temp_dir_7_list.push_back(new_dir_7);
                    temp_dir_8_list.push_back(new_dir_8);
                    temp_dir_index_1_list.push_back(new_dir_index_1);
                    temp_dir_index_2_list.push_back(new_dir_index_2);
                    temp_dir_index_3_list.push_back(new_dir_index_3);
                    temp_dir_index_4_list.push_back(new_dir_index_4);
                    temp_dir_index_5_list.push_back(new_dir_index_5);
                    temp_dir_index_6_list.push_back(new_dir_index_6);
                    temp_dir_index_7_list.push_back(new_dir_index_7);
                    temp_dir_index_8_list.push_back(new_dir_index_8);

                    temp_max_list.push_back(*max_element(max_collection.begin(), max_collection.end()));

                    temp_min_list.push_back(*min_element(min_collection.begin(), min_collection.end()));
                }

                // Fill the real info vectors
                for (unsigned long i(0) ; i < temp_limits_list.size() ; ++i)
                {
                    PackedGridInfo<D2Q9> new_info;
                    new_info.offset = temp_min_list[i];
                    new_info.limits = new DenseVector<unsigned long>(temp_limits_list[i].size());
                    new_info.types = new DenseVector<unsigned long>(temp_limits_list[i].size());
                    new_info.dir_1 = new DenseVector<unsigned long>(temp_dir_1_list[i].size());
                    new_info.dir_2 = new DenseVector<unsigned long>(temp_dir_2_list[i].size());
                    new_info.dir_3 = new DenseVector<unsigned long>(temp_dir_3_list[i].size());
                    new_info.dir_4 = new DenseVector<unsigned long>(temp_dir_4_list[i].size());
                    new_info.dir_5 = new DenseVector<unsigned long>(temp_dir_5_list[i].size());
                    new_info.dir_6 = new DenseVector<unsigned long>(temp_dir_6_list[i].size());
                    new_info.dir_7 = new DenseVector<unsigned long>(temp_dir_7_list[i].size());
                    new_info.dir_8 = new DenseVector<unsigned long>(temp_dir_8_list[i].size());
                    new_info.dir_index_1 = new DenseVector<unsigned long>(temp_dir_index_1_list[i].size());
                    new_info.dir_index_2 = new DenseVector<unsigned long>(temp_dir_index_2_list[i].size());
                    new_info.dir_index_3 = new DenseVector<unsigned long>(temp_dir_index_3_list[i].size());
                    new_info.dir_index_4 = new DenseVector<unsigned long>(temp_dir_index_4_list[i].size());
                    new_info.dir_index_5 = new DenseVector<unsigned long>(temp_dir_index_5_list[i].size());
                    new_info.dir_index_6 = new DenseVector<unsigned long>(temp_dir_index_6_list[i].size());
                    new_info.dir_index_7 = new DenseVector<unsigned long>(temp_dir_index_7_list[i].size());
                    new_info.dir_index_8 = new DenseVector<unsigned long>(temp_dir_index_8_list[i].size());
                    for (unsigned long j(0) ; j < temp_limits_list[i].size() ; ++j)
                    {
                        (*new_info.limits)[j] = temp_limits_list[i][j] - new_info.offset;
                        (*new_info.types)[j] = temp_types_list[i][j];
                    }
                    for (unsigned long j(0) ; j < temp_dir_1_list[i].size() ; ++j)
                    {
                        (*new_info.dir_1)[j] = temp_dir_1_list[i][j] - new_info.offset;
                        (*new_info.dir_index_1)[j * 2] = temp_dir_index_1_list[i][j * 2] - new_info.offset;
                        (*new_info.dir_index_1)[j * 2 + 1] = temp_dir_index_1_list[i][j * 2 + 1] - new_info.offset;
                    }
                    for (unsigned long j(0) ; j < temp_dir_2_list[i].size() ; ++j)
                    {
                        (*new_info.dir_2)[j] = temp_dir_2_list[i][j] - new_info.offset;
                        (*new_info.dir_index_2)[j * 2] = temp_dir_index_2_list[i][j * 2] - new_info.offset;
                        (*new_info.dir_index_2)[j * 2 + 1] = temp_dir_index_2_list[i][j * 2 + 1] - new_info.offset;
                    }
                    for (unsigned long j(0) ; j < temp_dir_3_list[i].size() ; ++j)
                    {
                        (*new_info.dir_3)[j] = temp_dir_3_list[i][j] - new_info.offset;
                        (*new_info.dir_index_3)[j * 2] = temp_dir_index_3_list[i][j * 2] - new_info.offset;
                        (*new_info.dir_index_3)[j * 2 + 1] = temp_dir_index_3_list[i][j * 2 + 1] - new_info.offset;
                    }
                    for (unsigned long j(0) ; j < temp_dir_4_list[i].size() ; ++j)
                    {
                        (*new_info.dir_4)[j] = temp_dir_4_list[i][j] - new_info.offset;
                        (*new_info.dir_index_4)[j * 2] = temp_dir_index_4_list[i][j * 2] - new_info.offset;
                        (*new_info.dir_index_4)[j * 2 + 1] = temp_dir_index_4_list[i][j * 2 + 1] - new_info.offset;
                    }
                    for (unsigned long j(0) ; j < temp_dir_5_list[i].size() ; ++j)
                    {
                        (*new_info.dir_5)[j] = temp_dir_5_list[i][j] - new_info.offset;
                        (*new_info.dir_index_5)[j * 2] = temp_dir_index_5_list[i][j * 2] - new_info.offset;
                        (*new_info.dir_index_5)[j * 2 + 1] = temp_dir_index_5_list[i][j * 2 + 1] - new_info.offset;
                    }
                    for (unsigned long j(0) ; j < temp_dir_6_list[i].size() ; ++j)
                    {
                        (*new_info.dir_6)[j] = temp_dir_6_list[i][j] - new_info.offset;
                        (*new_info.dir_index_6)[j * 2] = temp_dir_index_6_list[i][j * 2] - new_info.offset;
                        (*new_info.dir_index_6)[j * 2 + 1] = temp_dir_index_6_list[i][j * 2 + 1] - new_info.offset;
                    }
                    for (unsigned long j(0) ; j < temp_dir_7_list[i].size() ; ++j)
                    {
                        (*new_info.dir_7)[j] = temp_dir_7_list[i][j] - new_info.offset;
                        (*new_info.dir_index_7)[j * 2] = temp_dir_index_7_list[i][j * 2] - new_info.offset;
                        (*new_info.dir_index_7)[j * 2 + 1] = temp_dir_index_7_list[i][j * 2 + 1] - new_info.offset;
                    }
                    for (unsigned long j(0) ; j < temp_dir_8_list[i].size() ; ++j)
                    {
                        (*new_info.dir_8)[j] = temp_dir_8_list[i][j] - new_info.offset;
                        (*new_info.dir_index_8)[j * 2] = temp_dir_index_8_list[i][j * 2] - new_info.offset;
                        (*new_info.dir_index_8)[j * 2 + 1] = temp_dir_index_8_list[i][j * 2 + 1] - new_info.offset;
                    }
                    info_list.push_back(new_info);
                }

                // Fill the real data vectors
                for (unsigned long i(0) ; i < info_list.size() ; ++i)
                {
                    PackedGridData<D2Q9, DT_> new_data;
                    unsigned long data_size(temp_max_list[i] - temp_min_list[i]);
                    new_data.h = new DenseVector<DT_>(data_size + 1);
                    new_data.u = new DenseVector<DT_>(data_size + 1);
                    new_data.v = new DenseVector<DT_>(data_size + 1);
                    new_data.b = new DenseVector<DT_>(data_size + 1);
                    new_data.temp = new DenseVector<DT_>(data_size + 1);

                    new_data.f_0 = new DenseVector<DT_>(data_size + 1, DT_(0));
                    new_data.f_1 = new DenseVector<DT_>(data_size + 1, DT_(0));
                    new_data.f_2 = new DenseVector<DT_>(data_size + 1, DT_(0));
                    new_data.f_3 = new DenseVector<DT_>(data_size + 1, DT_(0));
                    new_data.f_4 = new DenseVector<DT_>(data_size + 1, DT_(0));
                    new_data.f_5 = new DenseVector<DT_>(data_size + 1, DT_(0));
                    new_data.f_6 = new DenseVector<DT_>(data_size + 1, DT_(0));
                    new_data.f_7 = new DenseVector<DT_>(data_size + 1, DT_(0));
                    new_data.f_8 = new DenseVector<DT_>(data_size + 1, DT_(0));

                    new_data.f_eq_0 = new DenseVector<DT_>(data_size + 1, DT_(0));
                    new_data.f_eq_1 = new DenseVector<DT_>(data_size + 1, DT_(0));
                    new_data.f_eq_2 = new DenseVector<DT_>(data_size + 1, DT_(0));
                    new_data.f_eq_3 = new DenseVector<DT_>(data_size + 1, DT_(0));
                    new_data.f_eq_4 = new DenseVector<DT_>(data_size + 1, DT_(0));
                    new_data.f_eq_5 = new DenseVector<DT_>(data_size + 1, DT_(0));
                    new_data.f_eq_6 = new DenseVector<DT_>(data_size + 1, DT_(0));
                    new_data.f_eq_7 = new DenseVector<DT_>(data_size + 1, DT_(0));
                    new_data.f_eq_8 = new DenseVector<DT_>(data_size + 1, DT_(0));

                    new_data.f_temp_0 = new DenseVector<DT_>(data_size + 1, DT_(0));
                    new_data.f_temp_1 = new DenseVector<DT_>(data_size + 1, DT_(0));
                    new_data.f_temp_2 = new DenseVector<DT_>(data_size + 1, DT_(0));
                    new_data.f_temp_3 = new DenseVector<DT_>(data_size + 1, DT_(0));
                    new_data.f_temp_4 = new DenseVector<DT_>(data_size + 1, DT_(0));
                    new_data.f_temp_5 = new DenseVector<DT_>(data_size + 1, DT_(0));
                    new_data.f_temp_6 = new DenseVector<DT_>(data_size + 1, DT_(0));
                    new_data.f_temp_7 = new DenseVector<DT_>(data_size + 1, DT_(0));
                    new_data.f_temp_8 = new DenseVector<DT_>(data_size + 1, DT_(0));

                    new_data.distribution_x = new DenseVector<DT_>(9ul, DT_(0));
                    new_data.distribution_y = new DenseVector<DT_>(9ul, DT_(0));

                    for (unsigned long j(0) ; j < temp_max_list[i] - temp_min_list[i]  + 1 ; ++j)
                    {
                        (*new_data.h)[j] = (*data.h)[j + info_list[i].offset];
                        (*new_data.u)[j] = (*data.u)[j + info_list[i].offset];
                        (*new_data.v)[j] = (*data.v)[j + info_list[i].offset];
                        (*new_data.b)[j] = (*data.b)[j + info_list[i].offset];
                    }
                    data_list.push_back(new_data);
                }

                // Compute fringes
                for (unsigned long i(0) ; i < info_list.size() ; ++i)
                {
                    PackedGridFringe<D2Q9> new_fringe;

                    _create_dir_fringe(i, *info_list[i].dir_index_1, *info_list[i].dir_1, info_list,
                            new_fringe.dir_index_1, new_fringe.dir_targets_1);
                    _create_dir_fringe(i, *info_list[i].dir_index_2, *info_list[i].dir_2, info_list,
                            new_fringe.dir_index_2, new_fringe.dir_targets_2);
                    _create_dir_fringe(i, *info_list[i].dir_index_3, *info_list[i].dir_3, info_list,
                            new_fringe.dir_index_3, new_fringe.dir_targets_3);
                    _create_dir_fringe(i, *info_list[i].dir_index_4, *info_list[i].dir_4, info_list,
                            new_fringe.dir_index_4, new_fringe.dir_targets_4);
                    _create_dir_fringe(i, *info_list[i].dir_index_5, *info_list[i].dir_5, info_list,
                            new_fringe.dir_index_5, new_fringe.dir_targets_5);
                    _create_dir_fringe(i, *info_list[i].dir_index_6, *info_list[i].dir_6, info_list,
                            new_fringe.dir_index_6, new_fringe.dir_targets_6);
                    _create_dir_fringe(i, *info_list[i].dir_index_7, *info_list[i].dir_7, info_list,
                            new_fringe.dir_index_7, new_fringe.dir_targets_7);
                    _create_dir_fringe(i, *info_list[i].dir_index_8, *info_list[i].dir_8, info_list,
                            new_fringe.dir_index_8, new_fringe.dir_targets_8);
                    _create_h_fringe(i, *info_list[i].limits, *data_list[i].h, info_list,
                            new_fringe.h_index, new_fringe.h_targets);
                    fringe_list.push_back(new_fringe);
                }
                // Add external_fringe data
                for (unsigned long i(0) ; i < info_list.size() ; ++i)
                {
                    _add_external_fringe(i, fringe_list);
                }

                /*std::cout<<std::endl;
                for (unsigned long i(0) ; i < temp_limits.size() ; ++i)
                {
                    std::cout<<temp_limits[i]<<" ";
                }
                std::cout<<std::endl;
                std::cout<<std::endl;
                for (unsigned long i(0) ; i < temp_limits.size() ; ++i)
                {
                    std::cout<<temp_dir_4[i]<<" ";
                }
                std::cout<<std::endl;*/
                /*for (unsigned long i(0) ; i < temp_limits.size() ; ++i)
                {
                    std::cout<<temp_types[i]<<" ";
                }
                std::cout<<std::endl;*/
                /*std::cout<<"barriers: ";
                for (unsigned long i(0) ; i < barriers.size() ; ++i)
                {
                    std::cout<<temp_limits[barriers[i]] <<"("<<barriers[i]<<")"<<" ";
                }
                std::cout<<std::endl;
                std::cout<<"neue limit vektoren: "<<std::endl;
                for (unsigned long i(0) ; i < temp_limits_list.size() ; ++i)
                {
                    for (unsigned long j(0); j < temp_limits_list[i].size() ; ++j)
                    {
                        std::cout<<temp_limits_list[i][j]<<" ";
                    }
                    std::cout<<std::endl;
                }
                std::cout<<std::endl;
                std::cout<<"min / max"<<std::endl;
                for (unsigned long i(0) ; i < temp_min_list.size() ; ++i)
                {
                    std::cout << temp_min_list[i] << "/";
                    std::cout << temp_max_list[i] << " ";
                    std::cout<<std::endl;
                }

                std::cout<<"echte limits vektoren: "<<std::endl;
                for (unsigned long i(0) ; i < info_list.size() ; ++i)
                {
                    for (unsigned long j(0); j < info_list[i].limits->size() ; ++j)
                    {
                        std::cout<<(*info_list[i].limits)[j]<<" ";
                    }
                    std::cout<<std::endl;
                }*/
            }

            static void recompose(PackedGridInfo<D2Q9> * info, PackedGridData<D2Q9, DT_> * data)
            {
                DenseVector<unsigned long> * tempul;

                tempul = new DenseVector<unsigned long>(info->limits->copy());
                delete info->limits;
                info->limits = tempul;

                tempul = new DenseVector<unsigned long>(info->types->copy());
                delete info->types;
                info->types = tempul;

                tempul = new DenseVector<unsigned long>(info->dir_1->copy());
                delete info->dir_1;
                info->dir_1 = tempul;

                tempul = new DenseVector<unsigned long>(info->dir_2->copy());
                delete info->dir_2;
                info->dir_2 = tempul;

                tempul = new DenseVector<unsigned long>(info->dir_3->copy());
                delete info->dir_3;
                info->dir_3 = tempul;

                tempul = new DenseVector<unsigned long>(info->dir_4->copy());
                delete info->dir_4;
                info->dir_4 = tempul;

                tempul = new DenseVector<unsigned long>(info->dir_5->copy());
                delete info->dir_5;
                info->dir_5 = tempul;

                tempul = new DenseVector<unsigned long>(info->dir_6->copy());
                delete info->dir_6;
                info->dir_6 = tempul;

                tempul = new DenseVector<unsigned long>(info->dir_7->copy());
                delete info->dir_7;
                info->dir_7 = tempul;

                tempul = new DenseVector<unsigned long>(info->dir_8->copy());
                delete info->dir_8;
                info->dir_8 = tempul;

                tempul = new DenseVector<unsigned long>(info->dir_index_1->copy());
                delete info->dir_index_1;
                info->dir_index_1 = tempul;

                tempul = new DenseVector<unsigned long>(info->dir_index_2->copy());
                delete info->dir_index_2;
                info->dir_index_2 = tempul;

                tempul = new DenseVector<unsigned long>(info->dir_index_3->copy());
                delete info->dir_index_3;
                info->dir_index_3 = tempul;

                tempul = new DenseVector<unsigned long>(info->dir_index_4->copy());
                delete info->dir_index_4;
                info->dir_index_4 = tempul;

                tempul = new DenseVector<unsigned long>(info->dir_index_5->copy());
                delete info->dir_index_5;
                info->dir_index_5 = tempul;

                tempul = new DenseVector<unsigned long>(info->dir_index_6->copy());
                delete info->dir_index_6;
                info->dir_index_6 = tempul;

                tempul = new DenseVector<unsigned long>(info->dir_index_7->copy());
                delete info->dir_index_7;
                info->dir_index_7 = tempul;

                tempul = new DenseVector<unsigned long>(info->dir_index_8->copy());
                delete info->dir_index_8;
                info->dir_index_8 = tempul;

                DenseVector<DT_> * temp;

                temp = new DenseVector<DT_>(data->h->copy());
                delete data->h;
                data->h = temp;

                temp = new DenseVector<DT_>(data->u->copy());
                delete data->u;
                data->u = temp;

                temp = new DenseVector<DT_>(data->v->copy());
                delete data->v;
                data->v = temp;

                temp = new DenseVector<DT_>(data->b->copy());
                delete data->b;
                data->b = temp;

                temp = new DenseVector<DT_>(data->f_0->copy());
                delete data->f_0;
                data->f_0 = temp;

                temp = new DenseVector<DT_>(data->f_1->copy());
                delete data->f_1;
                data->f_1 = temp;

                temp = new DenseVector<DT_>(data->f_2->copy());
                delete data->f_2;
                data->f_2 = temp;

                temp = new DenseVector<DT_>(data->f_3->copy());
                delete data->f_3;
                data->f_3 = temp;

                temp = new DenseVector<DT_>(data->f_4->copy());
                delete data->f_4;
                data->f_4 = temp;

                temp = new DenseVector<DT_>(data->f_5->copy());
                delete data->f_5;
                data->f_5 = temp;

                temp = new DenseVector<DT_>(data->f_6->copy());
                delete data->f_6;
                data->f_6 = temp;

                temp = new DenseVector<DT_>(data->f_7->copy());
                delete data->f_7;
                data->f_7 = temp;

                temp = new DenseVector<DT_>(data->f_8->copy());
                delete data->f_8;
                data->f_8 = temp;

                temp = new DenseVector<DT_>(data->f_eq_0->copy());
                delete data->f_eq_0;
                data->f_eq_0 = temp;

                temp = new DenseVector<DT_>(data->f_eq_1->copy());
                delete data->f_eq_1;
                data->f_eq_1 = temp;

                temp = new DenseVector<DT_>(data->f_eq_2->copy());
                delete data->f_eq_2;
                data->f_eq_2 = temp;

                temp = new DenseVector<DT_>(data->f_eq_3->copy());
                delete data->f_eq_3;
                data->f_eq_3 = temp;

                temp = new DenseVector<DT_>(data->f_eq_4->copy());
                delete data->f_eq_4;
                data->f_eq_4 = temp;

                temp = new DenseVector<DT_>(data->f_eq_5->copy());
                delete data->f_eq_5;
                data->f_eq_5 = temp;

                temp = new DenseVector<DT_>(data->f_eq_6->copy());
                delete data->f_eq_6;
                data->f_eq_6 = temp;

                temp = new DenseVector<DT_>(data->f_eq_7->copy());
                delete data->f_eq_7;
                data->f_eq_7 = temp;

                temp = new DenseVector<DT_>(data->f_eq_8->copy());
                delete data->f_eq_8;
                data->f_eq_8 = temp;

                temp = new DenseVector<DT_>(data->f_temp_0->copy());
                delete data->f_temp_0;
                data->f_temp_0 = temp;

                temp = new DenseVector<DT_>(data->f_temp_1->copy());
                delete data->f_temp_1;
                data->f_temp_1 = temp;

                temp = new DenseVector<DT_>(data->f_temp_2->copy());
                delete data->f_temp_2;
                data->f_temp_2 = temp;

                temp = new DenseVector<DT_>(data->f_temp_3->copy());
                delete data->f_temp_3;
                data->f_temp_3 = temp;

                temp = new DenseVector<DT_>(data->f_temp_4->copy());
                delete data->f_temp_4;
                data->f_temp_4 = temp;

                temp = new DenseVector<DT_>(data->f_temp_5->copy());
                delete data->f_temp_5;
                data->f_temp_5 = temp;

                temp = new DenseVector<DT_>(data->f_temp_6->copy());
                delete data->f_temp_6;
                data->f_temp_6 = temp;

                temp = new DenseVector<DT_>(data->f_temp_7->copy());
                delete data->f_temp_7;
                data->f_temp_7 = temp;

                temp = new DenseVector<DT_>(data->f_temp_8->copy());
                delete data->f_temp_8;
                data->f_temp_8 = temp;

                temp = new DenseVector<DT_>(data->distribution_x->copy());
                delete data->distribution_x;
                data->distribution_x = temp;

                temp = new DenseVector<DT_>(data->distribution_y->copy());
                delete data->distribution_y;
                data->distribution_y = temp;
            }
    };
}

#endif
