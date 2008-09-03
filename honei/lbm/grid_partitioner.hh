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
#include <vector>
#include <iostream>
#include <limits>

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
        public:

            static void compose(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data,
                    std::vector<PackedGridInfo<D2Q9> > & info_list, std::vector<PackedGridData<D2Q9, DT_> > & data_list)
            {
                /*for (unsigned long i(0); i < data_list[0].h->size(); ++i)
                {
                    (*data.h)[i] = (*data_list[0].h)[i];
                }
                for(unsigned long index(1); index < data_list.size(); ++index)
                {
                    for(unsigned long i(0); i < data_list[index].h->size(); ++i)
                    {
                        (*data.h)[i + info_list[index].offset] += (*data_list[index].h)[i];
                    }
                }*/

                for(unsigned long i(0); i < data.h->size(); ++i)
                {
                    unsigned long index(0);
                    while(info_list[index].offset + data_list[index].h->size() <= i)
                    {
                        ++index;
                    }
                    unsigned long start(index);
                    while(index < info_list.size() && info_list.at(index).offset <= i)
                    {
                        ++index;
                    }
                    unsigned long end(index);
                    (*data.h)[i] = (*data_list[start].h)[i - (info_list[start].offset)];
                }
            }

            static void synch(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data,
                    std::vector<PackedGridInfo<D2Q9> > & info_list, std::vector<PackedGridData<D2Q9, DT_> > & data_list)
            {
                for(unsigned long i(0); i < data.h->size(); ++i)
                {
                    unsigned long index(0);
                    while(info_list[index].offset + data_list[index].h->size() <= i)
                    {
                        ++index;
                    }
                    unsigned long start(index);
                    while(index < info_list.size() && info_list[index].offset <= i)
                    {
                        ++index;
                    }
                    unsigned long end(index);
                    DT_ values_h[end - start];
                    DT_ values_u[end - start];
                    DT_ values_v[end - start];

                    DT_ values_f_0[end - start];
                    DT_ values_f_1[end - start];
                    DT_ values_f_2[end - start];
                    DT_ values_f_3[end - start];
                    DT_ values_f_4[end - start];
                    DT_ values_f_5[end - start];
                    DT_ values_f_6[end - start];
                    DT_ values_f_7[end - start];
                    DT_ values_f_8[end - start];
                    DT_ values_f_9[end - start];

                    DT_ values_f_eq_0[end - start];
                    DT_ values_f_eq_1[end - start];
                    DT_ values_f_eq_2[end - start];
                    DT_ values_f_eq_3[end - start];
                    DT_ values_f_eq_4[end - start];
                    DT_ values_f_eq_5[end - start];
                    DT_ values_f_eq_6[end - start];
                    DT_ values_f_eq_7[end - start];
                    DT_ values_f_eq_8[end - start];

                    for(unsigned long l(start); l < end; ++l)
                    {
                        values_h[l - start] = (*data_list[l].h)[i - (info_list[l].offset)];
                        values_u[l - start] = (*data_list[l].u)[i - (info_list[l].offset)];
                        values_v[l - start] = (*data_list[l].v)[i - (info_list[l].offset)];

                        values_f_0[l - start] = (*data_list[l].f_0)[i - (info_list[l].offset)];
                        values_f_1[l - start] = (*data_list[l].f_1)[i - (info_list[l].offset)];
                        values_f_2[l - start] = (*data_list[l].f_2)[i - (info_list[l].offset)];
                        values_f_3[l - start] = (*data_list[l].f_3)[i - (info_list[l].offset)];
                        values_f_4[l - start] = (*data_list[l].f_4)[i - (info_list[l].offset)];
                        values_f_5[l - start] = (*data_list[l].f_5)[i - (info_list[l].offset)];
                        values_f_6[l - start] = (*data_list[l].f_6)[i - (info_list[l].offset)];
                        values_f_7[l - start] = (*data_list[l].f_7)[i - (info_list[l].offset)];
                        values_f_8[l - start] = (*data_list[l].f_8)[i - (info_list[l].offset)];

                        values_f_eq_0[l - start] = (*data_list[l].f_eq_0)[i - (info_list[l].offset)];
                        values_f_eq_1[l - start] = (*data_list[l].f_eq_1)[i - (info_list[l].offset)];
                        values_f_eq_2[l - start] = (*data_list[l].f_eq_2)[i - (info_list[l].offset)];
                        values_f_eq_3[l - start] = (*data_list[l].f_eq_3)[i - (info_list[l].offset)];
                        values_f_eq_4[l - start] = (*data_list[l].f_eq_4)[i - (info_list[l].offset)];
                        values_f_eq_5[l - start] = (*data_list[l].f_eq_5)[i - (info_list[l].offset)];
                        values_f_eq_6[l - start] = (*data_list[l].f_eq_6)[i - (info_list[l].offset)];
                        values_f_eq_7[l - start] = (*data_list[l].f_eq_7)[i - (info_list[l].offset)];
                        values_f_eq_8[l - start] = (*data_list[l].f_eq_8)[i - (info_list[l].offset)];
                    }

                    for(unsigned long k(start); k < end; ++k)
                    {
                        (*data_list[k].h)[i - (info_list[k].offset)] = values_h[0];
                        (*data_list[k].u)[i - (info_list[k].offset)] = values_u[0];
                        (*data_list[k].v)[i - (info_list[k].offset)] = values_v[0];

                        (*data_list[k].f_0)[i - (info_list[k].offset)] = values_f_0[0];
                        (*data_list[k].f_1)[i - (info_list[k].offset)] = values_f_1[0];
                        (*data_list[k].f_2)[i - (info_list[k].offset)] = values_f_2[0];
                        (*data_list[k].f_3)[i - (info_list[k].offset)] = values_f_3[0];
                        (*data_list[k].f_4)[i - (info_list[k].offset)] = values_f_4[0];
                        (*data_list[k].f_5)[i - (info_list[k].offset)] = values_f_5[0];
                        (*data_list[k].f_6)[i - (info_list[k].offset)] = values_f_6[0];
                        (*data_list[k].f_7)[i - (info_list[k].offset)] = values_f_7[0];
                        (*data_list[k].f_8)[i - (info_list[k].offset)] = values_f_8[0];


                        (*data_list[k].f_eq_0)[i - (info_list[k].offset)] = values_f_eq_0[0];
                        (*data_list[k].f_eq_1)[i - (info_list[k].offset)] = values_f_eq_1[0];
                        (*data_list[k].f_eq_2)[i - (info_list[k].offset)] = values_f_eq_2[0];
                        (*data_list[k].f_eq_3)[i - (info_list[k].offset)] = values_f_eq_3[0];
                        (*data_list[k].f_eq_4)[i - (info_list[k].offset)] = values_f_eq_4[0];
                        (*data_list[k].f_eq_5)[i - (info_list[k].offset)] = values_f_eq_5[0];
                        (*data_list[k].f_eq_6)[i - (info_list[k].offset)] = values_f_eq_6[0];
                        (*data_list[k].f_eq_7)[i - (info_list[k].offset)] = values_f_eq_7[0];
                        (*data_list[k].f_eq_8)[i - (info_list[k].offset)] = values_f_eq_8[0];

                        for(unsigned long j(1); j < end - start; ++j)
                        {
                            (*data_list[k].h)[i - (info_list[k].offset)] += values_h[j];
                            (*data_list[k].u)[i - (info_list[k].offset)] += values_u[j];
                            (*data_list[k].v)[i - (info_list[k].offset)] += values_v[j];

                            (*data_list[k].f_0)[i - (info_list[k].offset)] += values_f_0[j];
                            (*data_list[k].f_1)[i - (info_list[k].offset)] += values_f_1[j];
                            (*data_list[k].f_2)[i - (info_list[k].offset)] += values_f_2[j];
                            (*data_list[k].f_3)[i - (info_list[k].offset)] += values_f_3[j];
                            (*data_list[k].f_4)[i - (info_list[k].offset)] += values_f_4[j];
                            (*data_list[k].f_5)[i - (info_list[k].offset)] += values_f_5[j];
                            (*data_list[k].f_6)[i - (info_list[k].offset)] += values_f_6[j];
                            (*data_list[k].f_7)[i - (info_list[k].offset)] += values_f_7[j];
                            (*data_list[k].f_8)[i - (info_list[k].offset)] += values_f_8[j];

                            (*data_list[k].f_eq_0)[i - (info_list[k].offset)] += values_f_eq_0[j];
                            (*data_list[k].f_eq_1)[i - (info_list[k].offset)] += values_f_eq_1[j];
                            (*data_list[k].f_eq_2)[i - (info_list[k].offset)] += values_f_eq_2[j];
                            (*data_list[k].f_eq_3)[i - (info_list[k].offset)] += values_f_eq_3[j];
                            (*data_list[k].f_eq_4)[i - (info_list[k].offset)] += values_f_eq_4[j];
                            (*data_list[k].f_eq_5)[i - (info_list[k].offset)] += values_f_eq_5[j];
                            (*data_list[k].f_eq_6)[i - (info_list[k].offset)] += values_f_eq_6[j];
                            (*data_list[k].f_eq_7)[i - (info_list[k].offset)] += values_f_eq_7[j];
                            (*data_list[k].f_eq_8)[i - (info_list[k].offset)] += values_f_eq_8[j];
                        }
                    }
                }
            }

            static void partition_directions(std::vector<unsigned long> & dir_index, std::vector<unsigned long> & dir,
                    std::vector<unsigned long> & barriers, unsigned long start, unsigned long end)
            {
                CONTEXT("When partitioning one single direction:");
                /// todo bei abbruch dennoch barriers einbauen
                if (dir_index.size() == 0)
                    return;
                if (dir_index[dir_index.size() - 1] < start)
                    return;
                if (dir_index[0] > end)
                    return;

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
                if (dir_index[end_index] > end)
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
                    std::vector<PackedGridInfo<D2Q9> > & info_list, std::vector<PackedGridData<D2Q9, DT_> > & data_list)
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

                /*std::cout<<std::endl;
                for (unsigned long i(0) ; i < temp_limits.size() ; ++i)
                {
                    std::cout<<temp_limits[i]<<" ";
                }
                std::cout<<std::endl;*/
                /*for (unsigned long i(0) ; i < temp_limits.size() ; ++i)
                {
                    std::cout<<temp_types[i]<<" ";
                }
                std::cout<<std::endl;*/

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
                        std::vector<unsigned long>::iterator one_before(end_index);
                        --one_before;
                        // Insert dummy nodes
                        temp_limits.insert(end_index, end);

                        unsigned long one_before_offset(end - *one_before);
                        unsigned long i(0);
                        std::vector<unsigned long>::iterator it_types(temp_types.begin());
                        /*std::vector<unsigned long>::iterator it_dir_1(temp_dir_1.begin());
                        std::vector<unsigned long>::iterator it_dir_2(temp_dir_2.begin());
                        std::vector<unsigned long>::iterator it_dir_3(temp_dir_3.begin());
                        std::vector<unsigned long>::iterator it_dir_4(temp_dir_4.begin());
                        std::vector<unsigned long>::iterator it_dir_5(temp_dir_5.begin());
                        std::vector<unsigned long>::iterator it_dir_6(temp_dir_6.begin());
                        std::vector<unsigned long>::iterator it_dir_7(temp_dir_7.begin());
                        std::vector<unsigned long>::iterator it_dir_8(temp_dir_8.begin());*/
                        while (i != temp_limits.size() - 1 && temp_limits[i] != end)
                        {
                            ++i;
                            ++it_types;
                            /*++it_dir_1;
                            ++it_dir_2;
                            ++it_dir_3;
                            ++it_dir_4;
                            ++it_dir_5;
                            ++it_dir_6;
                            ++it_dir_7;
                            ++it_dir_8;*/
                        }
                        barriers.push_back(i);
                        temp_types.insert(it_types, *(--it_types));
                        /*temp_dir_1.insert(++it_dir_1, *(--it_dir_1) + one_before_offset);
                        temp_dir_2.insert(++it_dir_2, *(--it_dir_2) + one_before_offset);
                        temp_dir_3.insert(++it_dir_3, *(--it_dir_3) + one_before_offset);
                        temp_dir_4.insert(++it_dir_4, *(--it_dir_4) + one_before_offset);
                        temp_dir_5.insert(++it_dir_5, *(--it_dir_5) + one_before_offset);
                        temp_dir_6.insert(++it_dir_6, *(--it_dir_6) + one_before_offset);
                        temp_dir_7.insert(++it_dir_7, *(--it_dir_7) + one_before_offset);
                        temp_dir_8.insert(++it_dir_8, *(--it_dir_8) + one_before_offset);*/
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

                std::cout<<"barriers 0 vs 6: "<<barriers.size() <<" "<<barriers_6.size()<<std::endl;
                for (unsigned long i(0) ; i < barriers.size() ; ++i)
                    std::cout<<barriers[i]<<" ";
                std::cout<<std::endl;
                for (unsigned long i(0) ; i < barriers_6.size() ; ++i)
                    std::cout<<barriers_6[i]<<" ";
                std::cout<<std::endl;
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

                    for(unsigned long index(barriers[i]) ; index <= barriers[i + 1] ; ++index)
                    {
                        new_limits.push_back(temp_limits[index]);
                        new_types.push_back(temp_types[index]);
                        /*new_dir_1.push_back(temp_dir_1[index]);
                        new_dir_2.push_back(temp_dir_2[index]);
                        new_dir_3.push_back(temp_dir_3[index]);
                        new_dir_4.push_back(temp_dir_4[index]);
                        new_dir_5.push_back(temp_dir_5[index]);
                        new_dir_6.push_back(temp_dir_6[index]);
                        new_dir_7.push_back(temp_dir_7[index]);
                        new_dir_8.push_back(temp_dir_8[index]);*/
                    }
                    for (unsigned long index(barriers_1[i * 2]) ; index <= barriers_1[i * 2 + 1] ; ++index)
                    {
                        new_dir_index_1.push_back(temp_dir_index_1[index]);
                        if (index % 2 == 0)
                            new_dir_1.push_back(temp_dir_1[index / 2]);
                    }
                    for (unsigned long index(barriers_2[i * 2]) ; index <= barriers_2[i * 2 + 1] ; ++index)
                    {
                        new_dir_index_2.push_back(temp_dir_index_2[index]);
                        if (index % 2 == 0)
                            new_dir_2.push_back(temp_dir_2[index / 2]);
                    }
                    for (unsigned long index(barriers_3[i * 2]) ; index <= barriers_3[i * 2 + 1] ; ++index)
                    {
                        new_dir_index_3.push_back(temp_dir_index_3[index]);
                        if (index % 2 == 0)
                            new_dir_3.push_back(temp_dir_3[index / 2]);
                    }
                    for (unsigned long index(barriers_4[i * 2]) ; index <= barriers_4[i * 2 + 1] ; ++index)
                    {
                        new_dir_index_4.push_back(temp_dir_index_4[index]);
                        if (index % 2 == 0)
                            new_dir_4.push_back(temp_dir_4[index / 2]);
                    }
                    for (unsigned long index(barriers_5[i * 2]) ; index <= barriers_5[i * 2 + 1] ; ++index)
                    {
                        new_dir_index_5.push_back(temp_dir_index_5[index]);
                        if (index % 2 == 0)
                            new_dir_5.push_back(temp_dir_5[index / 2]);
                    }
                    for (unsigned long index(barriers_6[i * 2]) ; index <= barriers_6[i * 2 + 1] ; ++index)
                    {
                        new_dir_index_6.push_back(temp_dir_index_6[index]);
                        if (index % 2 == 0)
                            new_dir_6.push_back(temp_dir_6[index / 2]);
                    }
                    for (unsigned long index(barriers_7[i * 2]) ; index <= barriers_7[i * 2 + 1] ; ++index)
                    {
                        new_dir_index_7.push_back(temp_dir_index_7[index]);
                        if (index % 2 == 0)
                            new_dir_7.push_back(temp_dir_7[index / 2]);
                    }
                    for (unsigned long index(barriers_8[i * 2]) ; index <= barriers_8[i * 2 + 1] ; ++index)
                    {
                        new_dir_index_8.push_back(temp_dir_index_8[index]);
                        if (index % 2 == 0)
                            new_dir_8.push_back(temp_dir_8[index / 2]);
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

                    /// \todo max von dir muss jeweils noch + "offset" gerechnet werden
                    std::vector<unsigned long> max_collection;
                    max_collection.push_back(*max_element(new_dir_1.begin(), new_dir_1.end()));
                    max_collection.push_back(*max_element(new_dir_2.begin(), new_dir_2.end()));
                    max_collection.push_back(*max_element(new_dir_3.begin(), new_dir_3.end()));
                    max_collection.push_back(*max_element(new_dir_4.begin(), new_dir_4.end()));
                    max_collection.push_back(*max_element(new_dir_5.begin(), new_dir_5.end()));
                    max_collection.push_back(*max_element(new_dir_6.begin(), new_dir_6.end()));
                    max_collection.push_back(*max_element(new_dir_7.begin(), new_dir_7.end()));
                    max_collection.push_back(*max_element(new_dir_8.begin(), new_dir_8.end()));
                    temp_max_list.push_back(*max_element(max_collection.begin(), max_collection.end()));

                    std::vector<unsigned long> min_collection;
                    min_collection.push_back(*min_element(new_dir_1.begin(), new_dir_1.end()));
                    min_collection.push_back(*min_element(new_dir_2.begin(), new_dir_2.end()));
                    min_collection.push_back(*min_element(new_dir_3.begin(), new_dir_3.end()));
                    min_collection.push_back(*min_element(new_dir_4.begin(), new_dir_4.end()));
                    min_collection.push_back(*min_element(new_dir_5.begin(), new_dir_5.end()));
                    min_collection.push_back(*min_element(new_dir_6.begin(), new_dir_6.end()));
                    min_collection.push_back(*min_element(new_dir_7.begin(), new_dir_7.end()));
                    min_collection.push_back(*min_element(new_dir_8.begin(), new_dir_8.end()));
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
                        (*new_info.types)[j] = temp_types_list[i][j] - new_info.offset;
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
                    PackedGridData<D2Q9,DT_> new_data;
                    new_data.h = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1);
                    new_data.u = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1);
                    new_data.v = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1);

                    new_data.f_0 = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1, DT_(0));
                    new_data.f_1 = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1, DT_(0));
                    new_data.f_2 = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1, DT_(0));
                    new_data.f_3 = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1, DT_(0));
                    new_data.f_4 = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1, DT_(0));
                    new_data.f_5 = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1, DT_(0));
                    new_data.f_6 = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1, DT_(0));
                    new_data.f_7 = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1, DT_(0));
                    new_data.f_8 = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1, DT_(0));

                    new_data.f_eq_0 = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1, DT_(0));
                    new_data.f_eq_1 = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1, DT_(0));
                    new_data.f_eq_2 = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1, DT_(0));
                    new_data.f_eq_3 = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1, DT_(0));
                    new_data.f_eq_4 = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1, DT_(0));
                    new_data.f_eq_5 = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1, DT_(0));
                    new_data.f_eq_6 = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1, DT_(0));
                    new_data.f_eq_7 = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1, DT_(0));
                    new_data.f_eq_8 = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1, DT_(0));

                    new_data.f_temp_0 = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1, DT_(0));
                    new_data.f_temp_1 = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1, DT_(0));
                    new_data.f_temp_2 = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1, DT_(0));
                    new_data.f_temp_3 = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1, DT_(0));
                    new_data.f_temp_4 = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1, DT_(0));
                    new_data.f_temp_5 = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1, DT_(0));
                    new_data.f_temp_6 = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1, DT_(0));
                    new_data.f_temp_7 = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1, DT_(0));
                    new_data.f_temp_8 = new DenseVector<DT_>(temp_max_list[i] - temp_min_list[i] + 1, DT_(0));
                    new_data.distribution_x = new DenseVector<DT_>(9ul, DT_(0));
                    new_data.distribution_y = new DenseVector<DT_>(9ul, DT_(0));

                    for (unsigned long j(0) ; j < temp_max_list[i] - temp_min_list[i] + 1 ; ++j)
                    {
                        (*new_data.h)[j] = (*data.h)[j + info_list[i].offset];
                        (*new_data.u)[j] = (*data.u)[j + info_list[i].offset];
                        (*new_data.v)[j] = (*data.v)[j + info_list[i].offset];
                    }
                    data_list.push_back(new_data);
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
    };
}

#endif
