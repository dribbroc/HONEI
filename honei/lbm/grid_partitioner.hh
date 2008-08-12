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
#include <list>
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
            static void partition(unsigned long parts, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data,
                    std::list<PackedGridInfo<D2Q9> > & info_list, std::list<PackedGridData<D2Q9, DT_> > & data_list)
            {
                std::vector<unsigned long> temp_limits;
                std::vector<unsigned long> temp_types;
                std::vector<unsigned long> temp_dir_0;
                std::vector<unsigned long> temp_dir_1;
                std::vector<unsigned long> temp_dir_2;
                std::vector<unsigned long> temp_dir_3;
                std::vector<unsigned long> temp_dir_4;
                std::vector<unsigned long> temp_dir_5;
                std::vector<unsigned long> temp_dir_6;
                std::vector<unsigned long> temp_dir_7;
                std::vector<unsigned long> temp_dir_8;

                for (unsigned long i(0) ; i < info.limits->size() ; ++i)
                {
                    temp_limits.push_back((*info.limits)[i]);
                    temp_types.push_back((*info.types)[i]);
                    temp_dir_0.push_back((*info.dir_0)[i]);
                    temp_dir_1.push_back((*info.dir_1)[i]);
                    temp_dir_2.push_back((*info.dir_2)[i]);
                    temp_dir_3.push_back((*info.dir_3)[i]);
                    temp_dir_4.push_back((*info.dir_4)[i]);
                    temp_dir_5.push_back((*info.dir_5)[i]);
                    temp_dir_6.push_back((*info.dir_6)[i]);
                    temp_dir_7.push_back((*info.dir_7)[i]);
                    temp_dir_8.push_back((*info.dir_8)[i]);
                }
                /*std::cout<<std::endl;
                for (unsigned long i(0) ; i < temp_limits.size() ; ++i)
                {
                    std::cout<<temp_limits[i]<<" ";
                }
                std::cout<<std::endl;
                for (unsigned long i(0) ; i < temp_limits.size() ; ++i)
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
                std::vector<unsigned long> barriers;
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
                        std::vector<unsigned long>::iterator it_dir_0(temp_dir_0.begin());
                        std::vector<unsigned long>::iterator it_dir_1(temp_dir_1.begin());
                        std::vector<unsigned long>::iterator it_dir_2(temp_dir_2.begin());
                        std::vector<unsigned long>::iterator it_dir_3(temp_dir_3.begin());
                        std::vector<unsigned long>::iterator it_dir_4(temp_dir_4.begin());
                        std::vector<unsigned long>::iterator it_dir_5(temp_dir_5.begin());
                        std::vector<unsigned long>::iterator it_dir_6(temp_dir_6.begin());
                        std::vector<unsigned long>::iterator it_dir_7(temp_dir_7.begin());
                        std::vector<unsigned long>::iterator it_dir_8(temp_dir_8.begin());
                        while (i != temp_limits.size() - 1 && temp_limits[i] != end)
                        {
                            ++i;
                            ++it_types;
                            ++it_dir_0;
                            ++it_dir_1;
                            ++it_dir_2;
                            ++it_dir_3;
                            ++it_dir_4;
                            ++it_dir_5;
                            ++it_dir_6;
                            ++it_dir_7;
                            ++it_dir_8;
                        }
                        barriers.push_back(i);
                        temp_types.insert(it_types, *(--it_types));
                        temp_dir_0.insert(it_dir_0, *(--it_dir_0));
                        temp_dir_1.insert(it_dir_1, *(--it_dir_1));
                        temp_dir_2.insert(it_dir_2, *(--it_dir_2));
                        temp_dir_3.insert(it_dir_3, *(--it_dir_3));
                        temp_dir_4.insert(it_dir_4, *(--it_dir_4));
                        temp_dir_5.insert(it_dir_5, *(--it_dir_5));
                        temp_dir_6.insert(it_dir_6, *(--it_dir_6));
                        temp_dir_7.insert(it_dir_7, *(--it_dir_7));
                        temp_dir_8.insert(it_dir_8, *(--it_dir_8));
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
                    start = end;
                }

                // Create partitions
                std::vector<std::vector<unsigned long> > temp_limits_list;
                std::vector<std::vector<unsigned long> > temp_types_list;
                std::vector<std::vector<unsigned long> > temp_dir_0_list;
                std::vector<std::vector<unsigned long> > temp_dir_1_list;
                std::vector<std::vector<unsigned long> > temp_dir_2_list;
                std::vector<std::vector<unsigned long> > temp_dir_3_list;
                std::vector<std::vector<unsigned long> > temp_dir_4_list;
                std::vector<std::vector<unsigned long> > temp_dir_5_list;
                std::vector<std::vector<unsigned long> > temp_dir_6_list;
                std::vector<std::vector<unsigned long> > temp_dir_7_list;
                std::vector<std::vector<unsigned long> > temp_dir_8_list;

                for (unsigned long i(0) ; i < barriers.size() - 1; ++i)
                {
                    std::vector<unsigned long> new_limits;
                    std::vector<unsigned long> new_types;
                    std::vector<unsigned long> new_dir_0;
                    std::vector<unsigned long> new_dir_1;
                    std::vector<unsigned long> new_dir_2;
                    std::vector<unsigned long> new_dir_3;
                    std::vector<unsigned long> new_dir_4;
                    std::vector<unsigned long> new_dir_5;
                    std::vector<unsigned long> new_dir_6;
                    std::vector<unsigned long> new_dir_7;
                    std::vector<unsigned long> new_dir_8;
                    for(unsigned long index(barriers[i]) ; index <= barriers[i + 1] ; ++index)
                    {
                        new_limits.push_back(temp_limits[index]);
                        new_types.push_back(temp_types[index]);
                        new_dir_0.push_back(temp_dir_0[index]);
                        new_dir_1.push_back(temp_dir_1[index]);
                        new_dir_2.push_back(temp_dir_2[index]);
                        new_dir_3.push_back(temp_dir_3[index]);
                        new_dir_4.push_back(temp_dir_4[index]);
                        new_dir_5.push_back(temp_dir_5[index]);
                        new_dir_6.push_back(temp_dir_6[index]);
                        new_dir_7.push_back(temp_dir_7[index]);
                        new_dir_8.push_back(temp_dir_8[index]);
                    }
                    temp_limits_list.push_back(new_limits);
                    temp_types_list.push_back(new_types);
                    temp_dir_0_list.push_back(new_dir_0);
                    temp_dir_1_list.push_back(new_dir_1);
                    temp_dir_2_list.push_back(new_dir_2);
                    temp_dir_3_list.push_back(new_dir_3);
                    temp_dir_4_list.push_back(new_dir_4);
                    temp_dir_5_list.push_back(new_dir_5);
                    temp_dir_6_list.push_back(new_dir_6);
                    temp_dir_7_list.push_back(new_dir_7);
                    temp_dir_8_list.push_back(new_dir_8);
                }

                /*std::cout<<std::endl;
                for (unsigned long i(0) ; i < temp_limits.size() ; ++i)
                {
                    std::cout<<temp_limits[i]<<" ";
                }
                std::cout<<std::endl;
                for (unsigned long i(0) ; i < temp_limits.size() ; ++i)
                {
                    std::cout<<temp_types[i]<<" ";
                }
                std::cout<<std::endl;
                std::cout<<"barriers: ";
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
                std::cout<<"neue type vektoren: "<<std::endl;
                for (unsigned long i(0) ; i < temp_types_list.size() ; ++i)
                {
                    for (unsigned long j(0); j < temp_types_list[i].size() ; ++j)
                    {
                        std::cout<<temp_types_list[i][j]<<" ";
                    }
                    std::cout<<std::endl;
                }*/
            }
    };
}

#endif
