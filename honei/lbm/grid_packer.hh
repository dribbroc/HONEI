/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
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


#pragma once
#ifndef LBM_GUARD_GRID_PACKER_HH
#define LBM_GUARD_GRID_PACKER_HH 1

#include <honei/lbm/grid.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/dense_matrix.hh>
#include <vector>
#include <honei/util/attributes.hh>
#include <iostream>
/**
 * \file
 * Definition of LBM grid packer.
 *
 * \ingroup grpliblbm
 **/

using namespace honei;
using namespace lbm;
using namespace lbm_lattice_types;

namespace honei
{
    template <typename LatticeType_, typename BoundaryType_, typename DT_> struct GridPacker
    {
    };

    template <typename DT_> struct GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DT_>
    {
        private:
            static void _dir_post_pack(std::vector<unsigned long> & dir, std::vector<unsigned long> & dir_post,
                    std::vector<unsigned long> & dir_index, std::vector<unsigned long> & limits)
            {
                for (unsigned long i(0) ; i < limits.size() - 1; ++i)
                {
                    // if dir points not to our own cell (reflection)
                    if (dir[i] != limits[i])
                    {
                        dir_post.push_back(dir[i]);
                        dir_index.push_back(limits.at(i));
                        // search for the sequence end
                        unsigned long end(i + 1);
                        while(end < limits.size() - 1 && (dir[end] - dir[i]) == (limits[end] - limits[i]))
                        {
                            ++end;
                        }
                        dir_index.push_back(limits.at(end));
                        i = end - 1;
                    }
                }
            }

            static unsigned long _element_type(unsigned long i, unsigned long j, Grid<D2Q9, DT_> & grid)
            {
                unsigned long type(0);

                // Given DIR_N:
                // type &= 1<<(N-1)

                // Outer Boundaries
                if(i == 0)
                {
                    type |= 1<<1;
                    type |= 1<<2;
                    type |= 1<<3;
                }

                if(i == (*grid.obstacles).rows() - 1)
                {
                    type |= 1<<5;
                    type |= 1<<6;
                    type |= 1<<7;
                }

                if(j == 0)
                {
                    type |= 1<<3;
                    type |= 1<<4;
                    type |= 1<<5;
                }

                if(j == (*grid.obstacles).columns() - 1)
                {
                    type |= 1<<0;
                    type |= 1<<1;
                    type |= 1<<7;
                }


                // Inner Boundaries
                if(i > 0 && (*grid.obstacles)(i - 1, j))
                    type |= 1<<2;

                if(i <  (*grid.obstacles).rows() - 1 && (*grid.obstacles)(i + 1, j))
                    type |= 1<<6;

                if(j > 0 && (*grid.obstacles)(i, j - 1))
                    type |= 1<<4;

                if(j <  (*grid.obstacles).columns() - 1 && (*grid.obstacles)(i, j + 1))
                    type |= 1<<0;

                // Diagonal Boundaries
                if (i > 0 && j <  (*grid.obstacles).columns() - 1 && (*grid.obstacles)(i - 1, j + 1))
                    type |= 1<<1;

                if (i > 0 && j > 0 && (*grid.obstacles)(i - 1, j - 1))
                    type |= 1<<3;

                if (i < (*grid.obstacles).rows() - 1 && j > 0 && (*grid.obstacles)(i + 1, j - 1))
                    type |= 1<<5;

                if (i < (*grid.obstacles).rows() - 1 && j < (*grid.obstacles).columns() - 1 && (*grid.obstacles)(i + 1, j + 1))
                    type |= 1<<7;

                return type;
            }

            static void _element_direction(unsigned long packed_index, unsigned long i, unsigned long j, Grid<D2Q9, DT_> & grid,
                    std::vector<unsigned long> & dir_1,
                    std::vector<unsigned long> & dir_2,
                    std::vector<unsigned long> & dir_3,
                    std::vector<unsigned long> & dir_4,
                    std::vector<unsigned long> & dir_5,
                    std::vector<unsigned long> & dir_6,
                    std::vector<unsigned long> & dir_7,
                    std::vector<unsigned long> & dir_8)
            {
                /* D2Q9 according to Markus
                        3
                      4 . 2
                        .
                     5..0..1
                        .
                      6 . 8
                        7
                 */
                unsigned long type(_element_type(i, j, grid));

                // DIR_1
                if ((type & 1<<0) == 1<<0)
                    dir_1.push_back(packed_index);
                else
                    dir_1.push_back(packed_index + 1);


                // DIR_5
                if ((type & 1<<4) == 1<<4)
                    dir_5.push_back(packed_index);
                else
                    dir_5.push_back(packed_index - 1);


                // DIR_7 DIR_6 DIR_8
                if (i == grid.obstacles->rows() - 1)
                {
                    dir_7.push_back(packed_index);
                    dir_6.push_back(packed_index);
                    dir_8.push_back(packed_index);
                }
                else
                {
                    DenseMatrix<bool>::ElementIterator it(grid.obstacles->element_at(i * grid.obstacles->columns() + j));
                    unsigned long temp_index(packed_index);
                    while (it.index() != (i + 1) * grid.obstacles->columns() + j)
                    {
                        if (!*it)
                        {
                            ++temp_index;
                        }
                        ++it;
                    }

                    // DIR_7
                    if ((type & 1<<6) == 1<<6)
                        dir_7.push_back(packed_index);
                    else
                        dir_7.push_back(temp_index);

                    // DIR_6
                    if ((type & 1<<5) == 1<<5)
                        dir_6.push_back(packed_index);
                    else
                        dir_6.push_back(temp_index - 1);

                    // DIR_8
                    if ((type & 1<<7) == 1<<7)
                        dir_8.push_back(packed_index);
                    else
                    {
                        if ((type & 1<<6) == 1<<6)
                            dir_8.push_back(temp_index);
                        else
                            dir_8.push_back(temp_index + 1);
                    }
                }

                // DIR_3 DIR_2 DIR_4
                if (i == 0)
                {
                    dir_3.push_back(packed_index);
                    dir_2.push_back(packed_index);
                    dir_4.push_back(packed_index);
                }
                else
                {
                    /// \todo Use Iterator::operator--
                    DenseMatrix<bool>::ElementIterator it(grid.obstacles->element_at(i * grid.obstacles->columns() + j));
                    unsigned long it_index(it.index());
                    unsigned long temp_index(packed_index);
                    while (it.index() != (i - 1) * grid.obstacles->columns() + j)
                    {
                        if (!*it)
                        {
                            --temp_index;
                        }
                        --it_index;
                        it = grid.obstacles->element_at(it_index);
                    }

                    // DIR_3
                    if ((type & (1<<2)) == (1<<2))
                        dir_3.push_back(packed_index);
                    else
                        dir_3.push_back(temp_index);

                    // DIR_2
                    if ((type & 1<<1) == 1<<1)
                        dir_2.push_back(packed_index);
                    else
                        dir_2.push_back(temp_index + 1);

                    // DIR_4
                    if ((type & 1<<3) == 1<<3)
                        dir_4.push_back(packed_index);
                    else
                    {
                        if ((type & 1<<2) == 1<<2)
                            dir_4.push_back(temp_index);
                        else
                            dir_4.push_back(temp_index - 1);
                    }
                }
            }

            static void _expand_direction(DenseVector<unsigned long> * new_dir, DenseVector<unsigned long> * dir,
                    DenseVector<unsigned long> * dir_index)
            {
                if (dir_index->size() > 0)
                {
                    for (unsigned long begin(0), half(0) ; begin < dir_index->size() - 1; begin+=2, ++half)
                    {
                        for (unsigned long i((*dir_index)[begin]), offset(0) ; i < (*dir_index)[begin + 1] ; ++i, ++offset)
                        {
                            (*new_dir)[i] = (*dir)[half] + offset;
                        }
                    }
                }
            }


        public:
            static void cuda_pack(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data)
            {
                // expand types vector
                DenseVector<unsigned long> * new_types(new DenseVector<unsigned long>(data.h->size()));
                for (unsigned long begin(0) ; begin != info.limits->size() - 1 ; ++begin)
                {
                    for (unsigned long i((*info.limits)[begin]) ; i != (*info.limits)[begin + 1] ; ++i)
                    {
                        (*new_types)[i] = (*info.types)[begin];
                    }
                }
                info.cuda_types = new_types;

                // expand direction vectors
                unsigned long inf(-1);
                DenseVector<unsigned long> * new_dir_1(new DenseVector<unsigned long>(data.h->size(), inf));
                DenseVector<unsigned long> * new_dir_2(new DenseVector<unsigned long>(data.h->size(), inf));
                DenseVector<unsigned long> * new_dir_3(new DenseVector<unsigned long>(data.h->size(), inf));
                DenseVector<unsigned long> * new_dir_4(new DenseVector<unsigned long>(data.h->size(), inf));
                DenseVector<unsigned long> * new_dir_5(new DenseVector<unsigned long>(data.h->size(), inf));
                DenseVector<unsigned long> * new_dir_6(new DenseVector<unsigned long>(data.h->size(), inf));
                DenseVector<unsigned long> * new_dir_7(new DenseVector<unsigned long>(data.h->size(), inf));
                DenseVector<unsigned long> * new_dir_8(new DenseVector<unsigned long>(data.h->size(), inf));

                _expand_direction(new_dir_1, info.dir_1, info.dir_index_1);
                _expand_direction(new_dir_2, info.dir_2, info.dir_index_2);
                _expand_direction(new_dir_3, info.dir_3, info.dir_index_3);
                _expand_direction(new_dir_4, info.dir_4, info.dir_index_4);
                _expand_direction(new_dir_5, info.dir_5, info.dir_index_5);
                _expand_direction(new_dir_6, info.dir_6, info.dir_index_6);
                _expand_direction(new_dir_7, info.dir_7, info.dir_index_7);
                _expand_direction(new_dir_8, info.dir_8, info.dir_index_8);

                info.cuda_dir_1 = new_dir_1;
                info.cuda_dir_2 = new_dir_2;
                info.cuda_dir_3 = new_dir_3;
                info.cuda_dir_4 = new_dir_4;
                info.cuda_dir_5 = new_dir_5;
                info.cuda_dir_6 = new_dir_6;
                info.cuda_dir_7 = new_dir_7;
                info.cuda_dir_8 = new_dir_8;
            }

            static void pack(Grid<D2Q9, DT_> & grid, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data, bool alloc_all = true)
            {
                grid.obstacles->lock(lm_read_only);
                grid.h->lock(lm_read_only);
                grid.b->lock(lm_read_only);
                grid.u->lock(lm_read_only);
                grid.v->lock(lm_read_only);

                unsigned long fluid_count(0);
                for(unsigned long i(0); i < grid.obstacles->rows(); ++i)
                {
                    for(unsigned long j(0); j < grid.obstacles->columns(); ++j)
                    {
                        ///Zero is fluid, One means obstacle
                        if(!(*grid.obstacles)(i,j))
                        {
                            if( i > 0
                                    && i + 1 < grid.obstacles->rows()
                                    && j > 0
                                    && j + 1 < grid.obstacles->columns()
                                    && (*grid.obstacles)(i - 1, j)
                                    && (*grid.obstacles)(i, j - 1)
                                    && (*grid.obstacles)(i - 1, j - 1)
                                    && (*grid.obstacles)(i + 1, j)
                                    && (*grid.obstacles)(i, j + 1)
                                    && (*grid.obstacles)(i + 1, j + 1)
                                    && (*grid.obstacles)(i + 1, j - 1)
                                    && (*grid.obstacles)(i - 1, j + 1))
                            {
                                (*grid.obstacles)(i,j) = true;
                            }
                            else
                            {
                                ++fluid_count;
                            }
                        }
                    }
                }
                grid.h_index = new DenseMatrix<unsigned long>(grid.h->rows(), grid.h->columns(), grid.h->size());
                data.h = new DenseVector<DT_>(fluid_count);
                data.b = new DenseVector<DT_>(fluid_count);
                data.u = new DenseVector<DT_>(fluid_count);
                data.v = new DenseVector<DT_>(fluid_count);

                if (alloc_all)
                {
                    data.temp = new DenseVector<DT_>(fluid_count);

                    data.f_0 = new DenseVector<DT_>(fluid_count, DT_(0));
                    data.f_1 = new DenseVector<DT_>(fluid_count, DT_(0));
                    data.f_2 = new DenseVector<DT_>(fluid_count, DT_(0));
                    data.f_3 = new DenseVector<DT_>(fluid_count, DT_(0));
                    data.f_4 = new DenseVector<DT_>(fluid_count, DT_(0));
                    data.f_5 = new DenseVector<DT_>(fluid_count, DT_(0));
                    data.f_6 = new DenseVector<DT_>(fluid_count, DT_(0));
                    data.f_7 = new DenseVector<DT_>(fluid_count, DT_(0));
                    data.f_8 = new DenseVector<DT_>(fluid_count, DT_(0));

                    data.f_eq_0 = new DenseVector<DT_>(fluid_count, DT_(0));
                    data.f_eq_1 = new DenseVector<DT_>(fluid_count, DT_(0));
                    data.f_eq_2 = new DenseVector<DT_>(fluid_count, DT_(0));
                    data.f_eq_3 = new DenseVector<DT_>(fluid_count, DT_(0));
                    data.f_eq_4 = new DenseVector<DT_>(fluid_count, DT_(0));
                    data.f_eq_5 = new DenseVector<DT_>(fluid_count, DT_(0));
                    data.f_eq_6 = new DenseVector<DT_>(fluid_count, DT_(0));
                    data.f_eq_7 = new DenseVector<DT_>(fluid_count, DT_(0));
                    data.f_eq_8 = new DenseVector<DT_>(fluid_count, DT_(0));

                    data.f_temp_0 = new DenseVector<DT_>(fluid_count, DT_(0));
                    data.f_temp_1 = new DenseVector<DT_>(fluid_count, DT_(0));
                    data.f_temp_2 = new DenseVector<DT_>(fluid_count, DT_(0));
                    data.f_temp_3 = new DenseVector<DT_>(fluid_count, DT_(0));
                    data.f_temp_4 = new DenseVector<DT_>(fluid_count, DT_(0));
                    data.f_temp_5 = new DenseVector<DT_>(fluid_count, DT_(0));
                    data.f_temp_6 = new DenseVector<DT_>(fluid_count, DT_(0));
                    data.f_temp_7 = new DenseVector<DT_>(fluid_count, DT_(0));
                    data.f_temp_8 = new DenseVector<DT_>(fluid_count, DT_(0));
                    data.distribution_x = new DenseVector<DT_>(9ul, DT_(0));
                    data.distribution_y = new DenseVector<DT_>(9ul, DT_(0));
                }

                std::vector<unsigned long> temp_limits;
                std::vector<unsigned long> temp_types;
                std::vector<unsigned long> dir_1;
                std::vector<unsigned long> dir_2;
                std::vector<unsigned long> dir_3;
                std::vector<unsigned long> dir_4;
                std::vector<unsigned long> dir_5;
                std::vector<unsigned long> dir_6;
                std::vector<unsigned long> dir_7;
                std::vector<unsigned long> dir_8;
                std::vector<unsigned long> dir_post_1;
                std::vector<unsigned long> dir_post_2;
                std::vector<unsigned long> dir_post_3;
                std::vector<unsigned long> dir_post_4;
                std::vector<unsigned long> dir_post_5;
                std::vector<unsigned long> dir_post_6;
                std::vector<unsigned long> dir_post_7;
                std::vector<unsigned long> dir_post_8;
                std::vector<unsigned long> dir_index_1;
                std::vector<unsigned long> dir_index_2;
                std::vector<unsigned long> dir_index_3;
                std::vector<unsigned long> dir_index_4;
                std::vector<unsigned long> dir_index_5;
                std::vector<unsigned long> dir_index_6;
                std::vector<unsigned long> dir_index_7;
                std::vector<unsigned long> dir_index_8;

                unsigned long packed_index(0);

                for(unsigned long i(0); i < grid.obstacles->rows(); ++i)
                {
                    for(unsigned long j(0); j < grid.obstacles->columns(); ++j)
                    {
                        if((*grid.obstacles)(i, j))
                        {
                            // do nothig, ignore obstacle cell
                        }
                        else
                        {
                            if(j > 0)
                            {
                                if (_element_type(i, j, grid) != _element_type(i, j - 1, grid))
                                {
                                    // common case
                                    // insert current cell
                                    temp_limits.push_back(packed_index);
                                    temp_types.push_back(_element_type(i, j, grid));
                                    _element_direction(packed_index, i, j, grid, dir_1, dir_2, dir_3, dir_4, dir_5, dir_6, dir_7, dir_8);
                                }
                            }
                            else
                            {
                                // leftmost boundary cells
                                temp_limits.push_back(packed_index);
                                temp_types.push_back(_element_type(i, j, grid));
                                _element_direction(packed_index, i, j, grid, dir_1, dir_2, dir_3, dir_4, dir_5, dir_6, dir_7, dir_8);
                            }
                            (*grid.h_index)(i, j) = packed_index;
                            ++packed_index;
                        }
                    }
                }
                temp_limits.push_back(packed_index);
                temp_types.push_back(0);
                dir_1.push_back(packed_index - 1);
                dir_2.push_back(packed_index - 1);
                dir_3.push_back(packed_index - 1);
                dir_4.push_back(packed_index - 1);
                dir_5.push_back(packed_index - 1);
                dir_6.push_back(packed_index - 1);
                dir_7.push_back(packed_index - 1);
                dir_8.push_back(packed_index - 1);

                _dir_post_pack(dir_1, dir_post_1, dir_index_1, temp_limits);
                _dir_post_pack(dir_2, dir_post_2, dir_index_2, temp_limits);
                _dir_post_pack(dir_3, dir_post_3, dir_index_3, temp_limits);
                _dir_post_pack(dir_4, dir_post_4, dir_index_4, temp_limits);
                _dir_post_pack(dir_5, dir_post_5, dir_index_5, temp_limits);
                _dir_post_pack(dir_6, dir_post_6, dir_index_6, temp_limits);
                _dir_post_pack(dir_7, dir_post_7, dir_index_7, temp_limits);
                _dir_post_pack(dir_8, dir_post_8, dir_index_8, temp_limits);

                info.limits = new DenseVector<unsigned long>(temp_limits.size());
                info.types = new DenseVector<unsigned long>(temp_limits.size());
                info.dir_1 = new DenseVector<unsigned long>(dir_post_1.size());
                info.dir_2 = new DenseVector<unsigned long>(dir_post_2.size());
                info.dir_3 = new DenseVector<unsigned long>(dir_post_3.size());
                info.dir_4 = new DenseVector<unsigned long>(dir_post_4.size());
                info.dir_5 = new DenseVector<unsigned long>(dir_post_5.size());
                info.dir_6 = new DenseVector<unsigned long>(dir_post_6.size());
                info.dir_7 = new DenseVector<unsigned long>(dir_post_7.size());
                info.dir_8 = new DenseVector<unsigned long>(dir_post_8.size());
                info.dir_index_1 = new DenseVector<unsigned long>(dir_index_1.size());
                info.dir_index_2 = new DenseVector<unsigned long>(dir_index_2.size());
                info.dir_index_3 = new DenseVector<unsigned long>(dir_index_3.size());
                info.dir_index_4 = new DenseVector<unsigned long>(dir_index_4.size());
                info.dir_index_5 = new DenseVector<unsigned long>(dir_index_5.size());
                info.dir_index_6 = new DenseVector<unsigned long>(dir_index_6.size());
                info.dir_index_7 = new DenseVector<unsigned long>(dir_index_7.size());
                info.dir_index_8 = new DenseVector<unsigned long>(dir_index_8.size());

                unsigned long index2(0);
                for (std::vector<unsigned long>::iterator i(temp_limits.begin()), j(temp_types.begin()) ; i != temp_limits.end() ; ++i, ++j)
                {
                    (*info.limits)[index2] = *i;
                    (*info.types)[index2] = *j;
                    ++index2;
                }
                for (unsigned long i(0) ; i < dir_post_1.size() ; ++i)
                {
                    (*info.dir_1)[i] = dir_post_1[i];
                }
                for (unsigned long i(0) ; i < dir_post_2.size() ; ++i)
                {
                    (*info.dir_2)[i] = dir_post_2[i];
                }
                for (unsigned long i(0) ; i < dir_post_3.size() ; ++i)
                {
                    (*info.dir_3)[i] = dir_post_3[i];
                }
                for (unsigned long i(0) ; i < dir_post_4.size() ; ++i)
                {
                    (*info.dir_4)[i] = dir_post_4[i];
                }
                for (unsigned long i(0) ; i < dir_post_5.size() ; ++i)
                {
                    (*info.dir_5)[i] = dir_post_5[i];
                }
                for (unsigned long i(0) ; i < dir_post_6.size() ; ++i)
                {
                    (*info.dir_6)[i] = dir_post_6[i];
                }
                for (unsigned long i(0) ; i < dir_post_7.size() ; ++i)
                {
                    (*info.dir_7)[i] = dir_post_7[i];
                }
                for (unsigned long i(0) ; i < dir_post_8.size() ; ++i)
                {
                    (*info.dir_8)[i] = dir_post_8[i];
                }
                for (unsigned long i(0) ; i < dir_index_1.size() ; ++i)
                {
                    (*info.dir_index_1)[i] = dir_index_1[i];
                }
                for (unsigned long i(0) ; i < dir_index_2.size() ; ++i)
                {
                    (*info.dir_index_2)[i] = dir_index_2[i];
                }
                for (unsigned long i(0) ; i < dir_index_3.size() ; ++i)
                {
                    (*info.dir_index_3)[i] = dir_index_3[i];
                }
                for (unsigned long i(0) ; i < dir_index_4.size() ; ++i)
                {
                    (*info.dir_index_4)[i] = dir_index_4[i];
                }
                for (unsigned long i(0) ; i < dir_index_5.size() ; ++i)
                {
                    (*info.dir_index_5)[i] = dir_index_5[i];
                }
                for (unsigned long i(0) ; i < dir_index_6.size() ; ++i)
                {
                    (*info.dir_index_6)[i] = dir_index_6[i];
                }
                for (unsigned long i(0) ; i < dir_index_7.size() ; ++i)
                {
                    (*info.dir_index_7)[i] = dir_index_7[i];
                }
                for (unsigned long i(0) ; i < dir_index_8.size() ; ++i)
                {
                    (*info.dir_index_8)[i] = dir_index_8[i];
                }

                index2 = 0;
                for (DenseVector<unsigned long>::ConstElementIterator begin(info.limits->begin_elements()),end(info.limits->element_at(1))
                        ; end != info.limits->end_elements() ; ++begin, ++end)
                {
                    for (unsigned long i(*begin) ; i < *end ; ++i)
                    {
                        while ((*grid.obstacles)(index2 / grid.obstacles->columns(), index2 % grid.obstacles->columns()))
                        {
                            ++index2;
                        }
                        (*data.h)[i] = (*grid.h)(index2 / grid.obstacles->columns(), index2 % grid.obstacles->columns());
                        (*data.b)[i] = (*grid.b)(index2 / grid.obstacles->columns(), index2 % grid.obstacles->columns());
                        (*data.u)[i] = (*grid.u)(index2 / grid.obstacles->columns(), index2 % grid.obstacles->columns());
                        (*data.v)[i] = (*grid.v)(index2 / grid.obstacles->columns(), index2 % grid.obstacles->columns());
                        ++index2;
                    }
                }

                grid.obstacles->unlock(lm_read_only);
                grid.h->unlock(lm_read_only);
                grid.b->unlock(lm_read_only);
                grid.u->unlock(lm_read_only);
                grid.v->unlock(lm_read_only);
            }

            ///Generalized unpacking and extraction
            static void deflate(Grid<D2Q9, DT_> & grid, DenseVector<DT_> * from, DenseMatrix<DT_> * to)
            {
                grid.obstacles->lock(lm_read_only);
                to->lock(lm_write_only);
                from->lock(lm_read_only);
                unsigned long packed_index(0);

                for(unsigned long i(0); i < grid.obstacles->rows(); ++i)
                {
                    for(unsigned long j(0); j < grid.obstacles->columns(); ++j)
                    {
                        if((*grid.obstacles)(i, j))
                        {
                            (*to)(i, j) = DT_(0);
                        }
                        else
                        {
                            (*to)(i, j) = (*from)[packed_index];
                            ++packed_index;
                        }
                    }
                }
                grid.obstacles->unlock(lm_read_only);
                to->unlock(lm_write_only);
                from->unlock(lm_read_only);
            }
            ///END Generalized unpacking and extraction

            static void unpack(Grid<D2Q9, DT_> & grid, HONEI_UNUSED PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data)
            {
                grid.obstacles->lock(lm_read_only);
                grid.h->lock(lm_write_only);
                data.h->lock(lm_read_only);
                unsigned long packed_index(0);

                for(unsigned long i(0); i < grid.obstacles->rows(); ++i)
                {
                    for(unsigned long j(0); j < grid.obstacles->columns(); ++j)
                    {
                        if((*grid.obstacles)(i, j))
                        {
                            (*grid.h)(i, j) = DT_(0);
                        }
                        else
                        {
                            (*grid.h)(i, j) = (*data.h)[packed_index];
                            ++packed_index;
                        }
                    }
                }
                grid.obstacles->unlock(lm_read_only);
                grid.h->unlock(lm_write_only);
                data.h->unlock(lm_read_only);
            }
            static void unpack_u(Grid<D2Q9, DT_> & grid, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data)
            {
                grid.obstacles->lock(lm_read_only);
                grid.u->lock(lm_write_only);
                data.u->lock(lm_read_only);
                unsigned long packed_index(0);

                for(unsigned long i(0); i < grid.obstacles->rows(); ++i)
                {
                    for(unsigned long j(0); j < grid.obstacles->columns(); ++j)
                    {
                        if((*grid.obstacles)(i, j))
                        {
                            (*grid.u)(i, j) = DT_(0);
                        }
                        else
                        {
                            (*grid.u)(i, j) = (*data.u)[packed_index];
                            ++packed_index;
                        }
                    }
                }
                grid.obstacles->unlock(lm_read_only);
                grid.u->unlock(lm_write_only);
                data.u->unlock(lm_read_only);
            }


            static unsigned long h_index(Grid<D2Q9, DT_> & grid, unsigned long i, unsigned long j)
            {
                return (*grid.h_index)(i, j);
            }

            static DenseMatrix<DT_> extract_ftemp2(Grid<D2Q9, DT_> & grid, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data)
            {
                unsigned long packed_index(0);
                DenseMatrix<DT_> result(grid.obstacles->rows(), grid.obstacles->columns());

                for(unsigned long i(0); i < grid.obstacles->rows(); ++i)
                {
                    for(unsigned long j(0); j < grid.obstacles->columns(); ++j)
                    {
                        if((*grid.obstacles)(i, j))
                        {
                            result(i, j) = DT_(0);
                        }
                        else
                        {
                            result(i, j) = (*data.f_temp_2)[packed_index];
                            ++packed_index;
                        }
                    }
                }
                return result;
            }

            static DenseMatrix<DT_> extract_f2(Grid<D2Q9, DT_> & grid, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data)
            {
                unsigned long packed_index(0);
                DenseMatrix<DT_> result(grid.obstacles->rows(), grid.obstacles->columns());

                for(unsigned long i(0); i < grid.obstacles->rows(); ++i)
                {
                    for(unsigned long j(0); j < grid.obstacles->columns(); ++j)
                    {
                        if((*grid.obstacles)(i, j))
                        {
                            result(i, j) = DT_(0);
                        }
                        else
                        {
                            result(i, j) = (*data.f_2)[packed_index];
                            ++packed_index;
                        }
                    }
                }
                return result;
            }

            static DenseMatrix<DT_> extract_f6(Grid<D2Q9, DT_> & grid, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data)
            {
                unsigned long packed_index(0);
                DenseMatrix<DT_> result(grid.obstacles->rows(), grid.obstacles->columns());

                for(unsigned long i(0); i < grid.obstacles->rows(); ++i)
                {
                    for(unsigned long j(0); j < grid.obstacles->columns(); ++j)
                    {
                        if((*grid.obstacles)(i, j))
                        {
                            result(i, j) = DT_(0);
                        }
                        else
                        {
                            result(i, j) = (*data.f_6)[packed_index];
                            ++packed_index;
                        }
                    }
                }
                return result;
            }

            static DenseMatrix<DT_> extract_ftemp6(Grid<D2Q9, DT_> & grid, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data)
            {
                unsigned long packed_index(0);
                DenseMatrix<DT_> result(grid.obstacles->rows(), grid.obstacles->columns());

                for(unsigned long i(0); i < grid.obstacles->rows(); ++i)
                {
                    for(unsigned long j(0); j < grid.obstacles->columns(); ++j)
                    {
                        if((*grid.obstacles)(i, j))
                        {
                            result(i, j) = DT_(0);
                        }
                        else
                        {
                            result(i, j) = (*data.f_temp_6)[packed_index];
                            ++packed_index;
                        }
                    }
                }
                return result;
            }
    };

    ///FSI types:
    template <typename LatticeType_, typename BoundaryType_, typename DT_> struct GridPackerFSI
    {
    };

    template <typename DT_> struct GridPackerFSI<D2Q9, lbm_boundary_types::NOSLIP, DT_>
    {
        public:
            ///Generalized unpacking and extraction
            static void deflate(Grid<D2Q9, DT_> & grid, HONEI_UNUSED PackedSolidData<D2Q9, DT_> & data, DenseVector<bool> * from, DenseMatrix<bool> * to)
            {
                grid.obstacles->lock(lm_read_only);
                to->lock(lm_write_only);
                from->lock(lm_read_only);
                unsigned long packed_index(0);

                for(unsigned long i(0); i < grid.obstacles->rows(); ++i)
                {
                    for(unsigned long j(0); j < grid.obstacles->columns(); ++j)
                    {
                        if((*grid.obstacles)(i, j))
                        {
                            (*to)(i, j) = DT_(0);
                        }
                        else
                        {
                            (*to)(i, j) = (*from)[packed_index];
                            ++packed_index;
                        }
                    }
                }
                grid.obstacles->unlock(lm_read_only);
                to->unlock(lm_write_only);
                from->unlock(lm_read_only);
            }
            static void deflate(Grid<D2Q9, DT_> & grid, HONEI_UNUSED PackedSolidData<D2Q9, DT_> & data, DenseVector<DT_> * from, DenseMatrix<DT_> * to)
            {
                grid.obstacles->lock(lm_read_only);
                to->lock(lm_write_only);
                from->lock(lm_read_only);
                unsigned long packed_index(0);

                for(unsigned long i(0); i < grid.obstacles->rows(); ++i)
                {
                    for(unsigned long j(0); j < grid.obstacles->columns(); ++j)
                    {
                        if((*grid.obstacles)(i, j))
                        {
                            (*to)(i, j) = DT_(0);
                        }
                        else
                        {
                            (*to)(i, j) = (*from)[packed_index];
                            ++packed_index;
                        }
                    }
                }
                grid.obstacles->unlock(lm_read_only);
                to->unlock(lm_write_only);
                from->unlock(lm_read_only);
            }
            ///END Generalized unpacking and extraction

            static void allocate(PackedGridData<D2Q9, DT_> & data, PackedSolidData<D2Q9, DT_> & solids)
            {
                solids.boundary_flags = new DenseVector<bool>(data.h->size());
                solids.line_flags = new DenseVector<bool>(data.h->size());
                solids.solid_flags = new DenseVector<bool>(data.h->size());
                solids.solid_old_flags = new DenseVector<bool>(data.h->size());
                solids.solid_to_fluid_flags = new DenseVector<bool>(data.h->size());
                solids.stationary_flags = new DenseVector<bool>(data.h->size());


                solids.f_mea_1 = new DenseVector<DT_>(data.h->size());
                solids.f_mea_2 = new DenseVector<DT_>(data.h->size());
                solids.f_mea_3 = new DenseVector<DT_>(data.h->size());
                solids.f_mea_4 = new DenseVector<DT_>(data.h->size());
                solids.f_mea_5 = new DenseVector<DT_>(data.h->size());
                solids.f_mea_6 = new DenseVector<DT_>(data.h->size());
                solids.f_mea_7 = new DenseVector<DT_>(data.h->size());
                solids.f_mea_8 = new DenseVector<DT_>(data.h->size());
            }

            static void pack(Grid<D2Q9, DT_> & grid,
                    HONEI_UNUSED PackedGridData<D2Q9, DT_> & data,
                    PackedSolidData<D2Q9, DT_> & solids,
                    DenseMatrix<bool> & lines,
                    DenseMatrix<bool> & boundaries,
                    DenseMatrix<bool> & solid,
                    DenseMatrix<bool> & solid_to_fluid,
                    DenseMatrix<bool> & stationary)
            {

                solids.boundary_flags->lock(lm_write_only);
                solids.line_flags->lock(lm_write_only);
                solids.solid_flags->lock(lm_read_and_write);
                solids.solid_old_flags->lock(lm_write_only);
                solids.solid_to_fluid_flags->lock(lm_write_only);
                solids.stationary_flags->lock(lm_write_only);
                solids.f_mea_1->lock(lm_write_only);
                solids.f_mea_2->lock(lm_write_only);
                solids.f_mea_3->lock(lm_write_only);
                solids.f_mea_4->lock(lm_write_only);
                solids.f_mea_5->lock(lm_write_only);
                solids.f_mea_6->lock(lm_write_only);
                solids.f_mea_7->lock(lm_write_only);
                solids.f_mea_8->lock(lm_write_only);
                for (unsigned long i(0) ; i < lines.rows() ; ++i)
                {
                    for (unsigned long j(0) ; j < lines.columns() ; ++j)
                    {
                        unsigned long packed_index(GridPacker<D2Q9, lbm_boundary_types::NOSLIP, DT_>::h_index(grid, i, j));
                        (*solids.boundary_flags)[packed_index] = boundaries[i][j] ? true : false;
                        (*solids.line_flags)[packed_index] = lines[i][j] ? true : false;
                        (*solids.solid_flags)[packed_index] = solid[i][j] ? true : false;
                        (*solids.solid_old_flags)[packed_index] = (*solids.solid_flags)[packed_index]; //equals solid at start
                        (*solids.solid_to_fluid_flags)[packed_index] = solid_to_fluid[i][j] ? true : false;
                        (*solids.stationary_flags)[packed_index] = stationary[i][j] ? true : false;

                        (*solids.f_mea_1)[packed_index] = DT_(0);
                        (*solids.f_mea_2)[packed_index] = DT_(0);
                        (*solids.f_mea_3)[packed_index] = DT_(0);
                        (*solids.f_mea_4)[packed_index] = DT_(0);
                        (*solids.f_mea_5)[packed_index] = DT_(0);
                        (*solids.f_mea_6)[packed_index] = DT_(0);
                        (*solids.f_mea_7)[packed_index] = DT_(0);
                        (*solids.f_mea_8)[packed_index] = DT_(0);
                    }
                }
                solids.f_mea_1->unlock(lm_write_only);
                solids.f_mea_2->unlock(lm_write_only);
                solids.f_mea_3->unlock(lm_write_only);
                solids.f_mea_4->unlock(lm_write_only);
                solids.f_mea_5->unlock(lm_write_only);
                solids.f_mea_6->unlock(lm_write_only);
                solids.f_mea_7->unlock(lm_write_only);
                solids.f_mea_8->unlock(lm_write_only);
                solids.boundary_flags->unlock(lm_write_only);
                solids.line_flags->unlock(lm_write_only);
                solids.solid_flags->unlock(lm_read_and_write);
                solids.solid_old_flags->unlock(lm_write_only);
                solids.solid_to_fluid_flags->unlock(lm_write_only);
                solids.stationary_flags->unlock(lm_write_only);
            }
    };
}

#endif
