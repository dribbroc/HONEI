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


#ifndef LBM_GUARD_GRID_PACKER_HH
#define LBM_GUARD_GRID_PACKER_HH 1

#include <honei/lbm/grid.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/dense_matrix.hh>
#include <vector>
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
            static BoundaryTypeD2Q9 _element_type(unsigned long i, unsigned long j, Grid<D2Q9, DT_> & grid)
            {
                bool north(false);
                bool south(false);
                bool west(false);
                bool east(false);

                if(i == 0)
                    north = true;

                if(i == (*grid.obstacles).rows() - 1)
                    south = true;

                if(j == 0)
                    west = true;

                if(j == (*grid.obstacles).columns() - 1)
                    east = true;

                if( i > 0 && (*grid.obstacles)(i - 1, j))
                    north = true;

                if( i <  (*grid.obstacles).rows() - 1 && (*grid.obstacles)(i + 1, j))
                    south = true;

                if( j > 0 && (*grid.obstacles)(i, j - 1))
                    west = true;

                if( j <  (*grid.obstacles).columns() - 1 && (*grid.obstacles)(i, j + 1))
                    east = true;

                if((*grid.obstacles)(i, j) || north && south && west && east)
                    return bt_obstacle;

                if(!north && !south && !west && !east)
                    return bt_none;

                if(north && south && west && !east)
                    return bt_S_W_N;

                if(north && south && !west && east)
                    return bt_N_E_S;

                if(north && !south && west && east)
                    return bt_W_N_E;

                if(!north && south && west && east)
                    return bt_E_S_W;

                if(!north && !south && west && east)
                    return bt_E_W;

                if(north && south && !west && !east)
                    return bt_N_S;

                if(!north && south && !west && east)
                    return bt_E_S;

                if(!north && south && west && !east)
                    return bt_S_W;

                if(north && !south && west && !east)
                    return bt_W_N;

                if(north && !south && !west && east)
                    return bt_N_E;

                if(north && !south && !west && !east)
                    return bt_N;

                if(!north && south && !west && !east)
                    return bt_S;

                if(!north && !south && west && !east)
                    return bt_W;

                if(!north && !south && !west && east)
                    return bt_E;
            }

            static void _element_direction(unsigned long packed_index, unsigned long i, unsigned long j, Grid<D2Q9, DT_> & grid,
                    std::vector<unsigned long> & dir_0,
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
                BoundaryTypeD2Q9 boundary(_element_type(i, j, grid));

                // DIR_0
                dir_0.push_back(packed_index);

                // DIR_1
                if (boundary == bt_E || boundary == bt_N_E || boundary == bt_W_N_E || boundary == bt_E_W || boundary == bt_E_S ||
                        boundary == bt_E_S_W)
                    dir_1.push_back(packed_index);
                else
                    dir_1.push_back(packed_index + 1);

                // DIR_5
                if (boundary == bt_W || boundary == bt_W_N || boundary == bt_W_N_E || boundary == bt_E_W || boundary == bt_S_W ||
                        boundary == bt_E_S_W)
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
                    MutableMatrix<bool>::ElementIterator it(grid.obstacles->element_at(i * grid.obstacles->columns() + j));
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
                    if (boundary == bt_S || boundary == bt_S_W || boundary == bt_E_S || boundary == bt_N_S || boundary == bt_S_W_N ||
                            boundary == bt_E_S_W)
                        dir_7.push_back(packed_index);
                    else
                        dir_7.push_back(temp_index);

                    // DIR_6
                    if (j > 0 && i < grid.obstacles->rows() - 1)
                    {
                        if ((*grid.obstacles)(i + 1, j - 1))
                            dir_6.push_back(packed_index);
                        else
                            if (boundary == bt_S_W || boundary == bt_S_W_N || boundary == bt_E_S_W)
                                dir_6.push_back(packed_index);
                            else
                                dir_6.push_back(temp_index - 1);
                    }
                    else
                        dir_6.push_back(packed_index);

                    // DIR_8
                    if (j < grid.obstacles->columns() - 1 && i < grid.obstacles->rows() -1)
                    {
                        if ((*grid.obstacles)(i + 1, j + 1))
                            dir_8.push_back(packed_index);
                        else
                            if (boundary == bt_E_S || boundary == bt_N_E_S || boundary == bt_E_S_W)
                                dir_8.push_back(packed_index);
                            else
                                dir_8.push_back(temp_index + 1);
                    }
                    else
                        dir_8.push_back(packed_index);
                }

                // DIR_3 DIR_2 DIR_4
                if (i == 0)
                {
                    dir_3.push_back(0);
                    dir_2.push_back(0);
                    dir_4.push_back(0);
                }
                else
                {
                    /// \todo Use Iterator::operator--
                    MutableMatrix<bool>::ElementIterator it(grid.obstacles->element_at(i * grid.obstacles->columns() + j));
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
                    if (boundary == bt_N || boundary == bt_W_N || boundary == bt_N_E || boundary == bt_N_S || boundary == bt_S_W_N ||
                            boundary == bt_W_N_E)
                        dir_3.push_back(packed_index);
                    else
                        dir_3.push_back(temp_index);

                    // DIR_2
                    if (j < grid.obstacles->columns() - 1 && i > 0)
                    {
                        if ((*grid.obstacles)(i -1, j + 1))
                            dir_2.push_back(packed_index);
                        else
                            if (boundary == bt_N_E || boundary == bt_N_E_S || boundary == bt_W_N_E)
                                dir_2.push_back(packed_index);
                            else
                                dir_2.push_back(temp_index + 1);
                    }
                    else
                        dir_2.push_back(packed_index);

                    // DIR_4
                    if (j > 0 && i > 0)
                    {
                        if ((*grid.obstacles)(i -1, j - 1))
                            dir_4.push_back(packed_index);
                        else
                            if (boundary == bt_W_N || boundary == bt_S_W_N || boundary == bt_W_N_E)
                                dir_4.push_back(packed_index);
                            else
                                dir_4.push_back(temp_index -1);
                    }
                    else
                        dir_4.push_back(packed_index);
                }
            }


        public:
            static void pack(Grid<D2Q9, DT_> & grid, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data)
            {
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
                std::cout<<"Fluid count: "<<fluid_count<<std::endl;
                data.h = new DenseVector<DT_>(fluid_count);
                data.u = new DenseVector<DT_>(fluid_count);
                data.v = new DenseVector<DT_>(fluid_count);

                data.f_0 = new DenseVector<DT_>(fluid_count);
                data.f_1 = new DenseVector<DT_>(fluid_count);
                data.f_2 = new DenseVector<DT_>(fluid_count);
                data.f_3 = new DenseVector<DT_>(fluid_count);
                data.f_4 = new DenseVector<DT_>(fluid_count);
                data.f_5 = new DenseVector<DT_>(fluid_count);
                data.f_6 = new DenseVector<DT_>(fluid_count);
                data.f_7 = new DenseVector<DT_>(fluid_count);
                data.f_8 = new DenseVector<DT_>(fluid_count);

                data.f_eq_0 = new DenseVector<DT_>(fluid_count);
                data.f_eq_1 = new DenseVector<DT_>(fluid_count);
                data.f_eq_2 = new DenseVector<DT_>(fluid_count);
                data.f_eq_3 = new DenseVector<DT_>(fluid_count);
                data.f_eq_4 = new DenseVector<DT_>(fluid_count);
                data.f_eq_5 = new DenseVector<DT_>(fluid_count);
                data.f_eq_6 = new DenseVector<DT_>(fluid_count);
                data.f_eq_7 = new DenseVector<DT_>(fluid_count);
                data.f_eq_8 = new DenseVector<DT_>(fluid_count);

                data.f_temp_0 = new DenseVector<DT_>(fluid_count);
                data.f_temp_1 = new DenseVector<DT_>(fluid_count);
                data.f_temp_2 = new DenseVector<DT_>(fluid_count);
                data.f_temp_3 = new DenseVector<DT_>(fluid_count);
                data.f_temp_4 = new DenseVector<DT_>(fluid_count);
                data.f_temp_5 = new DenseVector<DT_>(fluid_count);
                data.f_temp_6 = new DenseVector<DT_>(fluid_count);
                data.f_temp_7 = new DenseVector<DT_>(fluid_count);
                data.f_temp_8 = new DenseVector<DT_>(fluid_count);

                std::vector<unsigned long> temp_limits;
                std::vector<unsigned long> dir_0;
                std::vector<unsigned long> dir_1;
                std::vector<unsigned long> dir_2;
                std::vector<unsigned long> dir_3;
                std::vector<unsigned long> dir_4;
                std::vector<unsigned long> dir_5;
                std::vector<unsigned long> dir_6;
                std::vector<unsigned long> dir_7;
                std::vector<unsigned long> dir_8;

                unsigned long packed_index(0);
                bool dirty_cell(false);

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
                                    // search for obstacles in diagonal directions DIR_4 and DIR_6
                                    /// \todo 2 directions too much in condition
                                    if(
                                            (i > 0 &&  (*grid.obstacles)(i - 1, j - 1)) ||
                                            (i > 0 &&  j < grid.obstacles->columns() - 1 && (*grid.obstacles)(i - 1, j + 1)) ||
                                            (i < grid.obstacles->rows() - 1 &&  (*grid.obstacles)(i + 1, j - 1)) ||
                                            (i < grid.obstacles->rows() - 1 &&  j < grid.obstacles->columns() - 1 && (*grid.obstacles)(i + 1, j + 1))
                                       )
                                    {
                                        dirty_cell = true;
                                    }

                                    // insert current cell
                                    temp_limits.push_back(packed_index);
                                    _element_direction(packed_index, i, j, grid, dir_0, dir_1, dir_2, dir_3, dir_4, dir_5, dir_6, dir_7, dir_8);
                                    std::cout<<packed_index<<" "<<_element_type(i, j, grid)<<std::endl;

                                }
                                // search for obstacles in diagonal directions DIR_2 and DIR_8
                                /// \todo 2 directions too much in condition
                                else if (
                                            (i > 0 &&  (*grid.obstacles)(i - 1, j - 1)) ||
                                            (i > 0 &&  j < grid.obstacles->columns() - 1 && (*grid.obstacles)(i - 1, j + 1)) ||
                                            (i < grid.obstacles->rows() - 1 &&  (*grid.obstacles)(i + 1, j - 1)) ||
                                            (i < grid.obstacles->rows() - 1 &&  j < grid.obstacles->columns() - 1 && (*grid.obstacles)(i + 1, j + 1))
                                        )
                                {
                                    // insert current cell
                                    temp_limits.push_back(packed_index);
                                    _element_direction(packed_index, i, j, grid, dir_0, dir_1, dir_2, dir_3, dir_4, dir_5, dir_6, dir_7, dir_8);
                                    std::cout<<packed_index<<" "<<_element_type(i, j, grid)<<std::endl;
                                    dirty_cell = false;
                                }
                                else if (dirty_cell)
                                {
                                    // insert current cell
                                    temp_limits.push_back(packed_index);
                                    _element_direction(packed_index, i, j, grid, dir_0, dir_1, dir_2, dir_3, dir_4, dir_5, dir_6, dir_7, dir_8);
                                    std::cout<<packed_index<<" "<<_element_type(i, j, grid)<<std::endl;
                                    dirty_cell = false;
                                }
                            }
                            else
                            {
                                // leftmost boundary cells
                                temp_limits.push_back(packed_index);
                                _element_direction(packed_index, i, j, grid, dir_0, dir_1, dir_2, dir_3, dir_4, dir_5, dir_6, dir_7, dir_8);
                                std::cout<<packed_index<<" "<<_element_type(i, j, grid)<<std::endl;
                                dirty_cell = false;
                            }
                            ++packed_index;
                        }
                    }
                }
                temp_limits.push_back(packed_index);

                info.limits = new DenseVector<unsigned long>(temp_limits.size());
                info.dir_0 = new DenseVector<unsigned long>(dir_0.size());
                info.dir_1 = new DenseVector<unsigned long>(dir_0.size());
                info.dir_2 = new DenseVector<unsigned long>(dir_0.size());
                info.dir_3 = new DenseVector<unsigned long>(dir_0.size());
                info.dir_4 = new DenseVector<unsigned long>(dir_0.size());
                info.dir_5 = new DenseVector<unsigned long>(dir_0.size());
                info.dir_6 = new DenseVector<unsigned long>(dir_0.size());
                info.dir_7 = new DenseVector<unsigned long>(dir_0.size());
                info.dir_8 = new DenseVector<unsigned long>(dir_0.size());

                unsigned long index2(0);
                for (std::vector<unsigned long>::iterator i(temp_limits.begin()) ; i != temp_limits.end() ; ++i)//, ++j)
                {
                    (*info.limits)[index2] = *i;
                    ++index2;
                }
                for (unsigned long i(0) ; i < dir_0.size() ; ++i)
                {
                    (*info.dir_0)[i] = dir_0[i];
                    (*info.dir_1)[i] = dir_1[i];
                    (*info.dir_2)[i] = dir_2[i];
                    (*info.dir_3)[i] = dir_3[i];
                    (*info.dir_4)[i] = dir_4[i];
                    (*info.dir_5)[i] = dir_5[i];
                    (*info.dir_6)[i] = dir_6[i];
                    (*info.dir_7)[i] = dir_7[i];
                    (*info.dir_8)[i] = dir_8[i];
                }

                index2 = 0;
                for (DenseVector<unsigned long>::ConstElementIterator begin(info.limits->begin_elements()),end(info.limits->element_at(1))
                        ; end != info.limits->end_elements() ; ++begin, ++end)
                {
                    for (unsigned long i(*begin) ; i < *end ; ++i)
                    {
                        while ((*grid.obstacles)(index2))
                        {
                            ++index2;
                        }
                        (*data.h)[i] = (*grid.h)(index2);
                        (*data.u)[i] = (*grid.u)(index2);
                        (*data.v)[i] = (*grid.v)(index2);
                        ++index2;
                    }
                }


            }
    };
}

#endif
