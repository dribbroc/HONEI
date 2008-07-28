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
#include <set>
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
    template <typename LatticeType_, typename DT_> struct GridPacker
    {
    };

    template <typename DT_> struct GridPacker<D2Q9, DT_>
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

                if(i == (*grid.obstacles).rows() -1)
                    south = true;

                if(j == 0)
                    west = true;

                if(j == (*grid.obstacles).columns() -1)
                    east = true;

                if( i > 0 && (*grid.obstacles)(i - 1, j))
                    north = true;

                if( i <  (*grid.obstacles).rows() - 1 && (*grid.obstacles)(i + 1, j))
                    south = true;

                if( j > 0 && (*grid.obstacles)(i, j - 1))
                    west = true;

                if( j <  (*grid.obstacles).columns() - 1 && (*grid.obstacles)(i, j + 1))
                    east = true;

                if(!north && !south && !west && !east)
                    return bt_none;

                if(north && south && west && east)
                    throw InternalError("Found inclosed fluid!");

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
                            if( i - 1 >= 0
                                    && i + 1 < grid.obstacles->rows()
                                    && j - 1 >= 0
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

                std::set<unsigned long> temp_limits;

                temp_limits.insert(0);
                for(unsigned long i(0); i < grid.obstacles->rows(); ++i)
                {
                    for(unsigned long j(0); j < grid.obstacles->columns(); ++j)
                    {
                        if(j > 0)
                        {
                            if (_element_type(i, j, grid) != _element_type(i, j - 1, grid))
                            {
                                temp_limits.insert(i * grid.obstacles->columns() + j);
                            }
                        }
                        else
                        {
                            temp_limits.insert(i * grid.obstacles->columns() + j);
                        }
                    }
                }
                temp_limits.insert(grid.obstacles->size());

                info.limits = new DenseVector<unsigned long>(temp_limits.size());
                info.types = new DenseVector<unsigned long>(temp_limits.size());

                unsigned long index(0);
                for (std::set<unsigned long>::iterator i(temp_limits.begin()) ; i != temp_limits.end() ; ++i)
                {
                    (*info.limits)[index] = *i;
                    (*info.types)[index] = _element_type(*i / grid.obstacles->columns(), *i % grid.obstacles->columns(), grid);
                    ++index;
                }
                for (unsigned long i(0) ; i < info.limits->size() ; ++i)
                {
                    std::cout<<(*info.limits)[i] << " "<<(*info.types)[i]<<std::endl;
                }

            }
    };
}

#endif
