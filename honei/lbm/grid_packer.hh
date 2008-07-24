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
        public:
            static void pack(Grid<D2Q9, DT_> & grid, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data)
            {
                unsigned long fluid_count(0);
                /*for(unsigned long i(0); i < grid.obstacles->rows(); ++i)
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
                }*/
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

                /*for(unsigned long i(0); i < grid.obstacles->rows(); ++i)
                {
                    for(unsigned long j(0); j < grid.obstacles->columns(); ++j)
                    {
                        if((*grid.obstacles)(i,j))
                            continue;
                        else
                        {

                        }
                    }
                }*/


            }
    };
}

#endif
