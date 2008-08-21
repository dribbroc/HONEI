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
                data.h = new DenseVector<DT_>(fluid_count);
                data.u = new DenseVector<DT_>(fluid_count);
                data.v = new DenseVector<DT_>(fluid_count);

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

                info.limits = new DenseVector<unsigned long>(temp_limits.size());
                info.types = new DenseVector<unsigned long>(temp_limits.size());
                info.dir_1 = new DenseVector<unsigned long>(dir_1.size());
                info.dir_2 = new DenseVector<unsigned long>(dir_1.size());
                info.dir_3 = new DenseVector<unsigned long>(dir_1.size());
                info.dir_4 = new DenseVector<unsigned long>(dir_1.size());
                info.dir_5 = new DenseVector<unsigned long>(dir_1.size());
                info.dir_6 = new DenseVector<unsigned long>(dir_1.size());
                info.dir_7 = new DenseVector<unsigned long>(dir_1.size());
                info.dir_8 = new DenseVector<unsigned long>(dir_1.size());

                unsigned long index2(0);
                for (std::vector<unsigned long>::iterator i(temp_limits.begin()), j(temp_types.begin()) ; i != temp_limits.end() ; ++i, ++j)
                {
                    (*info.limits)[index2] = *i;
                    (*info.types)[index2] = *j;
                    ++index2;
                }
                for (unsigned long i(0) ; i < dir_1.size() ; ++i)
                {
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
                        while ((*grid.obstacles)(index2 / grid.obstacles->columns(), index2 % grid.obstacles->columns()))
                        {
                            ++index2;
                        }
                        (*data.h)[i] = (*grid.h)(index2 / grid.obstacles->columns(), index2 % grid.obstacles->columns());
                        (*data.u)[i] = (*grid.u)(index2 / grid.obstacles->columns(), index2 % grid.obstacles->columns());
                        (*data.v)[i] = (*grid.v)(index2 / grid.obstacles->columns(), index2 % grid.obstacles->columns());
                        ++index2;
                    }
                }
            }

            static void unpack(Grid<D2Q9, DT_> & grid, PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT_> & data)
            {
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
}

#endif
