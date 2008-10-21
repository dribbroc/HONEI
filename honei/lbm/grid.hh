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


#ifndef LBM_GUARD_GRID_HH
#define LBM_GUARD_GRID_HH 1

#include <honei/lbm/tags.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/dense_matrix.hh>

/**
 * \file
 * Definition of LBM grid and related classes.
 *
 * \ingroup grpliblbm
 **/


using namespace honei;
using namespace lbm;
using namespace lbm_lattice_types;

namespace honei
{
    template <typename LatticeType_, typename DT_> class Grid
    {
    };

    template <typename LatticeType_> class PackedGridInfo
    {
    };

    template <typename LatticeType_, typename DT_> class PackedGridData
    {
    };

    template <typename LatticeType_> class PackedGridFringe
    {
    };

    template <typename DT_> class Grid<D2Q9, DT_>
    {
        public:
            Grid():
                obstacles(0),
                h(0),
                b(0),
                b_x(0),
                b_y(0),
                u(0),
                v(0)
            {
            }
            void destroy()
            {
                delete obstacles;
                delete h;
                delete b;
                delete b_x;
                delete b_y;
                delete u;
                delete v;

                obstacles = 0;
                h = 0;
                b = 0;
                b_x = 0;
                b_y = 0;
                u = 0;
                v = 0;
            }
            DenseMatrix<bool> * obstacles;
            DenseMatrix<DT_> * h;
            DenseMatrix<DT_> * b;
            DenseMatrix<DT_> * b_x;
            DenseMatrix<DT_> * b_y;
            DenseMatrix<DT_> * u;
            DenseMatrix<DT_> * v;
    };

    template <> class PackedGridInfo<D2Q9>
    {
        public:
            PackedGridInfo():
                offset(0),
                limits(0),
                types(0),
                dir_1(0),
                dir_2(0),
                dir_3(0),
                dir_4(0),
                dir_5(0),
                dir_6(0),
                dir_7(0),
                dir_8(0),
                dir_index_1(0),
                dir_index_2(0),
                dir_index_3(0),
                dir_index_4(0),
                dir_index_5(0),
                dir_index_6(0),
                dir_index_7(0),
                dir_index_8(0)
            {
            }

            void destroy()
            {
                delete limits;
                delete types;
                delete dir_1;
                delete dir_2;
                delete dir_3;
                delete dir_4;
                delete dir_5;
                delete dir_6;
                delete dir_7;
                delete dir_8;
                delete dir_index_1;
                delete dir_index_2;
                delete dir_index_3;
                delete dir_index_4;
                delete dir_index_5;
                delete dir_index_6;
                delete dir_index_7;
                delete dir_index_8;

                limits = 0;
                types = 0;
                dir_1 = 0;
                dir_2 = 0;
                dir_3 = 0;
                dir_4 = 0;
                dir_5 = 0;
                dir_6 = 0;
                dir_7 = 0;
                dir_8 = 0;
                dir_index_1 = 0;
                dir_index_2 = 0;
                dir_index_3 = 0;
                dir_index_4 = 0;
                dir_index_5 = 0;
                dir_index_6 = 0;
                dir_index_7 = 0;
                dir_index_8 = 0;
            }

            unsigned long offset;
            DenseVector<unsigned long> * limits;
            /// \todo Use std::bitset instead of unsigned long for types?
            DenseVector<unsigned long> * types;
            DenseVector<unsigned long> * dir_1;
            DenseVector<unsigned long> * dir_2;
            DenseVector<unsigned long> * dir_3;
            DenseVector<unsigned long> * dir_4;
            DenseVector<unsigned long> * dir_5;
            DenseVector<unsigned long> * dir_6;
            DenseVector<unsigned long> * dir_7;
            DenseVector<unsigned long> * dir_8;
            DenseVector<unsigned long> * dir_index_1;
            DenseVector<unsigned long> * dir_index_2;
            DenseVector<unsigned long> * dir_index_3;
            DenseVector<unsigned long> * dir_index_4;
            DenseVector<unsigned long> * dir_index_5;
            DenseVector<unsigned long> * dir_index_6;
            DenseVector<unsigned long> * dir_index_7;
            DenseVector<unsigned long> * dir_index_8;
    };

    template <typename DT_> class PackedGridData<D2Q9, DT_>
    {
        public:

            PackedGridData():
                h(0),
                b_x(0),
                b_y(0),
                u(0),
                v(0),
                f_0(0),
                f_1(0),
                f_2(0),
                f_3(0),
                f_4(0),
                f_5(0),
                f_6(0),
                f_7(0),
                f_8(0),
                f_eq_0(0),
                f_eq_1(0),
                f_eq_2(0),
                f_eq_3(0),
                f_eq_4(0),
                f_eq_5(0),
                f_eq_6(0),
                f_eq_7(0),
                f_eq_8(0),
                f_temp_0(0),
                f_temp_1(0),
                f_temp_2(0),
                f_temp_3(0),
                f_temp_4(0),
                f_temp_5(0),
                f_temp_6(0),
                f_temp_7(0),
                f_temp_8(0),

                distribution_x(0),
                distribution_y(0)
            {
            }

            void destroy()
            {
                delete h;
                delete b_x;
                delete b_y;
                delete u;
                delete v;
                delete f_0;
                delete f_1;
                delete f_2;
                delete f_3;
                delete f_4;
                delete f_5;
                delete f_6;
                delete f_7;
                delete f_8;
                delete f_eq_0;
                delete f_eq_1;
                delete f_eq_2;
                delete f_eq_3;
                delete f_eq_4;
                delete f_eq_5;
                delete f_eq_6;
                delete f_eq_7;
                delete f_eq_8;
                delete f_temp_0;
                delete f_temp_1;
                delete f_temp_2;
                delete f_temp_3;
                delete f_temp_4;
                delete f_temp_5;
                delete f_temp_6;
                delete f_temp_7;
                delete f_temp_8;

                delete distribution_x;
                delete distribution_y;

                h = 0;
                b_x = 0;
                b_y = 0;
                u = 0;
                v = 0;
                f_0 = 0;
                f_1 = 0;
                f_2 = 0;
                f_3 = 0;
                f_4 = 0;
                f_5 = 0;
                f_6 = 0;
                f_7 = 0;
                f_8 = 0;
                f_eq_0 = 0;
                f_eq_1 = 0;
                f_eq_2 = 0;
                f_eq_3 = 0;
                f_eq_4 = 0;
                f_eq_5 = 0;
                f_eq_6 = 0;
                f_eq_7 = 0;
                f_eq_8 = 0;
                f_temp_0 = 0;
                f_temp_1 = 0;
                f_temp_2 = 0;
                f_temp_3 = 0;
                f_temp_4 = 0;
                f_temp_5 = 0;
                f_temp_6 = 0;
                f_temp_7 = 0;
                f_temp_8 = 0;

                distribution_x = 0;
                distribution_y = 0;
            }
            DenseVector<DT_> * h;
            DenseVector<DT_> * b_x;
            DenseVector<DT_> * b_y;
            DenseVector<DT_> * u;
            DenseVector<DT_> * v;

            DenseVector<DT_> * f_0;
            DenseVector<DT_> * f_1;
            DenseVector<DT_> * f_2;
            DenseVector<DT_> * f_3;
            DenseVector<DT_> * f_4;
            DenseVector<DT_> * f_5;
            DenseVector<DT_> * f_6;
            DenseVector<DT_> * f_7;
            DenseVector<DT_> * f_8;

            DenseVector<DT_> * f_eq_0;
            DenseVector<DT_> * f_eq_1;
            DenseVector<DT_> * f_eq_2;
            DenseVector<DT_> * f_eq_3;
            DenseVector<DT_> * f_eq_4;
            DenseVector<DT_> * f_eq_5;
            DenseVector<DT_> * f_eq_6;
            DenseVector<DT_> * f_eq_7;
            DenseVector<DT_> * f_eq_8;

            DenseVector<DT_> * f_temp_0;
            DenseVector<DT_> * f_temp_1;
            DenseVector<DT_> * f_temp_2;
            DenseVector<DT_> * f_temp_3;
            DenseVector<DT_> * f_temp_4;
            DenseVector<DT_> * f_temp_5;
            DenseVector<DT_> * f_temp_6;
            DenseVector<DT_> * f_temp_7;
            DenseVector<DT_> * f_temp_8;

            DenseVector<DT_> * distribution_x;
            DenseVector<DT_> * distribution_y;
    };

    template <> class PackedGridFringe<D2Q9>
    {
        public:
            PackedGridFringe():
                h_index(0),
                h_targets(0),
                dir_index_1(0),
                dir_targets_1(0),
                dir_index_2(0),
                dir_targets_2(0),
                dir_index_3(0),
                dir_targets_3(0),
                dir_index_4(0),
                dir_targets_4(0),
                dir_index_5(0),
                dir_targets_5(0),
                dir_index_6(0),
                dir_targets_6(0),
                dir_index_7(0),
                dir_targets_7(0),
                dir_index_8(0),
                dir_targets_8(0)
            {
            }

            void destroy()
            {
                delete h_index;
                delete h_targets;
                delete dir_index_1;
                delete dir_targets_1;
                delete dir_index_2;
                delete dir_targets_2;
                delete dir_index_3;
                delete dir_targets_3;
                delete dir_index_4;
                delete dir_targets_4;
                delete dir_index_5;
                delete dir_targets_5;
                delete dir_index_6;
                delete dir_targets_6;
                delete dir_index_7;
                delete dir_targets_7;
                delete dir_index_8;
                delete dir_targets_8;

                h_index = 0;
                h_targets = 0;
                dir_index_1 = 0;
                dir_targets_1 = 0;
                dir_index_2 = 0;
                dir_targets_2 = 0;
                dir_index_3 = 0;
                dir_targets_3 = 0;
                dir_index_4 = 0;
                dir_targets_4 = 0;
                dir_index_5 = 0;
                dir_targets_5 = 0;
                dir_index_6 = 0;
                dir_targets_6 = 0;
                dir_index_7 = 0;
                dir_targets_7 = 0;
                dir_index_8 = 0;
                dir_targets_8 = 0;
            }

            DenseVector<unsigned long> * h_index;
            DenseVector<unsigned long> * h_targets;
            DenseVector<unsigned long> * dir_index_1;
            DenseVector<unsigned long> * dir_targets_1;
            DenseVector<unsigned long> * dir_index_2;
            DenseVector<unsigned long> * dir_targets_2;
            DenseVector<unsigned long> * dir_index_3;
            DenseVector<unsigned long> * dir_targets_3;
            DenseVector<unsigned long> * dir_index_4;
            DenseVector<unsigned long> * dir_targets_4;
            DenseVector<unsigned long> * dir_index_5;
            DenseVector<unsigned long> * dir_targets_5;
            DenseVector<unsigned long> * dir_index_6;
            DenseVector<unsigned long> * dir_targets_6;
            DenseVector<unsigned long> * dir_index_7;
            DenseVector<unsigned long> * dir_targets_7;
            DenseVector<unsigned long> * dir_index_8;
            DenseVector<unsigned long> * dir_targets_8;
    };
}

#endif
