/* vim: set number sw=4 sts=4 et nofoldenable : */

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


#ifndef LBM_GUARD_COLLIDE_STREAM_GRID_HH
#define LBM_GUARD_COLLIDE_STREAM_GRID_HH 1


/**
 * \file
 * Implementation of collision and streaming modules used by  LBM - (SWE, NavSto) solvers (using PackedGrid).
 *
 * \ingroup grpliblbm
 **/

#include <honei/lbm/tags.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/lbm/grid.hh>
#include <cmath>

using namespace honei::lbm;

namespace honei
{
    template <typename Tag_, typename Application_, typename BoundaryType_, typename LatticeType_>
    struct CollideStreamGrid
    {
    };

    /**
     * \brief Collision and streaming module for LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct CollideStreamGrid<Tag_, lbm_applications::LABSWE, lbm_boundary_types::NOSLIP, lbm_lattice_types::D2Q9>
    {
        /**
         * \name Collision and Streaming for direction 1..
         *
         * \brief Solves the LB equation.
         *
         * \param tau The relaxation time.
         */
        template <typename DT1_, typename DT2_>
        static void value(
                          PackedGridInfo<lbm_lattice_types::D2Q9> & info,
                          PackedGridData<lbm_lattice_types::D2Q9, DT1_> & data,
                          DT2_ tau)
        {
            CONTEXT("When performing collision and streaming:");

            for (unsigned long begin(0) ; begin != info.limits->size() - 1 ; ++begin)
            {
                for (unsigned long i((*info.limits)[begin]), offset(0) ; i != (*info.limits)[begin + 1] ; ++i, ++offset)
                {
                    (*data.f_temp_0)[i] = (*data.f_0)[i] - ((*data.f_0)[i] - (*data.f_eq_0)[i])/tau;
                }
            }

            for (unsigned long begin(0), half(0) ; begin < info.dir_index_1->size() - 1; begin+=2, ++half)
            {
                for (unsigned long i((*info.dir_index_1)[begin]), offset(0) ; i < (*info.dir_index_1)[begin + 1] ; ++i, ++offset)
                {
                    (*data.f_temp_1)[(*info.dir_1)[half] + offset] = (*data.f_1)[i] - ((*data.f_1)[i] - (*data.f_eq_1)[i])/tau;
                }
            }

            for (unsigned long begin(0), half(0) ; begin < info.dir_index_2->size() - 1; begin+=2, ++half)
            {
                for (unsigned long i((*info.dir_index_2)[begin]), offset(0) ; i < (*info.dir_index_2)[begin + 1] ; ++i, ++offset)
                {
                    (*data.f_temp_2)[(*info.dir_2)[half] + offset] = (*data.f_2)[i] - ((*data.f_2)[i] - (*data.f_eq_2)[i])/tau;
                }
            }

            for (unsigned long begin(0), half(0) ; begin < info.dir_index_3->size() - 1; begin+=2, ++half)
            {
                for (unsigned long i((*info.dir_index_3)[begin]), offset(0) ; i < (*info.dir_index_3)[begin + 1] ; ++i, ++offset)
                {
                    (*data.f_temp_3)[(*info.dir_3)[half] + offset] = (*data.f_3)[i] - ((*data.f_3)[i] - (*data.f_eq_3)[i])/tau;
                }
            }

            for (unsigned long begin(0), half(0) ; begin < info.dir_index_4->size() - 1; begin+=2, ++half)
            {
                for (unsigned long i((*info.dir_index_4)[begin]), offset(0) ; i < (*info.dir_index_4)[begin + 1] ; ++i, ++offset)
                {
                    (*data.f_temp_4)[(*info.dir_4)[half] + offset] = (*data.f_4)[i] - ((*data.f_4)[i] - (*data.f_eq_4)[i])/tau;
                }
            }

            for (unsigned long begin(0), half(0) ; begin < info.dir_index_5->size() - 1; begin+=2, ++half)
            {
                for (unsigned long i((*info.dir_index_5)[begin]), offset(0) ; i < (*info.dir_index_5)[begin + 1] ; ++i, ++offset)
                {
                    (*data.f_temp_5)[(*info.dir_5)[half] + offset] = (*data.f_5)[i] - ((*data.f_5)[i] - (*data.f_eq_5)[i])/tau;
                }
            }

            for (unsigned long begin(0), half(0) ; begin < info.dir_index_6->size() - 1; begin+=2, ++half)
            {
                for (unsigned long i((*info.dir_index_6)[begin]), offset(0) ; i < (*info.dir_index_6)[begin + 1] ; ++i, ++offset)
                {
                    (*data.f_temp_6)[(*info.dir_6)[half] + offset] = (*data.f_6)[i] - ((*data.f_6)[i] - (*data.f_eq_6)[i])/tau;
                }
            }

            for (unsigned long begin(0), half(0) ; begin < info.dir_index_7->size() - 1; begin+=2, ++half)
            {
                for (unsigned long i((*info.dir_index_7)[begin]), offset(0) ; i < (*info.dir_index_7)[begin + 1] ; ++i, ++offset)
                {
                    (*data.f_temp_7)[(*info.dir_7)[half] + offset] = (*data.f_7)[i] - ((*data.f_7)[i] - (*data.f_eq_7)[i])/tau;
                }
            }

            for (unsigned long begin(0), half(0) ; begin < info.dir_index_8->size() - 1; begin+=2, ++half)
            {
                for (unsigned long i((*info.dir_index_8)[begin]), offset(0) ; i < (*info.dir_index_8)[begin + 1] ; ++i, ++offset)
                {
                    (*data.f_temp_8)[(*info.dir_8)[half] + offset] = (*data.f_8)[i] - ((*data.f_8)[i] - (*data.f_eq_8)[i])/tau;
                }
            }
        }
    };
}
#endif
