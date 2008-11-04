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


#ifndef LBM_GUARD_FORCE_GRID_HH
#define LBM_GUARD_FORCE_GRID_HH 1

/**
 * \file
 * Implementation of source modules used by  LBM - (SWE) solvers using a
 * PackedGrid.
 *
 * \ingroup grpliblbm
 **/

#include <honei/lbm/tags.hh>
#include <honei/lbm/grid.hh>
#include <honei/la/dense_vector.hh>
#include <cmath>
using namespace honei::lbm;

namespace honei
{
   template <typename Tag_,
              typename App_,
              typename SourceType_,
              typename SourceScheme_>
    struct ForceGrid
    {
    };

    /**
     * \brief Simple force term for use with LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <>
    struct ForceGrid<tags::CPU, lbm_applications::LABSWE, lbm_source_types::CENTRED, lbm_source_schemes::CENTRALDIFF>
    {
        /**
         * \name Force term.
         *
         * \brief Computes a simple source term value.
         *
         * \param data Our packed grid data object.
         * \param info Our info object.
         * \param d_x Our delta x.
         * \param d_y Our delta y.
         * \param d_t Our delta t.
         * \param g The gravitational constant to be used.
         *
         */
        template<typename DT1_, typename DT2_>
            static void value(PackedGridData<D2Q9, DT1_> & data, PackedGridInfo<D2Q9> & info, DT2_ g, DT2_ d_x, DT2_ d_y, DT2_ d_t)
            {
                CONTEXT("When computing LABSWE force term:");
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_1->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_1)[begin]), offset(0) ; i < (*info.dir_index_1)[begin + 1] ; ++i, ++offset)
                    {
                        (*data.f_temp_1)[i] += d_t / (6 * d_x / d_t) * ((*data.distribution_x)[1]) *
                             (-g * (((*data.h)[(*info.dir_1)[half] + offset]) - ((*data.b_x)[i]))/ DT1_(2.) *
                             ((((*data.b_x)[(*info.dir_1)[half] + offset]) - ((*data.b_x)[i]))/ DT1_(2.)));

                        (*data.f_temp_1)[i] += d_t / (6 * d_y / d_t) * ((*data.distribution_y)[1]) *
                            (- g * (((*data.h)[(*info.dir_1)[half] + offset]) - ((*data.b_y)[i]))/ DT1_(2.) *
                             ((((*data.b_y)[(*info.dir_1)[half] + offset]) - ((*data.b_y)[i]))/ DT1_(2.)));
                    }
                }
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_2->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_2)[begin]), offset(0) ; i < (*info.dir_index_2)[begin + 1] ; ++i, ++offset)
                    {
                        (*data.f_temp_2)[i] += d_t / (6 * d_x / d_t) * ((*data.distribution_x)[2]) *
                             (-g * (((*data.h)[(*info.dir_2)[half] + offset]) - ((*data.b_x)[i]))/ DT1_(2.) *
                             ((((*data.b_x)[(*info.dir_2)[half] + offset]) - ((*data.b_x)[i]))/ DT1_(2.)));

                        (*data.f_temp_2)[i] += d_t / (6 * d_y / d_t) * ((*data.distribution_y)[2]) *
                            (- g * (((*data.h)[(*info.dir_2)[half] + offset]) - ((*data.b_y)[i]))/ DT1_(2.) *
                             ((((*data.b_y)[(*info.dir_2)[half] + offset]) - ((*data.b_y)[i]))/ DT1_(2.)));
                    }
                }
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_3->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_3)[begin]), offset(0) ; i < (*info.dir_index_3)[begin + 1] ; ++i, ++offset)
                    {
                        (*data.f_temp_3)[i] += d_t / (6 * d_x / d_t) * ((*data.distribution_x)[3]) *
                             (-g * (((*data.h)[(*info.dir_3)[half] + offset]) - ((*data.b_x)[i]))/ DT1_(2.) *
                             ((((*data.b_x)[(*info.dir_3)[half] + offset]) - ((*data.b_x)[i]))/ DT1_(2.)));

                        (*data.f_temp_3)[i] += d_t / (6 * d_y / d_t) * ((*data.distribution_y)[3]) *
                            (- g * (((*data.h)[(*info.dir_3)[half] + offset]) - ((*data.b_y)[i]))/ DT1_(2.) *
                             ((((*data.b_y)[(*info.dir_3)[half] + offset]) - ((*data.b_y)[i]))/ DT1_(2.)));
                    }
                }
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_4->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_4)[begin]), offset(0) ; i < (*info.dir_index_4)[begin + 1] ; ++i, ++offset)
                    {
                        (*data.f_temp_4)[i] += d_t / (6 * d_x / d_t) * ((*data.distribution_x)[4]) *
                             (-g * (((*data.h)[(*info.dir_4)[half] + offset]) - ((*data.b_x)[i]))/ DT1_(2.) *
                             ((((*data.b_x)[(*info.dir_4)[half] + offset]) - ((*data.b_x)[i]))/ DT1_(2.)));

                        (*data.f_temp_4)[i] += d_t / (6 * d_y / d_t) * ((*data.distribution_y)[4]) *
                            (- g * (((*data.h)[(*info.dir_4)[half] + offset]) - ((*data.b_y)[i]))/ DT1_(2.) *
                             ((((*data.b_y)[(*info.dir_4)[half] + offset]) - ((*data.b_y)[i]))/ DT1_(2.)));
                    }
                }
                for (unsigned long begin(0), half(0) ; begin < info.dir_index_5->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_5)[begin]), offset(0) ; i < (*info.dir_index_5)[begin + 1] ; ++i, ++offset)
                    {
                        (*data.f_temp_5)[i] += d_t / (6 * d_x / d_t) * ((*data.distribution_x)[5]) *
                             (-g * (((*data.h)[(*info.dir_5)[half] + offset]) - ((*data.b_x)[i]))/ DT1_(2.) *
                             ((((*data.b_x)[(*info.dir_5)[half] + offset]) - ((*data.b_x)[i]))/ DT1_(2.)));

                        (*data.f_temp_5)[i] += d_t / (6 * d_y / d_t) * ((*data.distribution_y)[5]) *
                            (- g * (((*data.h)[(*info.dir_5)[half] + offset]) - ((*data.b_y)[i]))/ DT1_(2.) *
                             ((((*data.b_y)[(*info.dir_5)[half] + offset]) - ((*data.b_y)[i]))/ DT1_(2.)));
                    }
                }

                for (unsigned long begin(0), half(0) ; begin < info.dir_index_6->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_6)[begin]), offset(0) ; i < (*info.dir_index_6)[begin + 1] ; ++i, ++offset)
                    {
                        (*data.f_temp_6)[i] += d_t / (6 * d_x / d_t) * ((*data.distribution_x)[6]) *
                             (-g * (((*data.h)[(*info.dir_6)[half] + offset]) - ((*data.b_x)[i]))/ DT1_(2.) *
                             ((((*data.b_x)[(*info.dir_6)[half] + offset]) - ((*data.b_x)[i]))/ DT1_(2.)));

                        (*data.f_temp_6)[i] += d_t / (6 * d_y / d_t) * ((*data.distribution_y)[6]) *
                            (- g * (((*data.h)[(*info.dir_6)[half] + offset]) - ((*data.b_y)[i]))/ DT1_(2.) *
                             ((((*data.b_y)[(*info.dir_6)[half] + offset]) - ((*data.b_y)[i]))/ DT1_(2.)));
                    }
                }

                for (unsigned long begin(0), half(0) ; begin < info.dir_index_7->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_7)[begin]), offset(0) ; i < (*info.dir_index_7)[begin + 1] ; ++i, ++offset)
                    {
                        (*data.f_temp_7)[i] += d_t / (6 * d_x / d_t) * ((*data.distribution_x)[7]) *
                             (-g * (((*data.h)[(*info.dir_7)[half] + offset]) - ((*data.b_x)[i]))/ DT1_(2.) *
                             ((((*data.b_x)[(*info.dir_7)[half] + offset]) - ((*data.b_x)[i]))/ DT1_(2.)));

                        (*data.f_temp_7)[i] += d_t / (6 * d_y / d_t) * ((*data.distribution_y)[7]) *
                            (- g * (((*data.h)[(*info.dir_7)[half] + offset]) - ((*data.b_y)[i]))/ DT1_(2.) *
                             ((((*data.b_y)[(*info.dir_7)[half] + offset]) - ((*data.b_y)[i]))/ DT1_(2.)));
                    }
                }

                for (unsigned long begin(0), half(0) ; begin < info.dir_index_8->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_8)[begin]), offset(0) ; i < (*info.dir_index_8)[begin + 1] ; ++i, ++offset)
                    {
                        (*data.f_temp_8)[i] += d_t / (6 * d_x / d_t) * ((*data.distribution_x)[8]) *
                             (-g * (((*data.h)[(*info.dir_8)[half] + offset]) - ((*data.b_x)[i]))/ DT1_(2.) *
                             ((((*data.b_x)[(*info.dir_8)[half] + offset]) - ((*data.b_x)[i]))/ DT1_(2.)));

                        (*data.f_temp_8)[i] += d_t / (6 * d_y / d_t) * ((*data.distribution_y)[8]) *
                            (- g * (((*data.h)[(*info.dir_8)[half] + offset]) - ((*data.b_y)[i]))/ DT1_(2.) *
                             ((((*data.b_y)[(*info.dir_8)[half] + offset]) - ((*data.b_y)[i]))/ DT1_(2.)));
                    }
                }
            }
    };
}
#endif
