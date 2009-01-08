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

                info.dir_index_1->lock(lm_read_only);
                info.dir_index_2->lock(lm_read_only);
                info.dir_index_3->lock(lm_read_only);
                info.dir_index_4->lock(lm_read_only);
                info.dir_index_5->lock(lm_read_only);
                info.dir_index_6->lock(lm_read_only);
                info.dir_index_7->lock(lm_read_only);
                info.dir_index_8->lock(lm_read_only);

                info.dir_1->lock(lm_read_only);
                info.dir_2->lock(lm_read_only);
                info.dir_3->lock(lm_read_only);
                info.dir_4->lock(lm_read_only);
                info.dir_5->lock(lm_read_only);
                info.dir_6->lock(lm_read_only);
                info.dir_7->lock(lm_read_only);
                info.dir_8->lock(lm_read_only);

                data.h->lock(lm_read_only);
                data.b->lock(lm_read_only);
                data.distribution_x->lock(lm_read_only);
                data.distribution_y->lock(lm_read_only);

                data.f_temp_1->lock(lm_read_and_write);
                data.f_temp_2->lock(lm_read_and_write);
                data.f_temp_3->lock(lm_read_and_write);
                data.f_temp_4->lock(lm_read_and_write);
                data.f_temp_5->lock(lm_read_and_write);
                data.f_temp_6->lock(lm_read_and_write);
                data.f_temp_7->lock(lm_read_and_write);
                data.f_temp_8->lock(lm_read_and_write);


                for (unsigned long begin(0), half(0) ; begin < info.dir_index_1->size() - 1; begin+=2, ++half)
                {
                    for (unsigned long i((*info.dir_index_1)[begin]), offset(0) ; i < (*info.dir_index_1)[begin + 1] ; ++i, ++offset)
                    {
                        /*(*data.f_temp_1)[i] += d_t / (6 * d_x / d_t) * ((*data.distribution_x)[1]) *
                             (-g * (((*data.h)[(*info.dir_1)[half] + offset]) - ((*data.b_x)[i]))/ DT1_(2.) *
                             ((((*data.b_x)[(*info.dir_1)[half] + offset]) - ((*data.b_x)[i]))/ DT1_(2.)));

                        (*data.f_temp_1)[i] += d_t / (6 * d_y / d_t) * ((*data.distribution_y)[1]) *
                            (- g * (((*data.h)[(*info.dir_1)[half] + offset]) - ((*data.b_y)[i]))/ DT1_(2.) *
                             ((((*data.b_y)[(*info.dir_1)[half] + offset]) - ((*data.b_y)[i]))/ DT1_(2.)));*/
                    }
                }

                info.dir_index_1->unlock(lm_read_only);
                info.dir_index_2->unlock(lm_read_only);
                info.dir_index_3->unlock(lm_read_only);
                info.dir_index_4->unlock(lm_read_only);
                info.dir_index_5->unlock(lm_read_only);
                info.dir_index_6->unlock(lm_read_only);
                info.dir_index_7->unlock(lm_read_only);
                info.dir_index_8->unlock(lm_read_only);

                info.dir_1->unlock(lm_read_only);
                info.dir_2->unlock(lm_read_only);
                info.dir_3->unlock(lm_read_only);
                info.dir_4->unlock(lm_read_only);
                info.dir_5->unlock(lm_read_only);
                info.dir_6->unlock(lm_read_only);
                info.dir_7->unlock(lm_read_only);
                info.dir_8->unlock(lm_read_only);

                data.h->unlock(lm_read_only);
                data.b->unlock(lm_read_only);
                data.distribution_x->unlock(lm_read_only);
                data.distribution_y->unlock(lm_read_only);

                data.f_temp_1->unlock(lm_read_and_write);
                data.f_temp_2->unlock(lm_read_and_write);
                data.f_temp_3->unlock(lm_read_and_write);
                data.f_temp_4->unlock(lm_read_and_write);
                data.f_temp_5->unlock(lm_read_and_write);
                data.f_temp_6->unlock(lm_read_and_write);
                data.f_temp_7->unlock(lm_read_and_write);
                data.f_temp_8->unlock(lm_read_and_write);
            }
    };

    template <>
    struct ForceGrid<tags::GPU::CUDA, lbm_applications::LABSWE, lbm_source_types::CENTRED, lbm_source_schemes::CENTRALDIFF>
    {
        static void value(PackedGridData<D2Q9, float> & data, PackedGridInfo<D2Q9> & info, float g, float d_x, float d_y, float d_t);
    };

    template <>
    struct ForceGrid<tags::CPU::SSE, lbm_applications::LABSWE, lbm_source_types::CENTRED, lbm_source_schemes::CENTRALDIFF>
    {
        static void value(PackedGridData<D2Q9, float> & data, PackedGridInfo<D2Q9> & info, float g, float d_x, float d_y, float d_t)
        {
            ForceGrid<tags::CPU, lbm_applications::LABSWE, lbm_source_types::CENTRED, lbm_source_schemes::CENTRALDIFF>::value
                (data, info, g, d_x, d_y, d_t);
        }
        static void value(PackedGridData<D2Q9, double> & data, PackedGridInfo<D2Q9> & info, double g, double d_x, double d_y, double d_t)
        {
            ForceGrid<tags::CPU, lbm_applications::LABSWE, lbm_source_types::CENTRED, lbm_source_schemes::CENTRALDIFF>::value
                (data, info, g, d_x, d_y, d_t);
        }
    };
}
#endif
