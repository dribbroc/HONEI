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


#ifndef LBM_GUARD_UPDATE_VELOCITY_DIRECTIONS_GRID_HH
#define LBM_GUARD_UPDATE_VELOCITY_DIRECTIONS_GRID_HH 1

/**
 * \file
 * Implementation of update velocity directions used by  LBM - (SWE) grid based solvers.
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
   template <typename Tag_,
              typename BoundaryType_>
    struct UpdateVelocityDirectionsGrid
    {
    };

    /**
     * \brief Temp update of velocity directions for use with LABSWE.
     *
     * \ingroup grplbmoperations
     */
   template <>
       struct UpdateVelocityDirectionsGrid<tags::CPU, lbm_boundary_types::NOSLIP>
       {
           /**
            * \name Velocity Update
            *
            * \brief Computes boundary velocity values.
            *
            */
           template<typename DT1_>
               static void value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, DT1_> & data)
               {
                   CONTEXT("When updating velocity directions:");

                   info.limits->lock(lm_read_only);
                   info.types->lock(lm_read_only);

                   data.f_temp_1->lock(lm_read_and_write);
                   data.f_temp_2->lock(lm_read_and_write);
                   data.f_temp_3->lock(lm_read_and_write);
                   data.f_temp_4->lock(lm_read_and_write);
                   data.f_temp_5->lock(lm_read_and_write);
                   data.f_temp_6->lock(lm_read_and_write);
                   data.f_temp_7->lock(lm_read_and_write);
                   data.f_temp_8->lock(lm_read_and_write);

                   for (unsigned long begin(0) ; begin < info.limits->size() - 1 ; ++begin)
                   {
                       if(((*info.types)[begin] & 1<<0) == 1<<0)
                           for (unsigned long i((*info.limits)[begin]) ; i != (*info.limits)[begin + 1] ; ++i)
                           {
                               (*data.f_temp_5)[i] = (*data.f_temp_1)[i];
                           }
                       if(((*info.types)[begin] & 1<<1) == 1<<1)
                           for (unsigned long i((*info.limits)[begin]) ; i != (*info.limits)[begin + 1] ; ++i)
                           {
                               (*data.f_temp_6)[i] = (*data.f_temp_2)[i];
                           }
                       if(((*info.types)[begin] & 1<<2) == 1<<2)
                           for (unsigned long i((*info.limits)[begin]) ; i != (*info.limits)[begin + 1] ; ++i)
                           {
                               (*data.f_temp_7)[i] = (*data.f_temp_3)[i];
                           }
                       if(((*info.types)[begin] & 1<<3) == 1<<3)
                           for (unsigned long i((*info.limits)[begin]) ; i != (*info.limits)[begin + 1] ; ++i)
                           {
                               (*data.f_temp_8)[i] = (*data.f_temp_4)[i];
                           }
                       if(((*info.types)[begin] & 1<<4) == 1<<4)
                           for (unsigned long i((*info.limits)[begin]) ; i != (*info.limits)[begin + 1] ; ++i)
                           {
                               (*data.f_temp_1)[i] = (*data.f_temp_5)[i];
                           }
                       if(((*info.types)[begin] & 1<<5) == 1<<5)
                           for (unsigned long i((*info.limits)[begin]) ; i != (*info.limits)[begin + 1] ; ++i)
                           {
                               (*data.f_temp_2)[i] = (*data.f_temp_6)[i];
                           }
                       if(((*info.types)[begin] & 1<<6) == 1<<6)
                           for (unsigned long i((*info.limits)[begin]) ; i != (*info.limits)[begin + 1] ; ++i)
                           {
                               (*data.f_temp_3)[i] = (*data.f_temp_7)[i];
                           }
                       if(((*info.types)[begin] & 1<<7) == 1<<7)
                           for (unsigned long i((*info.limits)[begin]) ; i != (*info.limits)[begin + 1] ; ++i)
                           {
                               (*data.f_temp_4)[i] = (*data.f_temp_8)[i];
                           }

                       // Corners
                       if(((*info.types)[begin] & 1<<2) == 1<<2 && ((*info.types)[begin] & 1<<4) == 1<<4)
                           for (unsigned long i((*info.limits)[begin]) ; i != (*info.limits)[begin + 1] ; ++i)
                           {
                               (*data.f_temp_2)[i] = (*data.f_temp_8)[i];
                               (*data.f_temp_6)[i] = (*data.f_temp_8)[i];
                           }
                       if(((*info.types)[begin] & 1<<4) == 1<<4 && ((*info.types)[begin] & 1<<6) == 1<<6)
                           for (unsigned long i((*info.limits)[begin]) ; i != (*info.limits)[begin + 1] ; ++i)
                           {
                               (*data.f_temp_4)[i] = (*data.f_temp_2)[i];
                               (*data.f_temp_8)[i] = (*data.f_temp_2)[i];
                           }
                       if(((*info.types)[begin] & 1<<0) == 1<<0 && ((*info.types)[begin] & 1<<6) == 1<<6)
                           for (unsigned long i((*info.limits)[begin]) ; i != (*info.limits)[begin + 1] ; ++i)
                           {
                               (*data.f_temp_2)[i] = (*data.f_temp_4)[i];
                               (*data.f_temp_6)[i] = (*data.f_temp_4)[i];
                           }
                       if(((*info.types)[begin] & 1<<0) == 1<<0 && ((*info.types)[begin] & 1<<2) == 1<<2)
                           for (unsigned long i((*info.limits)[begin]) ; i != (*info.limits)[begin + 1] ; ++i)
                           {
                               (*data.f_temp_4)[i] = (*data.f_temp_6)[i];
                               (*data.f_temp_8)[i] = (*data.f_temp_6)[i];
                           }
                   }

                   info.limits->unlock(lm_read_only);
                   info.types->unlock(lm_read_only);

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
       struct UpdateVelocityDirectionsGrid<tags::GPU::CUDA, lbm_boundary_types::NOSLIP>
       {
           static void value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data);
       };

    template <>
       struct UpdateVelocityDirectionsGrid<tags::CPU::SSE, lbm_boundary_types::NOSLIP>
       {
           static void value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, float> & data);
           static void value(PackedGridInfo<D2Q9> & info, PackedGridData<D2Q9, double> & data);
       };
}
#endif
