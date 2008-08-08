/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
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
   template <typename Tag_>
       struct UpdateVelocityDirectionsGrid<Tag_, lbm_boundary_types::NOSLIP>
       {
           /**
            * \name Velocity Update
            *
            * \brief Computes bla.
            *
            */
           template<typename DT1_>
               static void value(PackedGridData<D2Q9, DT1_> & data, PackedGridInfo<D2Q9> & info)
               {
                   CONTEXT("When updating velocity directions:");
                   for (unsigned long begin(0) ; begin != info.limits->size() - 1 ; ++begin)
                   {
                       for (unsigned long i((*info.limits)[begin]) ; i != (*info.limits)[begin + 1] ; ++i)
                       {
                           if(((*info.types)[begin] & 1<<0) == 1<<0)
                               (*data.f_temp_1)[i] = (*data.f_temp_5)[i];
                           if(((*info.types)[begin] & 1<<1) == 1<<0)
                               (*data.f_temp_2)[i] = (*data.f_temp_6)[i];
                           if(((*info.types)[begin] & 1<<2) == 1<<0)
                               (*data.f_temp_3)[i] = (*data.f_temp_7)[i];
                           if(((*info.types)[begin] & 1<<3) == 1<<0)
                               (*data.f_temp_4)[i] = (*data.f_temp_8)[i];
                           if(((*info.types)[begin] & 1<<4) == 1<<0)
                               (*data.f_temp_5)[i] = (*data.f_temp_1)[i];
                           if(((*info.types)[begin] & 1<<5) == 1<<0)
                               (*data.f_temp_6)[i] = (*data.f_temp_2)[i];
                           if(((*info.types)[begin] & 1<<6) == 1<<0)
                               (*data.f_temp_7)[i] = (*data.f_temp_3)[i];
                           if(((*info.types)[begin] & 1<<7) == 1<<0)
                               (*data.f_temp_8)[i] = (*data.f_temp_4)[i];
                       }
                   }
               }
       };
}
#endif
