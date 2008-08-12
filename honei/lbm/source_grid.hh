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


#ifndef LBM_GUARD_SOURCE_GRID_HH
#define LBM_GUARD_SOURCE_GRID_HH 1

/**
 * \file
 * Implementation of source modules used by  LBM - (SWE) solvers using a
 * PackedGrid.
 *
 * \ingroup grpliblbm
 **/

#include <honei/lbm/tags.hh>
#include <grid.hh>
#include <honei/la/dense_vector.hh>
#include <cmath>
using namespace honei::lbm;

namespace honei
{
   template <typename Tag_,
              typename App_,
              typename SourceType_,
              typename SourceScheme_>
    struct SourceGrid
    {
    };

    /**
     * \brief Simple source term for use with LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct SourceGrid<Tag_, lbm_applications::LABSWE, lbm_source_types::SIMPLE, lbm_source_schemes::BASIC>
    {
        /**
         * \name Source term.
         *
         * \brief Computes a simple source term value.
         *
         * \param result The destination matrix.
         * \param h The height matrix.
         * \param dbx The matrix containing the slope values.
         * \param g The gravitational constant to be used.
         *
         */
        template<typename DT1_, typename DT2_, typename DT3_>
            static void value(PackedGridData<D2Q9, DT1_> & data, DenseVector<DT1_>& result, DenseVector<DT2_>& dbx, DT3_ g)
            {
                CONTEXT("When computing LABSWE source term:");
                for(unsigned long i(0); i < data.h->size(); ++i)
                {
                    result[i] = - g * (*data.h)[i] * dbx[i];
                }
            }
    };
}
#endif
