/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the LBM C++ library. LibLBM is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibLBM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */


#ifndef LIBLBM_GUARD_SOURCE_HH
#define LIBLBM_GUARD_SOURCE_HH 1

/**
 * \file
 * Implementation of source modules used by  LBM - (SWE) solvers.
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
              typename App_,
              typename SourceType_,
              typename SourceScheme_>
    struct Source
    {
    };

    template <typename Tag_>
    struct Source<Tag_, lbm_applications::LABSWE, lbm_source_types::SIMPLE, lbm_source_schemes::BASIC>
    {
        template<typename DT1_, typename DT2_, typename DT3_>
            static void value(DenseMatrix<DT1_> & result, DenseMatrix<DT1_>& h, DenseMatrix<DT2_>& dbx, DT3_ g)
            {
                CONTEXT("When computing LABSWE source term.");
                for(unsigned long i(0); i < h.rows(); ++i)
                {
                    for(unsigned long j(0); j < h.columns(); ++j)
                    {
                        result(i,j) = - g * h(i,j) * dbx(i,j);
                    }
                }
            }
    };
}
#endif