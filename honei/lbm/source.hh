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


#ifndef LBM_GUARD_SOURCE_HH
#define LBM_GUARD_SOURCE_HH 1

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

    /**
     * \brief Simple source term for use with LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct Source<Tag_, lbm_applications::LABSWE, lbm_source_types::SIMPLE, lbm_source_schemes::BASIC>
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
            static void value(DenseMatrix<DT1_> & result, DenseMatrix<DT1_>& h, DenseMatrix<DT2_>& dbx, DT3_ g)
            {
                CONTEXT("When computing LABSWE source term:");
                for(unsigned long i(0); i < h.rows(); ++i)
                {
                    for(unsigned long j(0); j < h.columns(); ++j)
                    {
                        result(i,j) = - g * h(i,j) * dbx(i,j);
                    }
                }
            }
    };

    /**
     * \brief Simple source term for use with LABSWE.
     *
     * \ingroup grplbmoperations
     */
    template <typename Tag_>
    struct Source<Tag_, lbm_applications::LABSWE, lbm_source_types::CONSTANT, lbm_source_schemes::BASIC>
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
         * \param constant The constant to be applied.
         *
         */
        template<typename DT1_, typename DT2_>
            static void value(DenseMatrix<DT1_> & result, DT2_ constant)
            {
                CONTEXT("When computing LABSWE source term:");
                for(unsigned long i(0); i < result.rows(); ++i)
                {
                    for(unsigned long j(0); j < result.columns(); ++j)
                    {
                        result(i,j) = constant;
                    }
                }
            }
    };

}
#endif
