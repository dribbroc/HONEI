/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#pragma once
#ifndef LIBMATH_GUARD_PRECONDITIONING_HH
#define LIBMATH_GUARD_PRECONDITIONING_HH 1

#include <methods.hh>
#include <honei/util/attributes.hh>

namespace honei
{
    /**
     * \brief Preconditioning prototype
     *
     * \ingroup grpmatrixoperations
     */
    template <typename Tag_, typename Method_>
    struct Preconditioning
    {
    };

    template <typename Tag_>
    struct Preconditioning<Tag_, methods::NONE>
    {
        public:
            /**
            * \brief Computes approximate inverse
            *
            * \param target The approximate inverse of A.
            * \param source A.
            *
            */
            template <typename MatrixType_>
            static MatrixType_ value(HONEI_UNUSED MatrixType_ & source)
            {
                CONTEXT("When computing approximate inverse for preconditioning: ");
                //Does nothing - inverse has to be passed to this routine
            }
    };

    template <typename Tag_>
    struct Preconditioning<Tag_, methods::JAC>
    {
        public:
            /**
            * \brief Computes approximate inverse
            *
            * \param target The approximate inverse of A.
            * \param source A.
            *
            */
            template <typename MatrixType_>
            static MatrixType_ value(HONEI_UNUSED MatrixType_ & source)
            {
                CONTEXT("When computing approximate inverse for preconditioning with Jacobi: ");
                //Does nothing - inverse has to be passed to this routine
            }

            template <typename MatrixType_>
            static DenseVector<typename MatrixType_::DataType> value(HONEI_UNUSED MatrixType_ & source)
            {
                CONTEXT("When computing approximate inverse for preconditioning with Jacobi: ");
                //Does nothing - inverse has to be passed to this routine
            }
    };
}
#endif
