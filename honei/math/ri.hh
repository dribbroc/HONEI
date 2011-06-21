/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Markus Geveler <apryde@gmx.de>
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

#ifndef MATH_GUARD_RI_HH
#define MATH_GUARD_RI_HH 1

#include <honei/la/product.hh>
#include <honei/math/defect.hh>
#include <honei/la/sum.hh>
#include <honei/util/profiler.hh>

namespace honei
{

    /**
     * \brief  Smoothing with Richardson.
     *
     * \ingroup grpmatrixoperations
     * \ingroup grpvectoroperations
     */
    template <typename Tag_>
    struct RISmoother
    {
        public:
            /**
            * \brief Returns smoothed vector stemming from linear system with the preconditioned and damped Richardson iterative method.
            *
            */
            template <typename MatrixType_,
                      typename VectorType_,
                      typename PreconContType_>
            static void value(MatrixType_ & A,
                              PreconContType_ & P,
                              VectorType_ & b,
                              VectorType_ & x,
                              VectorType_ & temp_0,
                              VectorType_ & temp_1,
                              unsigned long max_iters)
            {
                CONTEXT("When smoothing with Richardson: ");
                PROFILER_START("RISmoother");

                for(unsigned long i(0) ; i < max_iters ; ++i)
                {
                    Defect<Tag_>::value(temp_0, b, A, x);
                    Product<Tag_>::value(temp_1, P, temp_0);
                    Sum<Tag_>::value(x, temp_1);
                }

                PROFILER_STOP("RISmoother");
            }
    };
}

#endif
