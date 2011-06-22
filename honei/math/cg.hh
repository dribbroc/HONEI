/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the MATH C++ library. LibMath is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibMath is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef MATH_GUARD_CG_HH
#define MATH_GUARD_CG_HH 1

#include <honei/util/tags.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/dense_vector.hh>
#include <honei/math/defect.hh>
#include <honei/la/product.hh>
#include <honei/la/sum.hh>
#include <honei/la/difference.hh>
#include <honei/la/dot_product.hh>
#include <honei/la/scale.hh>
#include <honei/la/scaled_sum.hh>
#include <honei/la/norm.hh>
#include <honei/math/methods.hh>
#include <honei/la/element_product.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/la/algorithm.hh>
#include <honei/util/profiler.hh>
#include <iostream>

namespace honei
{
    /**
     * \brief Solution of linear system with CG.
     *
     * \ingroup grpmatrixoperations
     * \ingroup grpvectoroperations
     */
    template <typename Tag_, typename PreconType_>
    struct CG
    {
    };

    /**
     * \brief Solution of linear system with CG. No preconditioning.
     *
     * \ingroup grpmatrixoperations
     * \ingroup grpvectoroperations
     */
    template <typename Tag_>
    struct CG<Tag_, methods::NONE>
    {
        public:
            template<typename DT_, typename MatrixType_, typename VectorType_>
            static inline VectorType_ & value(MatrixType_ & A,
                                              VectorType_ & b,
                                              VectorType_ & x,
                                              unsigned long max_iters,
                                              unsigned long & used_iters,
                                              DT_ eps_relative = 1e-8)
            {
                CONTEXT("When solving linear system with CG :");
                PROFILER_START("CGSolver NONE");

                VectorType_ p(b.size());
                VectorType_ r(b.size());
                VectorType_ v(b.size());

                DT_ alpha, alpha_old, lambda, initial_defect, current_defect(0);
                unsigned long iterations(0);

                Defect<Tag_>::value(r, b, A, x);
                copy<Tag_>(r, p);

                alpha = Norm<vnt_l_two, false, Tag_>::value(r);
                initial_defect = sqrt(alpha);
                while(iterations < max_iters)
                {
                    Product<Tag_>::value(v, A, p);
                    lambda = alpha / DotProduct<Tag_>::value(v, p);
                    ScaledSum<Tag_>::value(x, p, lambda);
                    ScaledSum<Tag_>::value(r, v, -lambda);
                    alpha_old = alpha;
                    alpha = Norm<vnt_l_two, false, Tag_>::value(r);

                    Scale<Tag_>::value(p, alpha / alpha_old);
                    Sum<Tag_>::value(p, r);

                    ++iterations;

                    current_defect = sqrt(alpha);
                    if(current_defect < eps_relative * initial_defect)
                    {
                        used_iters = iterations + 1;
                        break;
                    }
                    if(current_defect < eps_relative)
                    {
                        used_iters = iterations + 1;
                        break;
                    }
                    if(iterations == max_iters)
                    {
                        used_iters = iterations;
                    }
                }
                LOGMESSAGE(lc_solver, "CG(NONE) finished in " + stringify(used_iters) + " iterations with defect " + stringify(current_defect));
                PROFILER_STOP("CGSolver NONE");
                return x;
            }
    };
}
#endif
