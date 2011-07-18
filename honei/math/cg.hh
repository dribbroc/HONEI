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
            template<typename DT_, typename MatrixType_, typename VectorType_, typename PreconContType_>
            static inline VectorType_ & value(MatrixType_ & A,
                                              HONEI_UNUSED PreconContType_ & P,
                                              VectorType_ & b,
                                              VectorType_ & x,
                                              unsigned long max_iters,
                                              unsigned long & used_iters,
                                              DT_ eps_relative = 1e-8)
            {
                CONTEXT("When solving linear system with CG :");
                PROFILER_START("CGSolver NONE");

                std::cout << "MAX_ITERS IN CG" << max_iters << std::endl;

                VectorType_ p(b.size());
                VectorType_ r(b.size());
                VectorType_ v(b.size());

                DT_ alpha, alpha_old, lambda, initial_defect;
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

                    DT_ current_defect(sqrt(alpha));
                    if(current_defect < eps_relative * initial_defect)
                    {
                        used_iters = iterations;
                        break;
                    }
                    if(current_defect < eps_relative)
                    {
                        used_iters = iterations;
                        break;
                    }
                    if(iterations == max_iters)
                    {
                        used_iters = iterations;
                    }
                }

                PROFILER_STOP("CGSolver NONE");
                std::cout << used_iters << " " << max_iters << " " << iterations << std::endl;
                return x;
            }
    };

    /**
     * \brief Solution of linear system with CG. Variable preconditioning.
     *
     * \ingroup grpmatrixoperations
     * \ingroup grpvectoroperations
     */
    template <typename Tag_>
    struct CG<Tag_, methods::VAR>
    {
        public:
            template<typename DT_, typename PreconContType_, typename MatrixType_, typename VectorType_>
            static inline VectorType_ & value(MatrixType_ & A,
                                              PreconContType_ & P,
                                              VectorType_ & b,
                                              VectorType_ & x,
                                              unsigned long max_iters,
                                              unsigned long & used_iters,
                                              DT_ eps_relative = 1e-8)
            {
                CONTEXT("When solving linear system with CG :");
                PROFILER_START("CGSolver VAR");

                VectorType_ p(b.size());
                VectorType_ r(b.size());
                VectorType_ v(b.size());
                VectorType_ z(b.size());

                DT_ alpha, alpha_new, lambda, initial_defect;
                unsigned long iterations(0);

                Defect<Tag_>::value(r, b, A, x);

                Product<Tag_>::value(p, P, r);

                initial_defect = Norm<vnt_l_two, true, Tag_>::value(r);

                alpha_new = DotProduct<Tag_>::value(r, p);

                DT_ temp;

                while(iterations < max_iters)
                {
                    Product<Tag_>::value(v, A, p);

                    temp = DotProduct<Tag_>::value(v, p);
                    lambda = alpha_new / (fabs(temp) >= std::numeric_limits<DT_>::epsilon() ? temp : std::numeric_limits<DT_>::epsilon());

                    ++iterations;
                    ScaledSum<Tag_>::value(x, p, lambda);

                    ScaledSum<Tag_>::value(r, v, -lambda);

                    DT_ current_defect(Norm<vnt_l_two, true, Tag_>::value(r));
                    if(current_defect < eps_relative * initial_defect)
                    {
                        used_iters = iterations;
                        break;
                    }
                    if(current_defect < eps_relative)
                    {
                        used_iters = iterations;
                        break;
                    }
                    if(iterations == max_iters)
                    {
                        used_iters = iterations;
                        break;
                    }

                    Product<Tag_>::value(z, P, r);

                    alpha = alpha_new;

                    alpha_new = DotProduct<Tag_>::value(r, z);

                    Scale<Tag_>::value(p, alpha_new / (fabs(alpha) >= std::numeric_limits<DT_>::epsilon() ? alpha : std::numeric_limits<DT_>::epsilon()) );
                    Sum<Tag_>::value(p, z);
                }

                PROFILER_STOP("CGSolver VAR");
                std::cout << used_iters << " " << max_iters << " " << iterations << std::endl;
                return x;
            }
    };
}
#endif
