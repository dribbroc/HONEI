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

#pragma once
#ifndef MATH_GUARD_CG_HH
#define MATH_GUARD_CG_HH 1

#include <honei/util/tags.hh>
#include <honei/math/vectorpool.hh>
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
#include <math.h>

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

                //std::cout << "MAX_ITERS IN CG" << max_iters << std::endl;

                VectorType_ p(b.size());
                VectorType_ r(b.size());
                VectorType_ v(b.size());

                DT_ alpha, alpha_old, lambda, initial_defect, current_defect(0);
                unsigned long iterations(0);

                Defect<Tag_>::value(r, b, A, x);
                copy<Tag_>(r, p);

                DT_ temp;

                alpha = Norm<vnt_l_two, false, Tag_>::value(r);
                initial_defect = sqrt(alpha);
                while(iterations < max_iters)
                {
                    Product<Tag_>::value(v, A, p);
                    temp = DotProduct<Tag_>::value(v, p);
                    lambda = alpha / (std::abs(temp) > std::numeric_limits<DT_>::epsilon() ? temp : std::numeric_limits<DT_>::epsilon());
                    ScaledSum<Tag_>::value(x, p, lambda);
                    DT_ mlambda(-lambda);
                    ScaledSum<Tag_>::value(r, v, mlambda);
                    alpha_old = alpha;
                    alpha = Norm<vnt_l_two, false, Tag_>::value(r);

                    DT_ talpha(alpha / (std::abs(alpha_old) > std::numeric_limits<DT_>::epsilon() ? alpha_old : std::numeric_limits<DT_>::epsilon()));
                    Scale<Tag_>::value(p, talpha);
                    Sum<Tag_>::value(p, r);

                    ++iterations;

                    current_defect = sqrt(alpha);
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
                LOGMESSAGE(lc_solver, "CG(NONE) finished in " + stringify(used_iters) + " iterations with defect " + stringify(current_defect));
                PROFILER_STOP("CGSolver NONE");
                //std::cout << used_iters << " " << max_iters << " " << iterations << std::endl;
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
                    lambda = alpha_new / (std::abs(temp) > std::numeric_limits<DT_>::epsilon() ? temp : (std::numeric_limits<DT_>::epsilon()* (temp/std::abs(temp))));

                    ++iterations;
                    ScaledSum<Tag_>::value(x, p, lambda);

                    DT_ mlambda(-lambda);
                    ScaledSum<Tag_>::value(r, v, mlambda);

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

                    DT_ talpha_new(alpha_new / (std::abs(alpha) > std::numeric_limits<DT_>::epsilon() ? alpha : (std::numeric_limits<DT_>::epsilon() * (alpha / std::abs(alpha)))));
                    Scale<Tag_>::value(p, talpha_new);
                    Sum<Tag_>::value(p, z);
                }

                PROFILER_STOP("CGSolver VAR");
                //std::cout << used_iters << " " << max_iters << " " << iterations << std::endl;
                return x;
            }
    };

    /**
     * \brief Smoothing with PCG. Variable preconditioning.
     *
     * \ingroup grpmatrixoperations
     * \ingroup grpvectoroperations
     */
    template <typename Tag_>
    struct CGSmoother
    {
        public:
            static const unsigned long NUM_TEMPVECS = 4;

            template<typename MatrixType_, typename VectorType_, typename PreconContType_>
            static inline VectorType_ & value(MatrixType_ & A,
                                     PreconContType_ & P,
                                     VectorType_ & b,
                                     VectorType_ & x,
                                     std::vector<VectorType_> & temp_vecs,
                                     unsigned long max_iters)
            {
                CONTEXT("When smoothing with CG :");
                PROFILER_START("CGSmoother");

                VectorType_ p(temp_vecs.at(0));
                VectorType_ r(temp_vecs.at(1));
                VectorType_ v(temp_vecs.at(2));
                VectorType_ z(temp_vecs.at(3));

                typename VectorType_::DataType alpha, alpha_new, lambda, initial_defect;
                unsigned long iterations(0);

                Defect<Tag_>::value(r, b, A, x);

                Product<Tag_>::value(p, P, r);

                initial_defect = Norm<vnt_l_two, true, Tag_>::value(r);

                alpha_new = DotProduct<Tag_>::value(r, p);

                typename VectorType_::DataType temp;

                while(iterations < max_iters)
                {
                    Product<Tag_>::value(v, A, p);
                    temp = DotProduct<Tag_>::value(v, p);
                    lambda = alpha_new / (std::abs(temp) > std::numeric_limits<typename VectorType_::DataType>::epsilon() ? temp : (std::numeric_limits<typename VectorType_::DataType>::epsilon()* (temp/std::abs(temp))));

                    ++iterations;
                    ScaledSum<Tag_>::value(x, p, lambda);

                    typename VectorType_::DataType mlambda(-lambda);
                    ScaledSum<Tag_>::value(r, v, mlambda);

                    if(iterations == max_iters)
                    {
                        break;
                    }

                    Product<Tag_>::value(z, P, r);

                    alpha = alpha_new;

                    alpha_new = DotProduct<Tag_>::value(r, z);

                    typename VectorType_::DataType talpha_new(alpha_new / (std::abs(alpha) > std::numeric_limits<typename VectorType_::DataType>::epsilon() ? alpha : (std::numeric_limits<typename VectorType_::DataType>::epsilon() * (alpha / std::abs(alpha)))));
                    Scale<Tag_>::value(p, talpha_new);
                    Sum<Tag_>::value(p, z);
                }

                PROFILER_STOP("CGSmoother VAR");
                return x;
            }
    };
}
#endif
