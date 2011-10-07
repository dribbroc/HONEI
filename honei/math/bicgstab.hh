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
#ifndef MATH_GUARD_BICGSTAB_HH
#define MATH_GUARD_BICGSTAB_HH 1

#include <honei/util/tags.hh>
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
#include <honei/la/algorithm.hh>
#include <honei/util/profiler.hh>
#include <iostream>
#include <math.h>

namespace honei
{
    /**
     * \brief Solution of linear system with BiCGStab.
     *
     * \ingroup grpmatrixoperations
     * \ingroup grpvectoroperations
     */
    template <typename Tag_, typename PreconType_>
    struct BiCGStab
    {
    };

    /**
     * \brief Solution of linear system with BiCGStab. Variable preconditioning.
     *
     * \ingroup grpmatrixoperations
     * \ingroup grpvectoroperations
     */
    template <typename Tag_>
    struct BiCGStab<Tag_, methods::VAR>
    {
        public:
            template<typename DT_, typename MatrixType_, typename VectorType_, typename PreconContType_>
            static inline VectorType_ & value(MatrixType_ & A,
                                              PreconContType_ & P,
                                              VectorType_ & b,
                                              VectorType_ & x,
                                              unsigned long max_iters,
                                              unsigned long & used_iters,
                                              DT_ eps_relative = 1e-8)
            {
                CONTEXT("When solving linear system with BiCGStab :");
                PROFILER_START("BiCGStab VAR");

                double defnorm, defnorm_0, defnorm_00(1e14);
                //double kappa = 1.0;
                unsigned long iter = 0;
                double rho_tilde, rho_tilde_old, alpha_tilde, omega_tilde, beta_tilde, gamma_tilde;
                double nrm_r_tilde_0, nrm_tilde_00;
                bool early_exit = 0;
                bool num_restarts = 0;
                bool converged = 0;
                bool max_restarts = 50;

                DenseVector<DT_> r(b.size());
                DenseVector<DT_> r_tilde(b.size());
                DenseVector<DT_> r_tilde_0(b.size());
                DenseVector<DT_> p_tilde(b.size());
                DenseVector<DT_> v(b.size());
                DenseVector<DT_> v_tilde(b.size());
                DenseVector<DT_> s(b.size());
                DenseVector<DT_> s_tilde(b.size());
                DenseVector<DT_> t(b.size());
                DenseVector<DT_> t_tilde(b.size());

                do
                {
                    Defect<Tag_>::value(r, b, A, x);
                    defnorm_0 = Norm<vnt_l_two, false, Tag_>::value(r);
                    defnorm = defnorm_0;
                    Product<Tag_>::value(r_tilde_0, P, r);
                    nrm_r_tilde_0 = Norm<vnt_l_two, false, Tag_>::value(r_tilde_0);

                    if (num_restarts == 0)
                    {
                        defnorm_00 = defnorm_0;
                        nrm_tilde_00 = nrm_r_tilde_0;
                    }
                    copy<Tag_>(r_tilde_0, r_tilde);
                    copy<Tag_>(r_tilde_0, p_tilde);
                    rho_tilde = DotProduct<Tag_>::value(r_tilde_0, r_tilde_0);

                    // main BiCGStab loop
                    do
                    {
                        iter = iter + 1;

                        Product<Tag_>::value(v, A, p_tilde);
                        Product<Tag_>::value(v_tilde, P, v);

                        gamma_tilde = DotProduct<Tag_>::value(v_tilde, r_tilde_0);

                        if (fabs(gamma_tilde) < fabs(rho_tilde)*1e-14)
                        {
                            num_restarts++;
                            //std::cout << "Breakpoint 1" << std::endl;
                            break;
                        }

                        alpha_tilde = rho_tilde / gamma_tilde;

                        if ((fabs(alpha_tilde) * Norm<vnt_l_two, false, Tag_>::value(v_tilde)) / defnorm < 1e-5)
                        {
                            num_restarts++;
                            //std::cout << "Breakpoint 2" << std::endl;
                            //break;
                        }

                        ScaledSum<Tag_>::value(s, r, v, -alpha_tilde);

                        defnorm = Norm<vnt_l_two, false, Tag_>::value(s);
                        if (defnorm < eps_relative * defnorm_00)
                        {
                            ScaledSum<Tag_>::value(x, p_tilde, alpha_tilde);

                            early_exit = 1;
                            converged = 1;
                            //std::cout << "Breakpoint 3 (converged)" << std::endl;
                            break;
                        }
                        ScaledSum<Tag_>::value(s_tilde, r_tilde, v_tilde, -alpha_tilde);

                        Product<Tag_>::value(t, A, s_tilde);

                        Product<Tag_>::value(t_tilde, P, t);

                        gamma_tilde = DotProduct<Tag_>::value(t_tilde, t_tilde);
                        omega_tilde = DotProduct<Tag_>::value(t_tilde, s_tilde);

                        if (fabs(gamma_tilde) < fabs(omega_tilde) * 1e-14)
                        {
                            num_restarts++;
                            //std::cout << "Breakpoint 4" << std::endl;
                            break;
                        }
                        omega_tilde = omega_tilde / gamma_tilde;

                        ScaledSum<Tag_>::value(x, s_tilde, omega_tilde);
                        ScaledSum<Tag_>::value(x, p_tilde, alpha_tilde);

                        ScaledSum<Tag_>::value(r, s, t, -omega_tilde);

                        defnorm = Norm<vnt_l_two, false, Tag_>::value(r);
                        if (defnorm < eps_relative * defnorm_00)
                        {
                            converged = 1;
                            //std::cout << "Breakpoint 5 (converged)" << std::endl;
                            break;
                        }

                        ScaledSum<Tag_>::value(r_tilde, s_tilde, t_tilde, -omega_tilde);

                        rho_tilde_old = rho_tilde;
                        rho_tilde = DotProduct<Tag_>::value(r_tilde, r_tilde_0);

                        beta_tilde = (alpha_tilde / omega_tilde) * (rho_tilde / rho_tilde_old);

                        ScaledSum<Tag_>::value(p_tilde, v_tilde, -omega_tilde);
                        Scale<Tag_>::value(p_tilde, beta_tilde);
                        Sum<Tag_>::value(p_tilde, r_tilde);

                    } while (iter <= max_iters);

                } while (num_restarts < max_restarts  && converged == 0 && iter < max_iters);


                used_iters = iter + 1;
                //std::cout << "Norm: " << defnorm << std::endl;
                LOGMESSAGE(lc_solver, "BiCgStab(VAR) finished in " + stringify(used_iters) + " iterations with defect " + stringify(defnorm));
                PROFILER_STOP("CGSolver NONE");

                return x;
            }
    };
}
#endif
