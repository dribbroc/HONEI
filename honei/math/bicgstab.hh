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
#include <honei/math/vectorpool.hh>
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
    struct BiCGStabSolver
    {
    };

    /**
     * \brief Solution of linear system with BiCGStab. Variable preconditioning.
     *
     * \ingroup grpmatrixoperations
     * \ingroup grpvectoroperations
     */
    template <typename Tag_>
    struct BiCGStabSolver<Tag_, methods::VAR>
    {
        public:
            template<typename VT_>
            static void vectorpool(unsigned long size, std::vector<VT_> & result)
            {
                result = honei::create_vectorpool<VT_>(10, size);
            }

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

                DT_ defnorm, defnorm_0, defnorm_00(1e14);
                //DT_ kappa = 1.0;
                unsigned long iter = 0;
                DT_ rho_tilde, rho_tilde_old, alpha_tilde, omega_tilde, beta_tilde, gamma_tilde;
                DT_ nrm_r_tilde_0, nrm_tilde_00;
                bool early_exit = 0;
                bool restarted = false;
                bool converged = 0;

                VectorType_ r(b.size());
                VectorType_ r_tilde(b.size());
                VectorType_ r_tilde_0(b.size());
                VectorType_ p_tilde(b.size());
                VectorType_ v(b.size());
                VectorType_ v_tilde(b.size());
                VectorType_ s(b.size());
                VectorType_ s_tilde(b.size());
                VectorType_ t(b.size());
                VectorType_ t_tilde(b.size());

                do
                {
                    Defect<Tag_>::value(r, b, A, x);
                    defnorm_0 = Norm<vnt_l_two, false, Tag_>::value(r);
                    defnorm = defnorm_0;
                    Product<Tag_>::value(r_tilde_0, P, r);
                    nrm_r_tilde_0 = Norm<vnt_l_two, false, Tag_>::value(r_tilde_0);

                    if (restarted == false)
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

                        if (std::abs(gamma_tilde) < std::abs(rho_tilde)*1e-14)
                        {
                            restarted = true;
                            //std::cout << "Breakpoint 1" << std::endl;
                            break;
                        }

                        alpha_tilde = rho_tilde / gamma_tilde;

                        if ((std::abs(alpha_tilde) * Norm<vnt_l_two, false, Tag_>::value(v_tilde)) / defnorm < 1e-5)
                        {
                            restarted = true;;
                            //std::cout << "Breakpoint 2" << std::endl;
                            // \TODO warum ist das break hier nicht aktiv?
                            //break;
                        }

                        DT_ malpha_tilde(-alpha_tilde);
                        ScaledSum<Tag_>::value(s, r, v, malpha_tilde);

                        defnorm = Norm<vnt_l_two, false, Tag_>::value(s);
                        if (defnorm < eps_relative * defnorm_00)
                        {
                            ScaledSum<Tag_>::value(x, p_tilde, alpha_tilde);

                            early_exit = 1;
                            converged = 1;
                            //std::cout << "Breakpoint 3 (converged)" << std::endl;
                            break;
                        }
                        ScaledSum<Tag_>::value(s_tilde, r_tilde, v_tilde, malpha_tilde);

                        Product<Tag_>::value(t, A, s_tilde);

                        Product<Tag_>::value(t_tilde, P, t);

                        gamma_tilde = DotProduct<Tag_>::value(t_tilde, t_tilde);
                        omega_tilde = DotProduct<Tag_>::value(t_tilde, s_tilde);

                        if (std::abs(gamma_tilde) < std::abs(omega_tilde) * 1e-14)
                        {
                            restarted = true;
                            //std::cout << "Breakpoint 4" << std::endl;
                            break;
                        }
                        omega_tilde = omega_tilde / gamma_tilde;

                        ScaledSum<Tag_>::value(x, s_tilde, omega_tilde);
                        ScaledSum<Tag_>::value(x, p_tilde, alpha_tilde);

                        DT_ momega_tilde(-omega_tilde);
                        ScaledSum<Tag_>::value(r, s, t, momega_tilde);

                        defnorm = Norm<vnt_l_two, false, Tag_>::value(r);
                        if (defnorm < eps_relative * defnorm_00)
                        {
                            converged = 1;
                            //std::cout << "Breakpoint 5 (converged)" << std::endl;
                            break;
                        }

                        ScaledSum<Tag_>::value(r_tilde, s_tilde, t_tilde, momega_tilde);

                        rho_tilde_old = rho_tilde;
                        rho_tilde = DotProduct<Tag_>::value(r_tilde, r_tilde_0);

                        beta_tilde = (alpha_tilde / omega_tilde) * (rho_tilde / rho_tilde_old);

                        ScaledSum<Tag_>::value(p_tilde, v_tilde, momega_tilde);
                        Scale<Tag_>::value(p_tilde, beta_tilde);
                        Sum<Tag_>::value(p_tilde, r_tilde);

                    } while (iter <= max_iters);

                } while (converged == 0 && iter < max_iters);


                used_iters = iter + 1;
                //std::cout << "Norm: " << defnorm << std::endl;
                LOGMESSAGE(lc_solver, "BiCgStab(VAR) finished in " + stringify(used_iters) + " iterations with defect " + stringify(defnorm));
                PROFILER_STOP("CGSolver NONE");

                return x;
            }
    };

    /**
     * \brief Smoothing with BiCGStab. Variable preconditioning.
     *
     * \ingroup grpmatrixoperations
     * \ingroup grpvectoroperations
     */
    template <typename Tag_>
    struct BiCGStabSmoother
    {
        public:
            template<typename VT_>
            static void vectorpool(unsigned long size, std::vector<VT_> & result)
            {
                result = honei::create_vectorpool<VT_>(10, size);
            }

            template<typename MatrixType_, typename VectorType_, typename PreconContType_>
            static inline VectorType_ & value(MatrixType_ & A,
                                              PreconContType_ & P,
                                              VectorType_ & b,
                                              VectorType_ & x,
                                              std::vector<VectorType_> & temp_vecs,
                                              unsigned long max_iters)
            {
                CONTEXT("When smoothing with BiCGStab :");
                PROFILER_START("Smoothing BiCGStab");

                typename VectorType_::DataType defnorm, defnorm_0, defnorm_00(1e14);
                //typename VectorType_::DataType kappa = 1.0;
                unsigned long iter = 0;
                typename VectorType_::DataType rho_tilde, rho_tilde_old, alpha_tilde, omega_tilde, beta_tilde, gamma_tilde;
                typename VectorType_::DataType nrm_r_tilde_0, nrm_tilde_00;
                bool restarted = false;

                /*VectorType_ r(b.size());
                VectorType_ r_tilde(b.size());
                VectorType_ r_tilde_0(b.size());
                VectorType_ p_tilde(b.size());
                VectorType_ v(b.size());
                VectorType_ v_tilde(b.size());
                VectorType_ s(b.size());
                VectorType_ s_tilde(b.size());
                VectorType_ t(b.size());
                VectorType_ t_tilde(b.size());*/

                VectorType_ r(temp_vecs.at(0));
                VectorType_ r_tilde(temp_vecs.at(1));
                VectorType_ r_tilde_0(temp_vecs.at(2));
                VectorType_ p_tilde(temp_vecs.at(3));
                VectorType_ v(temp_vecs.at(4));
                VectorType_ v_tilde(temp_vecs.at(5));
                VectorType_ s(temp_vecs.at(6));
                VectorType_ s_tilde(temp_vecs.at(7));
                VectorType_ t(temp_vecs.at(8));
                VectorType_ t_tilde(temp_vecs.at(9));

                ///TODO optimise -> kill convergence control
                do
                {
                    Defect<Tag_>::value(r, b, A, x);
                    defnorm_0 = Norm<vnt_l_two, false, Tag_>::value(r);
                    defnorm = defnorm_0;
                    Product<Tag_>::value(r_tilde_0, P, r);
                    nrm_r_tilde_0 = Norm<vnt_l_two, false, Tag_>::value(r_tilde_0);

                    if (restarted == false)
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

                        if (std::abs(gamma_tilde) < std::abs(rho_tilde)*1e-14)
                        {
                            restarted = true;
                            //std::cout << "Breakpoint 1" << std::endl;
                            break;
                        }

                        alpha_tilde = rho_tilde / gamma_tilde;

                        if ((std::abs(alpha_tilde) * Norm<vnt_l_two, false, Tag_>::value(v_tilde)) / defnorm < 1e-5)
                        {
                            restarted = true;;
                            //std::cout << "Breakpoint 2" << std::endl;
                            // \TODO warum ist das break hier nicht aktiv?
                            //break;
                        }

                        typename VectorType_::DataType malpha_tilde(-alpha_tilde);
                        ScaledSum<Tag_>::value(s, r, v, malpha_tilde);

                        defnorm = Norm<vnt_l_two, false, Tag_>::value(s);
                        ScaledSum<Tag_>::value(s_tilde, r_tilde, v_tilde, malpha_tilde);

                        Product<Tag_>::value(t, A, s_tilde);

                        Product<Tag_>::value(t_tilde, P, t);

                        gamma_tilde = DotProduct<Tag_>::value(t_tilde, t_tilde);
                        omega_tilde = DotProduct<Tag_>::value(t_tilde, s_tilde);

                        if (std::abs(gamma_tilde) < std::abs(omega_tilde) * 1e-14)
                        {
                            restarted = true;
                            //std::cout << "Breakpoint 4" << std::endl;
                            break;
                        }
                        omega_tilde = omega_tilde / gamma_tilde;

                        ScaledSum<Tag_>::value(x, s_tilde, omega_tilde);
                        ScaledSum<Tag_>::value(x, p_tilde, alpha_tilde);

                        typename VectorType_::DataType momega_tilde(-omega_tilde);
                        ScaledSum<Tag_>::value(r, s, t, momega_tilde);

                        defnorm = Norm<vnt_l_two, false, Tag_>::value(r);

                        ScaledSum<Tag_>::value(r_tilde, s_tilde, t_tilde, momega_tilde);

                        rho_tilde_old = rho_tilde;
                        rho_tilde = DotProduct<Tag_>::value(r_tilde, r_tilde_0);

                        beta_tilde = (alpha_tilde / omega_tilde) * (rho_tilde / rho_tilde_old);

                        ScaledSum<Tag_>::value(p_tilde, v_tilde, momega_tilde);
                        Scale<Tag_>::value(p_tilde, beta_tilde);
                        Sum<Tag_>::value(p_tilde, r_tilde);

                    } while (iter <= max_iters);

                } while (iter < max_iters);


                //std::cout << "Norm: " << defnorm << std::endl;
                LOGMESSAGE(lc_solver, "BiCgStab(VAR) smoother finished in " + stringify(iters + 1) + " iterations with defect " + stringify(defnorm));
                PROFILER_STOP("BiCGStabSmoother");

                return x;
            }
    };
}
#endif
