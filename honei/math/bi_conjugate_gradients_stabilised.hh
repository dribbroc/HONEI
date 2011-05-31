/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef MATH_GUARD_BI_CONJUGATE_GRADIENTS_STABILISED_HH
#define MATH_GUARD_BI_CONJUGATE_GRADIENTS_STABILISED_HH 1

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
#include <iostream>

namespace honei
{
    template <typename Tag_, typename Method_>
        struct PBiCGStab
        {
        };
    /**
     * \brief Solution of linear system with preconditioned BiCGStab
     *
     * \ingroup grpmatrixoperations
     * \ingroup grpvectoroperations
     */
    template <typename Tag_>
        struct PBiCGStab<Tag_, methods::VAR>
        {
            public:
            template <typename DT_>
            static void value(SparseMatrixELL<DT_> & system_matrix,
                                           DenseVector<DT_> & right_hand_side,
                                           DenseVector<DT_> & x,
                                           SparseMatrixELL<DT_> & precon,
                                           unsigned long max_iters,
                                           unsigned long & used_iters,
                                           DT_ eps_relative = 1e-8)
            {
                CONTEXT("When solving sparse ELL linear system with PBiCGStab: ");
#if (defined SOLVER_VERBOSE_L2 || defined SOLVER_VERBOSE_L3)
                std::cout << "Calling CG solver, preconditioning=VAR, datalayout=ELL" << std::endl;
#endif
                bool converged(false);
                DenseVector<DT_> r(right_hand_side.size(), DT_(0));
                Defect<Tag_>::value(r, right_hand_side, system_matrix, x);

                DenseVector<DT_> v(right_hand_side.size(), DT_(0));
                DenseVector<DT_> p(right_hand_side.size());
                DenseVector<DT_> p_hat(right_hand_side.size());
                DenseVector<DT_> s(right_hand_side.size());
                DenseVector<DT_> s_hat(right_hand_side.size());

                DenseVector<DT_> t(right_hand_side.size());

                DenseVector<DT_> r_tilde(r.size());
                copy<Tag_>(r, r_tilde);
                DT_ norm_initial(Norm<vnt_l_two, true, Tag_>::value(r));
                std::cout << "initial_norm: " << norm_initial << std::endl;

                DT_ rho_old(1), rho_new(1), alpha(1), beta(1), omega(1);

                for(unsigned long i(0) ; i < max_iters ; ++i)
                {
                    used_iters = i + 1;
                    //std::cout << "iter:" << i << std::endl;
                    rho_new = DotProduct<Tag_>::value(r_tilde, r);
                    //if(fabs(rho_new) <= std::numeric_limits<DT_>::epsilon())
                    if(fabs(rho_new)  <= std::numeric_limits<DT_>::epsilon())
                    {
                        std::cout << "Search vectors orthogonal!" << std::endl;
                        DT_ norm_failed(Norm<vnt_l_two, true, Tag_>::value(r));
                        std::cout << "norm: " << norm_failed << " at iteration " << i << std::endl;
                        converged = true;
                        break;
                    }

                    if(i == 0)
                    {
                        copy<Tag_>(r, p);
                    }
                    else
                    {
                        beta = (rho_new / rho_old) * (alpha / omega);
                        Scale<Tag_>::value(v , -omega);
                        Sum<Tag_>::value(v, p);
                        ScaledSum<Tag_>::value(p, r, v, beta);
                    }

                    Product<Tag_>::value(p_hat, precon, p);
                    Product<Tag_>::value(v, system_matrix, p_hat);
                    alpha = DotProduct<Tag_>::value(r_tilde, v);
                    alpha = fabs(alpha) >= std::numeric_limits<DT_>::epsilon() ? rho_new / alpha : rho_new / std::numeric_limits<DT_>::epsilon();
                    ScaledSum<Tag_>::value(s, r, v, -alpha);

                    //main convergence-control
                    DT_ norm(Norm<vnt_l_two, true, Tag_>::value(s));
                    if(norm / norm_initial <= eps_relative)
                    {
                        ScaledSum<Tag_>::value(x, p_hat, alpha);
                        std::cout << "1st Converged with norm" << norm << " in iteration " << i << std::endl;
                        break;
                    }

                    Product<Tag_>::value(s_hat, precon, s);
                    Product<Tag_>::value(t, system_matrix, s_hat);

                    omega = DotProduct<Tag_>::value(t, t);
                    omega = fabs(omega) >= std::numeric_limits<DT_>::epsilon() ? DotProduct<Tag_>::value(t, s) / omega : DotProduct<Tag_>::value(t, s) / std::numeric_limits<DT_>::epsilon();
                    Scale<Tag_>::value(s_hat, omega);
                    Scale<Tag_>::value(p_hat, alpha);
                    Sum<Tag_>::value(s_hat, p_hat);
                    Sum<Tag_>::value(x, s_hat);
                    ScaledSum<Tag_>::value(r, s, t, -omega);

                    norm = Norm<vnt_l_two, true, Tag_>::value(r);

                    //std::cout << "rNORM: " << norm << std::endl;

                    if(norm / norm_initial <= eps_relative)
                    {
                        std::cout << "2nd Converged with norm" << norm << " in iteration " << i << std::endl;
                        converged = true;
                        break;
                    }

                    if(fabs(omega) < std::numeric_limits<DT_>::epsilon())
                    {
                        std::cout << "Omega near zero, aborting." << std::endl;
                        converged = true;
                        break;
                    }

                    rho_old = rho_new;
                }
                if(!converged)
                    std::cout << "No convergence after " << max_iters << " iterations. Norm: " << Norm<vnt_l_two, true, Tag_>::value(r) << std::endl;

            }
        };
}
#endif
