/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
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

#ifndef LIBMATH_GUARD_CONJUGATE_GRADIENTS_HH
#define LIBMATH_GUARD_CONJUGATE_GRADIENTS_HH 1

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

//#define SOLVER_VERBOSE_L2 1
//#define SOLVER_VERBOSE_L1 1

using namespace methods;
namespace honei
{

    template<typename Tag_, typename PreCon_>
    struct ConjugateGradients
    {
    };

    /**
     * \brief Solution of linear system with CG.
     *
     * \ingroup grpmatrixoperations
     * \ingroup grpvectoroperations
     */
    template <typename Tag_>
    struct ConjugateGradients<Tag_, methods::NONE>
    {
        private:
            template<typename DT1_, typename DT2_>
            static inline void cg_kernel(DenseMatrix<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, DenseVector<DT1_> & former_gradient, DenseVector<DT1_> & former_result, DenseVector<DT1_> & utility)
            {
                ///Compute x_i+1:
                DT1_ upper = DotProduct<Tag_>::value(former_gradient, former_gradient);
                DenseVector<DT1_> energy = Product<Tag_>::value(system_matrix, utility);
                DT1_ lower = DotProduct<Tag_>::value(energy, utility);
                DenseVector<DT1_> u_c(utility.copy());
                if(fabs(lower) >= std::numeric_limits<DT1_>::epsilon())
                {
                    Scale<Tag_>::value(u_c, DT1_(upper/lower));
                    energy = u_c;
                }
                else
                {
                    Scale<Tag_>::value(u_c, DT1_(upper/std::numeric_limits<DT1_>::epsilon()));
                    energy = u_c;

                }
                Sum<Tag_>::value(energy, former_result);
                ///Compute new gradient
                DenseVector<DT1_> new_gradient = Product<Tag_>::value(system_matrix, energy);
                Difference<Tag_>::value(new_gradient, right_hand_side);

                ///Compute new utility
                upper = DotProduct<Tag_>::value(new_gradient, new_gradient);
                lower = DotProduct<Tag_>::value(former_gradient, former_gradient);
                if(fabs(lower) >= std::numeric_limits<DT1_>::epsilon())
                {
                    Scale<Tag_>::value(utility, DT1_(upper/lower));
                }
                else
                {
                    Scale<Tag_>::value(utility, DT1_(upper/std::numeric_limits<DT1_>::epsilon()));
                }
                Difference<Tag_>::value(utility, new_gradient);

                ///Finishing:
                former_gradient = new_gradient;//.copy();
                former_result = energy;//.copy();

            }

            template<typename DT1_, typename DT2_>
            static inline void cg_kernel(BandedMatrixQ1<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, DenseVector<DT1_> & former_gradient, DenseVector<DT1_> & former_result, DenseVector<DT1_> & utility)
            {
                ///Compute x_i+1: (in energy)
                DT1_ upper = DotProduct<Tag_>::value(former_gradient, former_gradient);
                DenseVector<DT1_> energy = Product<Tag_>::value(system_matrix, utility);
                DT1_ lower = DotProduct<Tag_>::value(energy, utility);
                DenseVector<DT1_> u_c(utility.size());
                copy<Tag_>(utility, u_c);
                if(fabs(lower) >= std::numeric_limits<DT1_>::epsilon())
                {
                    Scale<Tag_>::value(u_c, DT1_(upper/lower));
                    energy = u_c;
                }
                else
                {
                    Scale<Tag_>::value(u_c, DT1_(upper/std::numeric_limits<DT1_>::epsilon()));
                    energy = u_c;
                }

                Sum<Tag_>::value(energy, former_result);
                ///Compute new gradient
                DenseVector<DT1_> new_gradient = Product<Tag_>::value(system_matrix, energy);
                Difference<Tag_>::value(new_gradient, right_hand_side);

                ///Compute new utility
                upper = DotProduct<Tag_>::value(new_gradient, new_gradient);
                lower = DotProduct<Tag_>::value(former_gradient, former_gradient);
                if(fabs(lower) >= std::numeric_limits<DT1_>::epsilon())
                {
                    Scale<Tag_>::value(utility, DT1_(upper/lower));
                }
                else
                {
                    Scale<Tag_>::value(utility, DT1_(upper/std::numeric_limits<DT1_>::epsilon()));
                }
                Difference<Tag_>::value(utility, new_gradient);

                ///Finishing:
                former_gradient = new_gradient;
                former_result = energy;

            }

            //NEW ELL type:
            template<typename DT1_>
            static inline void cg_kernel(SparseMatrixELL<DT1_> & system_matrix, DenseVector<DT1_> & right_hand_side, DenseVector<DT1_> & x, DenseVector<DT1_> & r, DenseVector<DT1_> & p, DenseVector<DT1_> & v, unsigned long max_iters)
            {
                //MARKE
                DT1_ alpha, alpha_old, lambda, initial_defect;
                unsigned long iterations(0);
                //r=p=b-Ax
                Defect<Tag_>::value(r, right_hand_side, system_matrix, x);
                copy<Tag_>(r, p);
                alpha = Norm<vnt_l_two, false, Tag_>::value(r);
                initial_defect = sqrt(alpha);
                std::cout << "Initial defect NORM: " << initial_defect << std::endl;
                std::cout << "Initial alpha: " << alpha << std::endl;
                while(iterations < max_iters)
                {
                    Product<Tag_>::value(v, system_matrix, p);
                    lambda = alpha / DotProduct<Tag_>::value(v, p);
                    ScaledSum<Tag_>::value(x, p, lambda);
                    ScaledSum<Tag_>::value(r, v, -lambda);
                    alpha_old = alpha;
                    alpha = Norm<vnt_l_two, false, Tag_>::value(r);

                    Scale<Tag_>::value(p, alpha / alpha_old);
                    Sum<Tag_>::value(p, r);

                    ++iterations;

                    DT1_ current_defect(sqrt(alpha));
                    if(current_defect < 1e-08 * initial_defect)
                    {
                        std::cout << "Final defect NORM: " << current_defect << " after " << iterations << " iterations." << std::endl;
                        break;
                    }
                    if(current_defect < 1e-08)
                    {
                        std::cout << "ABORT. Final defect NORM: " << current_defect << " after " << iterations << " iterations." << std::endl;
                        break;
                    }
                    if(iterations == max_iters)
                        std::cout << "NO CONVERGENCE after " << max_iters << " iterations! Norm: " << current_defect << std::endl;
                }

            }

            template<typename DT1_, typename DT2_>
            static inline void cg_kernel(BandedMatrix<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, DenseVector<DT1_> & former_gradient, DenseVector<DT1_> & former_result, DenseVector<DT1_> & utility)
            {
                ///Compute x_i+1: (in energy)
                DT1_ upper = DotProduct<Tag_>::value(former_gradient, former_gradient);
                DenseVector<DT1_> energy = Product<Tag_>::value(system_matrix, utility);
                DT1_ lower = DotProduct<Tag_>::value(energy, utility);
                DenseVector<DT1_> u_c(utility.copy());
                if(fabs(lower) >= std::numeric_limits<DT1_>::epsilon())
                {
                    Scale<Tag_>::value(u_c, DT1_(upper/lower));
                    energy = u_c;
                }
                else
                {
                    Scale<Tag_>::value(u_c, DT1_(upper/std::numeric_limits<DT1_>::epsilon()));
                    energy = u_c;
                }
                Sum<Tag_>::value(energy, former_result);
                ///Compute new gradient
                DenseVector<DT1_> new_gradient = Product<Tag_>::value(system_matrix, energy);
                Difference<Tag_>::value(new_gradient, right_hand_side);

                ///Compute new utility
                upper = DotProduct<Tag_>::value(new_gradient, new_gradient);
                lower = DotProduct<Tag_>::value(former_gradient, former_gradient);
                if(fabs(lower) >= std::numeric_limits<DT1_>::epsilon())
                {
                    Scale<Tag_>::value(utility, DT1_(upper/lower));
                }
                else
                {
                    Scale<Tag_>::value(utility, DT1_(upper/std::numeric_limits<DT1_>::epsilon()));
                }
                Difference<Tag_>::value(utility, new_gradient);

                ///Finishing:
                former_gradient = new_gradient;//.copy();
                former_result = energy;//.copy();

            }

            template<typename DT1_, typename DT2_>
            static inline void cg_kernel(SparseMatrix<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, DenseVector<DT1_> & former_gradient, DenseVector<DT1_> & former_result, DenseVector<DT1_> & utility)
            {
                ///Compute x_i+1: (in energy)
                DT1_ upper = DotProduct<Tag_>::value(former_gradient, former_gradient);
                DenseVector<DT1_> energy = Product<Tag_>::value(system_matrix, utility);
                DT1_ lower = DotProduct<Tag_>::value(energy, utility);
                DenseVector<DT1_> u_c(utility.copy());
                if(fabs(lower) >= std::numeric_limits<DT1_>::epsilon())
                {
                    Scale<Tag_>::value(u_c, DT1_(upper/lower));
                    energy = u_c;
                }
                else
                {
                    Scale<Tag_>::value(u_c, DT1_(upper/std::numeric_limits<DT1_>::epsilon()));
                    energy = u_c;
                }
                Sum<Tag_>::value(energy, former_result);
                ///Compute new gradient
                DenseVector<DT1_> new_gradient = Product<Tag_>::value(system_matrix, energy);
                Difference<Tag_>::value(new_gradient, right_hand_side);

                ///Compute new utility
                upper = DotProduct<Tag_>::value(new_gradient, new_gradient);
                lower = DotProduct<Tag_>::value(former_gradient, former_gradient);
                if(fabs(lower) >= std::numeric_limits<DT1_>::epsilon())
                {
                    Scale<Tag_>::value(utility, DT1_(upper/lower));
                }
                else
                {
                    Scale<Tag_>::value(utility, DT1_(upper/std::numeric_limits<DT1_>::epsilon()));
                }
                Difference<Tag_>::value(utility, new_gradient);

                ///Finishing:
                former_gradient = new_gradient;//.copy();
                former_result = energy;//.copy();

            }

        public:
            /**
            * \brief Returns solution of linear system given by a DenseMatrix and a Vector.
            *
            * \param system_matrix The system matrix.
            * \param right_hand_side The right hand side of the system.
            * \param x The solution and initial guess.
            * \param iter_number The fixed number of iterations.
            *
           */
            template <typename DT1_, typename DT2_>
            static void value(DenseMatrix<DT1_> & system_matrix,
                    DenseVector<DT2_> & right_hand_side,
                    DenseVector<DT2_> & x,
                    unsigned long iter_number)
            {
                CONTEXT("When solving dense linear system with CG (fixed # iterations):");
#ifdef SOLVER_VERBOSE_L2
                std::cout << "Calling CG solver, preconditioning=NONE, datalayout=DENSE" << std::endl;
#endif

                DenseVector<DT1_> g = Product<Tag_>::value(system_matrix, x);
                Difference<Tag_>::value(g, right_hand_side);
                DenseVector<DT1_> g_c(g.copy());
                Scale<Tag_>::value(g_c, DT1_(-1.));
                DenseVector<DT1_> u(g_c.copy());
                for(unsigned long i = 0; i<iter_number; ++i)
                {
                    cg_kernel(system_matrix, right_hand_side, g, x, u);
                }
            }

            /**
            * \brief Returns solution of linear system given by a DenseMatrix and a Vector.
            *
            * \param system_matrix The system matrix.
            * \param right_hand_side The right hand side of the system.
            * \param x The solution and initial guess.
            * \param konv_rad The parameter for convergence control.
            *
            */
            template <typename DT1_, typename DT2_>
            static void value(DenseMatrix<DT1_> & system_matrix,
                    DenseVector<DT2_> & right_hand_side,
                    DenseVector<DT2_> & x,
                    double konv_rad)
            {
                CONTEXT("When solving dense linear system with CG (with given convergence parameter):");
#ifdef SOLVER_VERBOSE_L2
                std::cout << "Calling CG solver, preconditioning=NONE, datalayout=DENSE" << std::endl;
#endif

                DenseVector<DT1_> g = Product<Tag_>::value(system_matrix, x);
                Difference<Tag_>::value(g, right_hand_side);
                DenseVector<DT1_> g_c(g.copy());
                Scale<Tag_>::value(g_c, DT1_(-1.));
                DenseVector<DT1_> u(g_c.copy());
                DenseVector<DT1_> x_last(x.copy());

                DT1_ norm_x_last = DT1_(0);
                DT1_ norm_x = DT1_(1);

                while(fabs(norm_x - norm_x_last) > konv_rad)
                {

                    cg_kernel(system_matrix, right_hand_side, g, x, u);
                    norm_x = Norm<vnt_l_two, false, Tag_>::value(x);
                    norm_x_last = Norm<vnt_l_two, false, Tag_>::value(x_last);
                    x_last = x.copy();
                }
            }

            /**
            * \brief Returns solution of linear system given by a BandedMatrix and a Vector.
            *
            * \param system_matrix The system matrix.
            * \param right_hand_side The right hand side of the system.
            * \param x The solution and initial guess.
            * \param iter_number The fixed number of iterations.
            *
            */
            template <typename DT1_, typename DT2_>
            static void value(BandedMatrix<DT1_> & system_matrix,
                    DenseVector<DT2_> & right_hand_side,
                    DenseVector<DT2_> & x,
                    unsigned long iter_number)
            {
                CONTEXT("When solving banded linear system with CG (with fixed # of iterations):");
#ifdef SOLVER_VERBOSE_L2
                std::cout << "Calling CG solver, preconditioning=NONE, datalayout=BANDED" << std::endl;
#endif
                DenseVector<DT1_> g = Product<Tag_>::value(system_matrix, x);
                Difference<Tag_>::value(g, right_hand_side);
                DenseVector<DT1_> g_c(g.copy());
                Scale<Tag_>::value(g_c, DT1_(-1.));
                DenseVector<DT1_> u(g_c.copy());
                for(unsigned long i = 0; i<iter_number; ++i)
                {
                    cg_kernel(system_matrix, right_hand_side, g, x, u);
                }
            }

            /**
            * \brief Returns solution of linear system given by a BandedMatrix and a Vector.
            *
            * \param system_matrix The system matrix.
            * \param right_hand_side The right hand side of the system.
            * \param x The solution and initial guess.
            * \param konv_rad The parameter for convergence control.
            *
            */
            template <typename DT1_, typename DT2_>
            static void value(BandedMatrix<DT1_> & system_matrix,
                    DenseVector<DT2_> & right_hand_side,
                    DenseVector<DT2_> & x,
                    double konv_rad)
            {
                CONTEXT("When solving banded linear system with CG (with given convergence parameter):");
#ifdef SOLVER_VERBOSE_L2
                std::cout << "Calling CG solver, preconditioning=NONE, datalayout=BANDED" << std::endl;
#endif
                DenseVector<DT1_> g = Product<Tag_>::value(system_matrix, x);
                Difference<Tag_>::value(g, right_hand_side);
                DenseVector<DT1_> g_c(g.copy());
                Scale<Tag_>::value(g_c, DT1_(-1.));
                DenseVector<DT1_> u(g_c.copy());
                DenseVector<DT1_> x_last(x.copy());

                DT1_ norm_x_last = DT1_(0);
                DT1_ norm_x = DT1_(1);

                while(norm_x - norm_x_last > konv_rad)
                {

                    cg_kernel(system_matrix, right_hand_side, g, x, u);
                    norm_x = Norm<vnt_l_two, false, Tag_>::value(x);
                    norm_x_last = Norm<vnt_l_two, false, Tag_>::value(x_last);
                    x_last = x.copy();
                }
            }

            /**
            * \brief Returns solution of linear system given by a BandedMatrix and a Vector.
            *
            * \param system_matrix The system matrix.
            * \param right_hand_side The right hand side of the system.
            * \param x The solution and initial guess.
            * \param konv_rad The parameter for convergence control.
            *
            */
            template <typename DT1_, typename DT2_>
            static void value(BandedMatrixQ1<DT1_> & system_matrix,
                    DenseVector<DT2_> & right_hand_side,
                    DenseVector<DT2_> & x,
                    double konv_rad)
            {
                CONTEXT("When solving banded Q1 linear system with CG (with given convergence parameter):");
#ifdef SOLVER_VERBOSE_L2
                std::cout << "Calling CG solver, preconditioning=NONE, datalayout=Q1, MULTIGRID" << std::endl;
#endif

                fill<Tag_>(x, DT1_(0));
                DenseVector<DT1_> g = Product<Tag_>::value(system_matrix, x);
                Difference<Tag_>::value(g, right_hand_side);

                DenseVector<DT1_> g_c(g.size());
                copy<Tag_>(g, g_c);
                Scale<Tag_>::value(g_c, DT1_(-1.));
                DenseVector<DT1_> u(g_c.size());
                copy<Tag_>(g_c, u);
                DenseVector<DT1_> x_last(x.size());
                copy<Tag_>(x, x_last);
                DT1_ norm_x_last = DT1_(0);
                DT1_ norm_x = DT1_(1);

                while(norm_x - norm_x_last > konv_rad)
                {
                    cg_kernel(system_matrix, right_hand_side, g, x, u);
                    norm_x = Norm<vnt_l_two, false, Tag_>::value(x);
                    norm_x_last = Norm<vnt_l_two, false, Tag_>::value(x_last);
                    copy<Tag_>(x, x_last);
                }
                DenseVector<DT1_> r_def(Defect<Tag_>::value(right_hand_side, system_matrix, x));
                DT1_ norm_bla = Norm<vnt_l_two, false, Tag_>::value(r_def);
                std::cout << norm_bla << std::endl;
            }

            template <typename DT1_, typename DT2_>
            static void value(BandedMatrixQ1<DT1_> & system_matrix,
                    DenseVector<DT2_> & right_hand_side,
                    DenseVector<DT2_> & x,
                    unsigned long iters)
            {
                CONTEXT("When solving banded Q1 linear system with CG (with given convergence parameter):");
#ifdef SOLVER_VERBOSE_L2
                std::cout << "Calling CG solver, preconditioning=NONE, datalayout=Q1" << std::endl;
#endif

                fill<Tag_>(x, DT1_(0));
                DenseVector<DT1_> g = Product<Tag_>::value(system_matrix, x);
                Difference<Tag_>::value(g, right_hand_side);
                DenseVector<DT1_> g_c(g.size());
                copy<Tag_>(g, g_c);
                Scale<Tag_>::value(g_c, DT1_(-1.));
                DenseVector<DT1_> u(g_c.size());
                copy<Tag_>(g_c, u);
                DenseVector<DT1_> x_last(x.size());
                copy<Tag_>(x, x_last);
                //DT1_ norm_x_last = DT1_(0);
                //DT1_ norm_x = DT1_(1);

                DT1_ initial_norm(Norm<vnt_l_two, false, Tag_>::value(g));

                unsigned long i(0);
                while(i < iters)
                {
                    ++i;
                    cg_kernel(system_matrix, right_hand_side, g, x, u);
                    //norm_x = Norm<vnt_l_two, false, Tag_>::value(x);
                    //norm_x_last = Norm<vnt_l_two, false, Tag_>::value(x_last);
                    //copy<Tag_>(x, x_last);
                    DT1_ current_norm(Norm<vnt_l_two, false, Tag_>::value(Defect<Tag_>::value(right_hand_side, system_matrix, x)));
                    if(current_norm < initial_norm * 1e-08)
                    {
                        std::cout << "converged after " << i << " iterations. NORM: " << current_norm << std::endl;
                        break;
                    }
                }
                std::cout << "NO CONVERGENCE after " << i << " iterations!" << std::endl;
            }

            ///NEW ELL type:
            template <typename DT1_>
            static void value(SparseMatrixELL<DT1_> & system_matrix,
                    DenseVector<DT1_> & right_hand_side,
                    DenseVector<DT1_> & x,
                    unsigned long max_iters)
            {
                CONTEXT("When solving sparse ELL linear system with CG :");
#ifdef SOLVER_VERBOSE_L2
                std::cout << "Calling CG solver, preconditioning=NONE, datalayout=ELLPACK" << std::endl;
#endif
                DenseVector<DT1_> p(right_hand_side.size());
                DenseVector<DT1_> r(right_hand_side.size());
                DenseVector<DT1_> v(right_hand_side.size());
                cg_kernel(system_matrix, right_hand_side, x,  r, p, v, max_iters);
            }

            template <typename DT1_, typename DT2_>
            static void value(SparseMatrix<DT1_> & system_matrix,
                                               DenseVector<DT2_> & right_hand_side,
                                               DenseVector<DT2_> & x,
                                               double konv_rad)
            {
                CONTEXT("When solving sparse linear system with CG (with given convergence parameter):");
#ifdef SOLVER_VERBOSE_L2
                std::cout << "Calling CG solver, preconditioning=NONE, datalayout=SPARSE" << std::endl;
#endif
                DenseVector<DT1_> g = Product<Tag_>::value(system_matrix, x);
                Difference<Tag_>::value(g, right_hand_side);
                DenseVector<DT1_> g_c(g.copy());
                Scale<Tag_>::value(g_c, DT1_(-1.));
                DenseVector<DT1_> u(g_c.copy());
                DenseVector<DT1_> x_last(x.copy());

                DT1_ norm_x_last = DT1_(0);
                DT1_ norm_x = DT1_(1);

                while(norm_x - norm_x_last > konv_rad)
                {

                    cg_kernel(system_matrix, right_hand_side, g, x, u);
                    norm_x = Norm<vnt_l_two, false, Tag_>::value(x);
                    norm_x_last = Norm<vnt_l_two, false, Tag_>::value(x_last);
                    x_last = x.copy();
                }
            }

            ///Mixed precision implementations:
            template <typename DT1_, typename DT2_>
            static void value(BandedMatrix<DT1_> & system_matrix,
                                           DenseVector<DT2_> & right_hand_side,
                                           DenseVector<DT2_> & x,
                                           double konv_rad,
                                           int mixed_prec_iter_num)
            {
                CONTEXT("When solving banded linear system with CG (with given convergence parameter), MIXEDPREC:");
#ifdef SOLVER_VERBOSE_L2
                std::cout << "Calling CG solver, preconditioning=NONE, datalayout=BANDED, MIXEDPREC (k-version)" << std::endl;
#endif
                DenseVector<DT1_> g = Product<Tag_>::value(system_matrix, x);
                Difference<Tag_>::value(g, right_hand_side);
                DenseVector<DT1_> g_c(g.copy());
                Scale<Tag_>::value(g_c, DT1_(-1.));
                DenseVector<DT1_> u(g_c.copy());
                DenseVector<DT1_> x_last(x.copy());

                DT1_ norm_x_last = DT1_(0);
                DT1_ norm_x = DT1_(1);
                unsigned long iteration_count(0);

                BandedMatrix<double> bm_c(system_matrix.size());
                DenseVector<double> b_c(right_hand_side.size());

                while(norm_x - norm_x_last > konv_rad)
                {
                    if(iteration_count % mixed_prec_iter_num != 0)
                    {
                        cg_kernel(system_matrix, right_hand_side, g, x, u);
                        norm_x = Norm<vnt_l_two, false, Tag_>::value(x);
                        norm_x_last = Norm<vnt_l_two, false, Tag_>::value(x_last);
                        x_last = x.copy();
                    }
                    else
                    {
                        DenseVector<double> x_c(right_hand_side.size());
                        DenseVector<double> x_last_c(right_hand_side.size());
                        DenseVector<double> g_c(right_hand_side.size());
                        DenseVector<double> u_c(right_hand_side.size());

                        convert(bm_c, system_matrix);
                        convert(b_c, right_hand_side);
                        convert(x_c, x);
                        convert(x_last_c, x_last);
                        convert(g_c, g);
                        convert(u_c, u);
                        cg_kernel(bm_c, b_c, g_c, x_c, u_c);
                        norm_x = Norm<vnt_l_two, false, Tag_>::value(x_c);
                        norm_x_last = Norm<vnt_l_two, false, Tag_>::value(x_last_c);
                        x_last_c = x_c.copy();

                        convert(x, x_c);
                        convert(x_last, x_last_c);
                        convert(g, g_c);
                        convert(u, u_c);
                    }
                    ++iteration_count;
                }
            }
    };


//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /**
     * \brief Solution of linear system with PCG - Jacobi.
     *
     * \ingroup grpmatrixoperations
     * \ingroup grpvectoroperations
     */
    template <typename Tag_>
        struct ConjugateGradients<Tag_, methods::JAC>
        {
            private:
                template<typename DT1_, typename DT2_>
                    static inline void cg_kernel(DenseMatrix<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, DenseVector<DT1_> & former_gradient, DenseVector<DT1_> & former_result, DenseVector<DT1_> & utility, DenseVector<DT1_> & former_residual, DenseVector<DT1_> & diag_inverted)
                    {
                        CONTEXT("When calling CG kernel, preconditioner=JACOBI, datalayout=DENSE.");
#ifdef SOLVER_VERBOSE_L1
                        std::cout << "Calling single iteration of CG, preconditioner=JACOBI, datalayout=DENSE." << std::endl;
#endif

                        DT1_ alpha, beta, upper, lower;
                        upper = DotProduct<Tag_>::value(former_residual, former_gradient);
                        DenseVector<DT1_> temp = Product<Tag_>::value(system_matrix, utility);
                        lower = DotProduct<Tag_>::value(temp, utility);
                        if(fabs(lower) >= std::numeric_limits<DT1_>::epsilon())
                        {
                            alpha = upper/lower;
                        }
                        else
                        {
                            alpha = upper/ std::numeric_limits<DT1_>::epsilon();
                        }
                        ///Compute new result:
                        Scale<Tag_>::value(utility, alpha);
                        DenseVector<DT1_> temp2(utility.copy());
                        Sum<Tag_>::value(former_result, temp2);

                        ///Compute new residual:
                        Scale<Tag_>::value(temp, alpha);
                        DenseVector<DT1_> temp3(temp.copy());
                        Difference<Tag_>::value(former_residual, temp3);
                        ///Compute new gradient:
                        DenseVector<DT1_> r_c = former_residual.copy();
                        ElementProduct<Tag_>::value(r_c, diag_inverted);
                        former_gradient = r_c;

                        ///Compute new utility:
                        DT1_ upper_2;
                        upper_2 = DotProduct<Tag_>::value(former_gradient, former_residual);
                        lower = upper;
                        if(fabs(lower) >= std::numeric_limits<DT1_>::epsilon())
                        {
                            beta = upper_2/lower;
                        }
                        else
                        {
                            beta = upper_2/ std::numeric_limits<DT1_>::epsilon();
                        }

                        Scale<Tag_>::value(utility, beta);
                        Sum<Tag_>::value(utility, former_gradient);
            }

            template<typename DT1_, typename DT2_>
            static inline void cg_kernel(BandedMatrix<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, DenseVector<DT1_> & former_gradient, DenseVector<DT1_> & former_result, DenseVector<DT1_> & utility, DenseVector<DT1_> & former_residual, DenseVector<DT1_> & diag_inverted)
            {
                CONTEXT("When calling CG kernel, preconditioner=JACOBI, datalayout=BANDED.");
#ifdef SOLVER_VERBOSE_L1
                std::cout << "Calling single iteration of CG, preconditioner=JACOBI, datalayout=BANDED." << std::endl;
#endif

                DT1_ alpha, beta, upper, lower;
                upper = DotProduct<Tag_>::value(former_residual, former_gradient);
                DenseVector<DT1_> temp = Product<Tag_>::value(system_matrix, utility);
                lower = DotProduct<Tag_>::value(temp, utility);
                if(fabs(lower) >= std::numeric_limits<DT1_>::epsilon())
                {
                    alpha = upper/lower;
                }
                else
                {
                    alpha = upper/ std::numeric_limits<DT1_>::epsilon();
                }
                ///Compute new result:
                Scale<Tag_>::value(utility, alpha);
                DenseVector<DT1_> temp2(utility.copy());
                Sum<Tag_>::value(former_result, temp2);

                ///Compute new residual:
                Scale<Tag_>::value(temp, alpha);
                DenseVector<DT1_> temp3(temp.copy());
                Difference<Tag_>::value(former_residual, temp3);
                ///Compute new gradient:
                DenseVector<DT1_> r_c = former_residual.copy();
                ElementProduct<Tag_>::value(r_c, diag_inverted);
                former_gradient = r_c;

                ///Compute new utility:
                DT1_ upper_2;
                upper_2 = DotProduct<Tag_>::value(former_gradient, former_residual);
                lower = upper;
                if(fabs(lower) >= std::numeric_limits<DT1_>::epsilon())
                {
                    beta = upper_2/lower;
                }
                else
                {
                    beta = upper_2/ std::numeric_limits<DT1_>::epsilon();
                }

                Scale<Tag_>::value(utility, beta);
                Sum<Tag_>::value(utility, former_gradient);
            }

            template<typename DT1_, typename DT2_>
            static inline void cg_kernel(SparseMatrix<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, DenseVector<DT1_> & former_gradient, DenseVector<DT1_> & former_result, DenseVector<DT1_> & utility, DenseVector<DT1_> & former_residual, DenseVector<DT1_> & diag_inverted)
            {
                CONTEXT("When calling CG kernel, preconditioner=JACOBI, datalayout=SPARSE.");
#ifdef SOLVER_VERBOSE_L1
                std::cout << "Calling single iteration of CG, preconditioner=JACOBI, datalayout=SPARSE." << std::endl;
#endif

                DT1_ alpha, beta, upper, lower;
                upper = DotProduct<Tag_>::value(former_residual, former_gradient);
                DenseVector<DT1_> temp = Product<Tag_>::value(system_matrix, utility);
                lower = DotProduct<Tag_>::value(temp, utility);
                if(fabs(lower) >= std::numeric_limits<DT1_>::epsilon())
                {
                    alpha = upper/lower;
                }
                else
                {
                    alpha = upper/ std::numeric_limits<DT1_>::epsilon();
                }
                ///Compute new result:
                Scale<Tag_>::value(utility, alpha);
                DenseVector<DT1_> temp2(utility.copy());
                Sum<Tag_>::value(former_result, temp2);

                ///Compute new residual:
                Scale<Tag_>::value(temp, alpha);
                DenseVector<DT1_> temp3(temp.copy());
                Difference<Tag_>::value(former_residual, temp3);
                ///Compute new gradient:
                DenseVector<DT1_> r_c = former_residual.copy();
                ElementProduct<Tag_>::value(r_c, diag_inverted);
                former_gradient = r_c;

                ///Compute new utility:
                DT1_ upper_2;
                upper_2 = DotProduct<Tag_>::value(former_gradient, former_residual);
                lower = upper;
                if(fabs(lower) >= std::numeric_limits<DT1_>::epsilon())
                {
                    beta = upper_2/lower;
                }
                else
                {
                    beta = upper_2/ std::numeric_limits<DT1_>::epsilon();
                }

                Scale<Tag_>::value(utility, beta);
                Sum<Tag_>::value(utility, former_gradient);
            }

            ///SMELL type:
            template<typename DT_>
            static inline void cg_kernel(SparseMatrixELL<DT_> & A,
                                         DenseVector<DT_> & r,
                                         DenseVector<DT_> & z,
                                         DenseVector<DT_> & d,
                                         DenseVector<DT_> & x,
                                         DenseVector<DT_> & dd_inverted,
                                         DenseVector<DT_> & temp_0)
            {
                CONTEXT("When calling CG kernel, preconditioner=JACOBI, datalayout=ELLPACK.");
#ifdef SOLVER_VERBOSE_L1
                std::cout << "Calling single iteration of CG, preconditioner=JACOBI, datalayout=ELL." << std::endl;
#endif

                DT_ alpha, beta;
                beta = DotProduct<Tag_>::value(r, z);
                Product<Tag_>::value(temp_0, A, d);
                alpha = DotProduct<Tag_>::value(d, temp_0);
                alpha = std::fabs(alpha) >= std::numeric_limits<DT_>::epsilon() ? beta / alpha : beta / std::numeric_limits<DT_>::epsilon();
                ScaledSum<Tag_>::value(x, d, alpha);
                ScaledSum<Tag_>::value(r, temp_0, -alpha);

                ElementProduct<Tag_>::value(z, dd_inverted, r);
                beta = std::fabs(beta) >= std::numeric_limits<DT_>::epsilon() ? DotProduct<Tag_>::value(z, r) / beta : DotProduct<Tag_>::value(z, r) / std::numeric_limits<DT_>::epsilon();

                copy<Tag_>(d, temp_0);
                ScaledSum<Tag_>::value(d, z, temp_0, beta);

            }
            //end SMELL type

        public:
            /**
            * \brief Returns solution of linear system given by a DenseMatrix and a Vector.
            *
            * \param system_matrix The system matrix.
            * \param right_hand_side The right hand side of the system.
            * \param x The solution (and initial guess).
            * \param konv_rad The parameter for convergence control.
            *
            */
            template <typename DT1_, typename DT2_>
            static void value(DenseMatrix<DT1_> & system_matrix,
                                           DenseVector<DT2_> & right_hand_side,
                                           DenseVector<DT2_> & x,
                                           double konv_rad)
            {
                CONTEXT("When solving dense linear system with PCG-Jacobi (with given convergence parameter):");
#ifdef SOLVER_VERBOSE_L2
                std::cout << "Calling CG solver, preconditioning=JACOBI, datalayout=DENSE" << std::endl;
#endif

                DenseVector<DT1_> r = Product<Tag_>::value(system_matrix, x);
                Difference<Tag_>::value(r, right_hand_side);
                Scale<Tag_>::value(r, DT1_(-1.));

                ///PCG - Jacobi part:
                //TODO: remove diag, we dont need it really
                DenseVector<DT1_> diag(right_hand_side.size(), DT1_(0));
                DenseVector<DT1_> diag_inverted(right_hand_side.size(), DT1_(0));
                ///Create diagonal, invert.
                for(unsigned long i =0; i < diag.size(); ++i)
                {

                    diag[i] = system_matrix[i][i];
                    if(fabs(diag[i]) >= std::numeric_limits<DT1_>::epsilon())
                    {
                        diag_inverted[i] = DT1_(1) / diag[i];
                    }
                    else
                    {
                        diag_inverted[i] = DT1_(1) / std::numeric_limits<DT1_>::epsilon();
                    }
                }
                ///Now, we have computed "C¯1" alias diag_inverted.
                ///g has to be multiplied elementwisely by diag_inverted:
                DenseVector<DT1_> temp(r.copy());
                ElementProduct<Tag_>::value(temp,diag_inverted);
                DenseVector<DT1_> g(temp);
                DenseVector<DT1_> u(g.copy());

                ///End of the PCG part, we will give diag_inverted to the kernel, in order not to be forced to compute it again.

                DenseVector<DT1_> x_last(x.copy());
                DT1_ norm_x_last = DT1_(0);
                DT1_ norm_x = DT1_(1);

                while(fabs(norm_x - norm_x_last) > konv_rad)
                {

                    cg_kernel(system_matrix, right_hand_side, g, x, u, r, diag_inverted);
                    norm_x = Norm<vnt_l_two, false, Tag_>::value(x);
                    norm_x_last = Norm<vnt_l_two, false, Tag_>::value(x_last);
                    x_last = x.copy();
                }
            }

            /**
            * \brief Returns solution of linear system given by a BandedMatrix and a Vector.
            *
            * \param system_matrix The system matrix.
            * \param right_hand_side The right hand side of the system.
            * \param x The solution (and initial guess).
            * \param konv_rad The parameter for convergence control.
            *
            */
            template <typename DT1_, typename DT2_>
            static void value(BandedMatrix<DT1_> & system_matrix,
                              DenseVector<DT2_> & right_hand_side,
                              DenseVector<DT2_> & x,
                              double konv_rad)
            {
                CONTEXT("When solving banded linear system with PCG-Jacobi (with given convergence parameter):");
#ifdef SOLVER_VERBOSE_L2
                std::cout << "Calling CG solver, preconditioning=JACOBI, datalayout=BANDED" << std::endl;
#endif

                DenseVector<DT1_> r = Product<Tag_>::value(system_matrix, x);
                Difference<Tag_>::value(r, right_hand_side);
                Scale<Tag_>::value(r, DT1_(-1.));

                ///PCG - Jacobi part:
                //TODO: remove diag, we dont need it really
                DenseVector<DT1_> diag(system_matrix.band(0));
                DenseVector<DT1_> diag_inverted(right_hand_side.size(), DT1_(0));
                ///Create diagonal, invert.
                for(unsigned long i =0; i < diag.size(); ++i)
                {

                    if(fabs(diag[i]) >= std::numeric_limits<DT1_>::epsilon())
                    {
                        diag_inverted[i] = DT1_(1) / diag[i];
                    }
                    else
                    {
                        diag_inverted[i] = DT1_(1) / std::numeric_limits<DT1_>::epsilon();
                    }
                }
                ///Now, we have computed "C¯1" alias diag_inverted.
                ///g has to be multiplied elementwisely by diag_inverted:
                DenseVector<DT1_> temp(r.copy());
                ElementProduct<Tag_>::value(temp,diag_inverted);
                DenseVector<DT1_> g(temp);
                DenseVector<DT1_> u(g.copy());

                ///End of the PCG part, we will give diag_inverted to the kernel, in order not to be forced to compute it again.

                DenseVector<DT1_> x_last(x.copy());
                DT1_ norm_x_last = DT1_(0);
                DT1_ norm_x = DT1_(1);

                while(fabs(norm_x - norm_x_last) > konv_rad)
                {

                    cg_kernel(system_matrix, right_hand_side, g, x, u, r, diag_inverted);
                    norm_x = Norm<vnt_l_two, false, Tag_>::value(x);
                    norm_x_last = Norm<vnt_l_two, false, Tag_>::value(x_last);
                    x_last = x.copy();
                }
            }

            /**
            * \brief Returns solution of linear system given by a SparseMatrix and a Vector.
            *
            * \param system_matrix The system matrix.
            * \param right_hand_side The right hand side of the system.
            * \param x The solution (and initial guess).
            * \param konv_rad The parameter for convergence control.
            *
            */
            template <typename DT1_, typename DT2_>
            static void value(SparseMatrix<DT1_> & system_matrix,
                         DenseVector<DT2_> & right_hand_side,
                         DenseVector<DT2_> & x,
                         double konv_rad)
            {
                CONTEXT("When solving sparse linear system with PCG-Jacobi (with given convergence parameter):");
#ifdef SOLVER_VERBOSE_L2
                std::cout << "Calling CG solver, preconditioning=JACOBI, datalayout=SPARSE" << std::endl;
#endif

                DenseVector<DT1_> r = Product<Tag_>::value(system_matrix, x);
                Difference<Tag_>::value(r, right_hand_side);
                Scale<Tag_>::value(r, DT1_(-1.));

                ///PCG - Jacobi part:
                //TODO: remove diag, we dont need it really
                DenseVector<DT1_> diag(right_hand_side.size(), DT1_(0));
                DenseVector<DT1_> diag_inverted(right_hand_side.size(), DT1_(0));
                ///Create diagonal, invert.
                for(unsigned long i =0; i < diag.size(); ++i)
                {
                    diag[i] = system_matrix[i][i];
                    if(fabs(diag[i]) >= std::numeric_limits<DT1_>::epsilon())
                    {
                        diag_inverted[i] = DT1_(1) / diag[i];
                    }
                    else
                    {
                        diag_inverted[i] = DT1_(1) / std::numeric_limits<DT1_>::epsilon();
                    }
                }
                ///Now, we have computed "C¯1" alias diag_inverted.
                ///g has to be multiplied elementwisely by diag_inverted:
                DenseVector<DT1_> temp(r.copy());
                ElementProduct<Tag_>::value(temp,diag_inverted);
                DenseVector<DT1_> g(temp);
                DenseVector<DT1_> u(g.copy());

                ///End of the PCG part, we will give diag_inverted to the kernel, in order not to be forced to compute it again.

                DenseVector<DT1_> x_last(x.copy());
                DT1_ norm_x_last = DT1_(0);
                DT1_ norm_x = DT1_(1);

                while(fabs(norm_x - norm_x_last) > konv_rad)
                {

                    cg_kernel(system_matrix, right_hand_side, g, x, u, r, diag_inverted);
                    norm_x = Norm<vnt_l_two, false, Tag_>::value(x);
                    norm_x_last = Norm<vnt_l_two, false, Tag_>::value(x_last);
                    x_last = x.copy();
                }
            }


            template <typename DT_>
            static void value(SparseMatrixELL<DT_> & system_matrix,
                                           DenseVector<DT_> & right_hand_side,
                                           DenseVector<DT_> & x,
                                           DenseVector<DT_> & dd_inverted,
                                           unsigned long max_iters)
            {
                CONTEXT("When solving sparse ELL linear system with PCG-Jacobi: ");
#ifdef SOLVER_VERBOSE_L2
                std::cout << "Calling CG solver, preconditioning=JACOBI, datalayout=ELL" << std::endl;
#endif

                DenseVector<DT_> t_0(right_hand_side.size());

                DenseVector<DT_> r(right_hand_side.size());
                Defect<Tag_>::value(r, right_hand_side, system_matrix, x);
                DT_ initial_defect_norm(Norm<vnt_l_two, true, Tag_>::value(r));
                std::cout << "Initial defect NORM: " << initial_defect_norm << std::endl;

                DenseVector<DT_> z(right_hand_side.size());
                ElementProduct<Tag_>::value(z, dd_inverted, r);
                DenseVector<DT_> d(z.size());
                copy<Tag_>(z , d);

                DenseVector<DT_> t_1(right_hand_side.size());
                for(unsigned long i(0) ; i < max_iters ; ++i)
                {
                    cg_kernel(system_matrix, r, z, d, x, dd_inverted, t_0);
                    DT_ current_defect_norm(Norm<vnt_l_two, true, Tag_>::value(r));

                    if(current_defect_norm < initial_defect_norm * 1e-8)
                    {
                        std::cout << "Converged after " << i + 1 << " iterations: NORM: " << current_defect_norm << std::endl;
                        break;
                    }
                    if(i == max_iters - 1)
                    {
                        std::cout << "NO convergence after " << i + 1 << " iterations: NORM: " << current_defect_norm << std::endl;
                    }
                }

            }

        };
}

#endif
