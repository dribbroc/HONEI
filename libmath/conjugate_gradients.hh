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

#include <libutil/tags.hh>
#include <libla/dense_matrix.hh>
#include <libla/dense_vector.hh>
#include <libla/product.hh>
#include <libla/sum.hh>
#include <libla/difference.hh>
#include <libla/dot_product.hh>
#include <libla/scale.hh>
#include <libla/norm.hh>
#include <iostream>
#include <libmath/methods.hh>
#include <libla/element_product.hh>
/**
 * \file
 *
 * Templatized definition and implementation of a Conjugate Gradients Solver.
 *
 * \ingroup grpoperations
 */

using namespace std;
using namespace methods;
namespace honei
{

    template<typename Tag_, typename PreCon_>
    struct ConjugateGradients
    {
    };
    /**
     * \brief Solution of LES with CG.
     *
     * The referenced containers are invariant under this operation.
     * In every case, a new object is created and returned.
     *
     * \ingroup grpmatrixoperations
     * \ingroup grpvectoroperations
     */

    template <>
    struct ConjugateGradients<tags::CPU, methods::NONE>
    {
        private:
            template<typename DT1_, typename DT2_>
            static inline void cg_kernel(DenseMatrix<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, DenseVector<DT1_> & former_gradient, DenseVector<DT1_> & former_result, DenseVector<DT1_> & utility)
            {
                ///Compute x_i+1:
                DT1_ upper = DotProduct<>::value(former_gradient, former_gradient);
                DenseVector<DT1_> energy = Product<>::value(system_matrix, utility);
                DT1_ lower = DotProduct<>::value(energy, utility);
                DenseVector<DT1_> u_c(utility.copy());
                if(fabs(lower) >= std::numeric_limits<DT1_>::epsilon())
                {
                    energy = Scale<>::value(DT1_(upper/lower), u_c);
                }
                else
                {
                    energy = Scale<>::value(DT1_(upper/std::numeric_limits<DT1_>::epsilon()), u_c);

                }
                Sum<>::value(energy, former_result);
                ///Compute new gradient
                DenseVector<DT1_> new_gradient = Product<>::value(system_matrix, energy);
                Difference<>::value(new_gradient, right_hand_side);

                ///Compute new utility
                upper = DotProduct<>::value(new_gradient, new_gradient);
                lower = DotProduct<>::value(former_gradient, former_gradient);
                if(fabs(lower) >= std::numeric_limits<DT1_>::epsilon())
                {
                    Scale<>::value(DT1_(upper/lower), utility);
                }
                else
                {
                    Scale<>::value(DT1_(upper/std::numeric_limits<DT1_>::epsilon()), utility);
                }
                Difference<>::value(utility, new_gradient);

                ///Finishing:
                former_gradient = new_gradient.copy();
                former_result = energy.copy();

            }

            template<typename DT1_, typename DT2_>
            static inline void cg_kernel(BandedMatrix<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, DenseVector<DT1_> & former_gradient, DenseVector<DT1_> & former_result, DenseVector<DT1_> & utility)
            {
                ///Compute x_i+1: (in energy)
                DT1_ upper = DotProduct<>::value(former_gradient, former_gradient);
                DenseVector<DT1_> energy = Product<>::value(system_matrix, utility);
                DT1_ lower = DotProduct<>::value(energy, utility);
                DenseVector<DT1_> u_c(utility.copy());
                if(fabs(lower) >= std::numeric_limits<DT1_>::epsilon())
                {
                    energy = Scale<>::value(DT1_(upper/lower), u_c);
                }
                else
                {
                    energy = Scale<>::value(DT1_(upper/std::numeric_limits<DT1_>::epsilon()), u_c);
                }
                Sum<>::value(energy, former_result);
                ///Compute new gradient
                DenseVector<DT1_> new_gradient = Product<>::value(system_matrix, energy);
                Difference<>::value(new_gradient, right_hand_side);

                ///Compute new utility
                upper = DotProduct<>::value(new_gradient, new_gradient);
                lower = DotProduct<>::value(former_gradient, former_gradient);
                if(fabs(lower) >= std::numeric_limits<DT1_>::epsilon())
                {
                    Scale<>::value(DT1_(upper/lower), utility);
                }
                else
                {
                    Scale<>::value(DT1_(upper/std::numeric_limits<DT1_>::epsilon()), utility);
                }
                Difference<>::value(utility, new_gradient);

                ///Finishing:
                former_gradient = new_gradient.copy();
                former_result = energy.copy();

            }

            template<typename DT1_, typename DT2_>
            static inline void cg_kernel(SparseMatrix<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, DenseVector<DT1_> & former_gradient, DenseVector<DT1_> & former_result, DenseVector<DT1_> & utility)
            {
                ///Compute x_i+1: (in energy)
                DT1_ upper = DotProduct<>::value(former_gradient, former_gradient);
                SparseVector<DT1_> energy = Product<>::value(system_matrix, utility);
                DT1_ lower = DotProduct<>::value(energy, utility);
                DenseVector<DT1_> u_c(utility.copy());
                if(fabs(lower) >= std::numeric_limits<DT1_>::epsilon())
                {
                    energy = Scale<>::value(DT1_(upper/lower), u_c);
                }
                else
                {
                    energy = Scale<>::value(DT1_(upper/std::numeric_limits<DT1_>::epsilon()), u_c);
                }
                Sum<>::value(energy, former_result);
                ///Compute new gradient
                SparseVector<DT1_> new_gradient = Product<>::value(system_matrix, energy);
                Difference<>::value(new_gradient, right_hand_side);

                ///Compute new utility
                upper = DotProduct<>::value(new_gradient, new_gradient);
                lower = DotProduct<>::value(former_gradient, former_gradient);
                if(fabs(lower) >= std::numeric_limits<DT1_>::epsilon())
                {
                    Scale<>::value(DT1_(upper/lower), utility);
                }
                else
                {
                    Scale<>::value(DT1_(upper/std::numeric_limits<DT1_>::epsilon()), utility);
                }
                Difference<>::value(utility, new_gradient);

                ///Finishing:
                former_gradient = new_gradient.copy();
                former_result = energy.copy();

            }




        public:
            /**
            * \brief Returns solution of LES given by a DenseMatrix and a Vector.
            *
            * \param system_matrix The system matrix.
            * \param right_hand_side The right hand side of the system.
            * \param iter_number The fixed number of iterations.
            *
           */

            /// \{


            template <typename DT1_, typename DT2_>
            static DenseVector<DT1_> value(DenseMatrix<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side,long iter_number)
            {
                CONTEXT("When solving dense linear system with CG (fixed # iterations):");


                DenseVector<DT1_> x(right_hand_side.size(), DT1_(0));
                DenseVector<DT1_> g = Product<>::value(system_matrix, x);
                Difference<>::value(g, right_hand_side);
                DenseVector<DT1_> g_c(g.copy());
                DenseVector<DT1_> u = Scale<>::value(DT1_(-1.), g_c);
                for(unsigned long i = 0; i<iter_number; ++i)
                {
                    cg_kernel(system_matrix, right_hand_side, g, x, u);
                }
                return x;
            }
            /**
            * \brief Returns solution of LES given by a DenseMatrix and a Vector.
            *
            * \param system_matrix The system matrix.
            * \param right_hand_side The right hand side of the system.
            * \param konv_rad The parameter for convergence control.
            *
            */

            /// \{


            template <typename DT1_, typename DT2_>
            static DenseVector<DT1_> value(DenseMatrix<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, double konv_rad)
            {
                CONTEXT("When solving dense linear system with CG (with given convergence parameter):");


                DenseVector<DT1_> x(right_hand_side.size(), DT1_(0));
                DenseVector<DT1_> g = Product<>::value(system_matrix, x);
                Difference<>::value(g, right_hand_side);
                DenseVector<DT1_> g_c(g.copy());
                DenseVector<DT1_> u = Scale<>::value(DT1_(-1.), g_c);
                DenseVector<DT1_> x_last(x.copy());

                DT1_ norm_x_last = DT1_(0);
                DT1_ norm_x = DT1_(1);

                while(norm_x - norm_x_last > konv_rad)
                {

                    cg_kernel(system_matrix, right_hand_side, g, x, u);
                    norm_x = Norm<vnt_l_two, false, tags::CPU>::value(x);
                    norm_x_last = Norm<vnt_l_two, false, tags::CPU>::value(x_last);
                    x_last = x.copy();
                }
                return x;


            }
            /**
            * \brief Returns solution of LES given by a BandedMatrix and a Vector.
            *
            * \param system_matrix The system matrix.
            * \param right_hand_side The right hand side of the system.
            * \param iter_number The fixed number of iterations.
            *
            */

            /// \{



            template <typename DT1_, typename DT2_>
            static DenseVector<DT1_> value(BandedMatrix<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side,long iter_number)
            {
                CONTEXT("When solving banded linear system with CG (with fixed # of iterations):");



                DenseVector<DT1_> x(right_hand_side.size(), DT1_(0));
                DenseVector<DT1_> g = Product<>::value(system_matrix, x);
                Difference<>::value(g, right_hand_side);
                DenseVector<DT1_> g_c(g.copy());
                DenseVector<DT1_> u = Scale<>::value(DT1_(-1.), g_c);
                for(unsigned long i = 0; i<iter_number; ++i)
                {
                    cg_kernel(system_matrix, right_hand_side, g, x, u);
                }
                return x;
            }
            /**
            * \brief Returns solution of LES given by a BandedMatrix and a Vector.
            *
            * \param system_matrix The system matrix.
            * \param right_hand_side The right hand side of the system.
            * \param konv_rad The parameter for convergence control.
            *
            */

            /// \{


            template <typename DT1_, typename DT2_>
            static DenseVector<DT1_> value(BandedMatrix<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, double konv_rad)
            {
                CONTEXT("When solving banded linear system with CG (with given convergence parameter):");


                DenseVector<DT1_> x(right_hand_side.size(), DT1_(0));
                DenseVector<DT1_> g = Product<>::value(system_matrix, x);
                Difference<>::value(g, right_hand_side);
                DenseVector<DT1_> g_c(g.copy());
                DenseVector<DT1_> u = Scale<>::value(DT1_(-1.), g_c);
                DenseVector<DT1_> x_last(x.copy());

                DT1_ norm_x_last = DT1_(0);
                DT1_ norm_x = DT1_(1);

                while(norm_x - norm_x_last > konv_rad)
                {

                    cg_kernel(system_matrix, right_hand_side, g, x, u);
                    norm_x = Norm<vnt_l_two, false, tags::CPU>::value(x);
                    norm_x_last = Norm<vnt_l_two, false, tags::CPU>::value(x_last);
                    x_last = x.copy();
                }
                return x;


            }
   };

    /**
     * \brief Solution of LES with PCG - Jacobi.
     *
     * The referenced containers are invariant under this operation.
     * In every case, a new object is created and returned.
     *
     * \ingroup grpmatrixoperations
     * \ingroup grpvectoroperations
     */

    template <>
    struct ConjugateGradients<tags::CPU, methods::JAC>
    {
        private:
            template<typename DT1_, typename DT2_>
            static inline void cg_kernel(DenseMatrix<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, DenseVector<DT1_> & former_gradient, DenseVector<DT1_> & former_result, DenseVector<DT1_> & utility, DenseVector<DT1_> & former_residual, DenseVector<DT1_> & diag_inverted)
            {

                DT1_ alpha, beta, upper, lower;
                upper = DotProduct<>::value(former_residual, former_gradient);
                DenseVector<DT1_> temp = Product<>::value(system_matrix, utility);
                lower = DotProduct<>::value(temp, utility);
                if(fabs(lower) >= std::numeric_limits<DT1_>::epsilon())
                {
                    alpha = upper/lower;
                }
                else
                {
                    alpha = upper/ std::numeric_limits<DT1_>::epsilon();
                }
                ///Compute new result:
                DenseVector<DT1_> temp2 = Scale<>::value(alpha, utility);
                Sum<>::value(former_result, temp2);

                ///Compute new residual:
                DenseVector<DT1_> temp3 = Scale<>::value(alpha, temp);
                Difference<>::value(former_residual, temp3);
                ///Compute new gradient:
                DenseVector<DT1_> r_c = former_residual.copy();
                former_gradient = ElementProduct<>::value(r_c, diag_inverted);

                ///Compute new utility:
                DT1_ upper_2;
                upper_2 = DotProduct<>::value(former_gradient, former_residual);
                lower = upper;
                if(fabs(lower) >= std::numeric_limits<DT1_>::epsilon())
                {
                    beta = upper_2/lower;
                }
                else
                {
                    beta = upper_2/ std::numeric_limits<DT1_>::epsilon();
                }

                Scale<>::value(beta, utility);
                Sum<>::value(utility, former_gradient);
            }

            template<typename DT1_, typename DT2_>
            static inline void cg_kernel(BandedMatrix<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, DenseVector<DT1_> & former_gradient, DenseVector<DT1_> & former_result, DenseVector<DT1_> & utility, DenseVector<DT1_> & former_residual, DenseVector<DT1_> & diag_inverted)
            {

                DT1_ alpha, beta, upper, lower;
                upper = DotProduct<>::value(former_residual, former_gradient);
                DenseVector<DT1_> temp = Product<>::value(system_matrix, utility);
                lower = DotProduct<>::value(temp, utility);
                if(fabs(lower) >= std::numeric_limits<DT1_>::epsilon())
                {
                    alpha = upper/lower;
                }
                else
                {
                    alpha = upper/ std::numeric_limits<DT1_>::epsilon();
                }
                ///Compute new result:
                DenseVector<DT1_> temp2 = Scale<>::value(alpha, utility);
                Sum<>::value(former_result, temp2);

                ///Compute new residual:
                DenseVector<DT1_> temp3 = Scale<>::value(alpha, temp);
                Difference<>::value(former_residual, temp3);
                ///Compute new gradient:
                DenseVector<DT1_> r_c = former_residual.copy();
                former_gradient = ElementProduct<>::value(r_c, diag_inverted);

                ///Compute new utility:
                DT1_ upper_2;
                upper_2 = DotProduct<>::value(former_gradient, former_residual);
                lower = upper;
                if(fabs(lower) >= std::numeric_limits<DT1_>::epsilon())
                {
                    beta = upper_2/lower;
                }
                else
                {
                    beta = upper_2/ std::numeric_limits<DT1_>::epsilon();
                }

                Scale<>::value(beta, utility);
                Sum<>::value(utility, former_gradient);
            }

        public:
            /**
            * \brief Returns solution of LES given by a DenseMatrix and a Vector.
            *
            * \param system_matrix The system matrix.
            * \param right_hand_side The right hand side of the system.
            * \param konv_rad The parameter for convergence control.
            *
            */

            /// \{

            template <typename DT1_, typename DT2_>
            static DenseVector<DT1_> value(DenseMatrix<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, double konv_rad)
            {
                CONTEXT("When solving dense linear system with PCG-Jacobi (with given convergence parameter):");


                DenseVector<DT1_> x(right_hand_side.size(), DT1_(0));
                DenseVector<DT1_> r = Product<>::value(system_matrix, x);
                Difference<>::value(r, right_hand_side);
                Scale<>::value(DT1_(-1.), r);

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
                DenseVector<DT1_> g(ElementProduct<tags::CPU>::value(temp,diag_inverted));
                DenseVector<DT1_> u(g.copy());

                ///End of the PCG part, we will give diag_inverted to the kernel, in order not to be forced to compute it again.

                DenseVector<DT1_> x_last(x.copy());
                DT1_ norm_x_last = DT1_(0);
                DT1_ norm_x = DT1_(1);

                while(fabs(norm_x - norm_x_last) > konv_rad)
                {

                    cg_kernel(system_matrix, right_hand_side, g, x, u, r, diag_inverted);
                    norm_x = Norm<vnt_l_two, false, tags::CPU>::value(x);
                    norm_x_last = Norm<vnt_l_two, false, tags::CPU>::value(x_last);
                    x_last = x.copy();
                }
                return x;


            }

            /**
            * \brief Returns solution of LES given by a BandedMatrix and a Vector.
            *
            * \param system_matrix The system matrix.
            * \param right_hand_side The right hand side of the system.
            * \param konv_rad The parameter for convergence control.
            *
            */

            /// \{

            template <typename DT1_, typename DT2_>
            static DenseVector<DT1_> value(BandedMatrix<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, double konv_rad)
            {
                CONTEXT("When solving banded linear system with PCG-Jacobi (with given convergence parameter):");


                DenseVector<DT1_> x(right_hand_side.size(), DT1_(0));
                DenseVector<DT1_> r = Product<>::value(system_matrix, x);
                Difference<>::value(r, right_hand_side);
                Scale<>::value(DT1_(-1.), r);

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
                DenseVector<DT1_> g(ElementProduct<tags::CPU>::value(temp,diag_inverted));
                DenseVector<DT1_> u(g.copy());

                ///End of the PCG part, we will give diag_inverted to the kernel, in order not to be forced to compute it again.

                DenseVector<DT1_> x_last(x.copy());
                DT1_ norm_x_last = DT1_(0);
                DT1_ norm_x = DT1_(1);

                while(fabs(norm_x - norm_x_last) > konv_rad)
                {

                    cg_kernel(system_matrix, right_hand_side, g, x, u, r, diag_inverted);
                    norm_x = Norm<vnt_l_two, false, tags::CPU>::value(x);
                    norm_x_last = Norm<vnt_l_two, false, tags::CPU>::value(x_last);
                    x_last = x.copy();
                }
                return x;


            }

    };
}

#endif
