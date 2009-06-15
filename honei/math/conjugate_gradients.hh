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
            static inline void cg_kernel(SparseMatrixELL<DT1_> & system_matrix, DenseVector<DT1_> & right_hand_side, DenseVector<DT1_> & diag_inverted, DenseVector<DT1_> & former_gradient, DenseVector<DT1_> & former_result, DenseVector<DT1_> & utility)
            {
                ///Compute x_i+1: (in energy)
                DT1_ upper = DotProduct<Tag_>::value(former_gradient, former_gradient);
                DenseVector<DT1_> energy(former_result.size());
                Product<Tag_>::value(energy, system_matrix, utility);
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
                DenseVector<DT1_> new_gradient(energy.size());
                Product<Tag_>::value(new_gradient, system_matrix, energy);
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
            * \brief Returns solution of LES given by a DenseMatrix and a Vector.
            *
            * \param system_matrix The system matrix.
            * \param right_hand_side The right hand side of the system.
            * \param iter_number The fixed number of iterations.
            *
           */

            /// \{


            template <typename DT1_, typename DT2_>
            static DenseVector<DT1_> value(DenseMatrix<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, unsigned long iter_number)
            {
                CONTEXT("When solving dense linear system with CG (fixed # iterations):");


                DenseVector<DT1_> x(right_hand_side.size(), DT1_(0));
                DenseVector<DT1_> g = Product<Tag_>::value(system_matrix, x);
                Difference<Tag_>::value(g, right_hand_side);
                DenseVector<DT1_> g_c(g.copy());
                Scale<Tag_>::value(g_c, DT1_(-1.));
                DenseVector<DT1_> u(g_c.copy());
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
            static DenseVector<DT1_> value(BandedMatrix<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, unsigned long iter_number)
            {
                CONTEXT("When solving banded linear system with CG (with fixed # of iterations):");



                DenseVector<DT1_> x(right_hand_side.size(), DT1_(0));
                DenseVector<DT1_> g = Product<Tag_>::value(system_matrix, x);
                Difference<Tag_>::value(g, right_hand_side);
                DenseVector<DT1_> g_c(g.copy());
                Scale<Tag_>::value(g_c, DT1_(-1.));
                DenseVector<DT1_> u(g_c.copy());
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
            static DenseVector<DT1_> value(BandedMatrixQ1<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, double konv_rad)
            {
                CONTEXT("When solving banded Q1 linear system with CG (with given convergence parameter):");

                DenseVector<DT1_> x(right_hand_side.size());
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
                return x;
            }


            template <typename DT1_, typename DT2_>
            static DenseVector<DT1_> value(BandedMatrixQ1<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, unsigned long iters)
            {
                CONTEXT("When solving banded Q1 linear system with CG (with given convergence parameter):");

                DenseVector<DT1_> x(right_hand_side.size());
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

                unsigned long i(0);
                while(i < iters)
                {
                    ++i;
                    cg_kernel(system_matrix, right_hand_side, g, x, u);
                    norm_x = Norm<vnt_l_two, false, Tag_>::value(x);
                    norm_x_last = Norm<vnt_l_two, false, Tag_>::value(x_last);
                    copy<Tag_>(x, x_last);
                }
                return x;
            }

            ///NEW ELL type:
            template <typename DT1_>
            static DenseVector<DT1_> value(SparseMatrixELL<DT1_> & system_matrix, DenseVector<DT1_> & right_hand_side, DenseVector<DT1_> & x, DenseVector<DT1_> & diag_inverted, unsigned long iters)
            {
                CONTEXT("When solving banded Q1 linear system with CG (with given convergence parameter):");

                //DenseVector<DT1_> x(right_hand_side.size());
                //fill<Tag_>(x, DT1_(0));
                DenseVector<DT1_> g(right_hand_side.size());
                Product<Tag_>::value(g, system_matrix, x);
                Difference<Tag_>::value(g, right_hand_side);
                DenseVector<DT1_> g_c(g.size());
                copy<Tag_>(g, g_c);
                Scale<Tag_>::value(g_c, DT1_(-1.));
                DenseVector<DT1_> u(g_c.size());
                copy<Tag_>(g_c, u);
                DenseVector<DT1_> x_last(x.size());
                copy<Tag_>(x, x_last);

                unsigned long i(0);
                while(i < iters)
                {
                    ++i;
                    cg_kernel(system_matrix, right_hand_side, diag_inverted, g, x, u);
                    copy<Tag_>(x, x_last);
                }
                return x;
            }

            ///Mixed precision implementations:
            template <typename DT1_, typename DT2_>
                static DenseVector<DT1_> value(BandedMatrix<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, double konv_rad, int mixed_prec_iter_num)
            {
                CONTEXT("When solving banded linear system with CG (with given convergence parameter), MIXEDPREC:");


                DenseVector<DT1_> x(right_hand_side.size(), DT1_(0));
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
                return x;
            }


            template <typename DT1_, typename DT2_>
                static DenseVector<DT1_> value(SparseMatrix<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, double konv_rad)
                {
                    CONTEXT("When solving sparse linear system with CG (with given convergence parameter):");


                    DenseVector<DT1_> x(right_hand_side.size(), DT1_(0));
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

    template <typename Tag_>
        struct ConjugateGradients<Tag_, methods::JAC>
        {
            private:
                template<typename DT1_, typename DT2_>
                    static inline void cg_kernel(DenseMatrix<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, DenseVector<DT1_> & former_gradient, DenseVector<DT1_> & former_result, DenseVector<DT1_> & utility, DenseVector<DT1_> & former_residual, DenseVector<DT1_> & diag_inverted)
                    {

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
                DT_ alpha, beta;
                beta = DotProduct<Tag_>::value(r, z);
                Product<Tag_>::value(temp_0, A, d);

                alpha = beta / DotProduct<Tag_>::value(d, temp_0);
                ScaledSum<Tag_>::value(x, d, alpha);
                ScaledSum<Tag_>::value(r, temp_0, -alpha);

                ElementProduct<Tag_>::value(z, dd_inverted, r);
                beta = DotProduct<Tag_>::value(z, r) / beta;

                copy<Tag_>(d, temp_0);
                ScaledSum<Tag_>::value(d, z, temp_0, beta);

                /*
                DT_ alpha, beta;

                Product<Tag_>::value(temp_0, A, d);

                //alpha = z_k^T * r_k / d_k^T A d_k
                z.lock(lm_read_only);
                r.lock(lm_read_only);
                temp_0.lock(lm_read_only);
                d.lock(lm_read_only);
                DT_ pre_alpha(DotProduct<Tag_>::value(temp_0, d));

                if(pre_alpha > std::numeric_limits<DT_>::epsilon())
                    alpha = DotProduct<Tag_>::value(z, r) / pre_alpha;
                else
                    alpha = DotProduct<Tag_>::value(z, r) / std::numeric_limits<DT_>::epsilon();
                z.unlock(lm_read_only);
                r.unlock(lm_read_only);
                temp_0.unlock(lm_read_only);
                d.unlock(lm_read_only);

                //x_k+1 = x_k + alpha z_k
                ScaledSum<Tag_>::value(x, d, alpha);                                               //STATUS:: STORE x_k+1

                z.lock(lm_read_only);
                r.lock(lm_read_only);
                beta = DotProduct<Tag_>::value(z, r);
                z.unlock(lm_read_only);
                r.unlock(lm_read_only);

                //r_k+1 = r_k - alpha A d_k
                ScaledSum<Tag_>::value(r, temp_0, -alpha);                                         //STATUS: STORE r_k+1

                //Preconditioner:
                //z_k+1 = C^-1 r_k+1
                copy<Tag_>(dd_inverted, temp_0);
                //temp_0 = dd_inverted.copy();
                ElementProduct<Tag_>::value(temp_0, r);                                            //STATUS: z_k+1 stored in temp_0

                //beta = z_k+1^T r_k+t / z_k^T r_k
                r.lock(lm_read_only);
                temp_0.lock(lm_read_only);
                if(beta > std::numeric_limits<DT_>::epsilon())
                    beta = DotProduct<Tag_>::value(temp_0, r) / beta;
                else
                    beta = DotProduct<Tag_>::value(temp_0, r) / std::numeric_limits<DT_>::epsilon();

                r.unlock(lm_read_only);
                temp_0.unlock(lm_read_only);
                copy<Tag_>(temp_0, z);                                                             //STATUS: STORE z_k+1
                //z = temp_0.copy();

                //d_k+1 = z_k+1 + beta d_k
                ScaledSum<Tag_>::value(temp_0, d, beta);
                copy<Tag_>(temp_0, d);                                                             //STATUS: STORE d_k+1
                //d = temp_0.copy();
                */
            }
            //end SMELL type

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
                return x;


            }
            /**
            * \brief Returns solution of LES given by a SparseMatrix and a Vector.
            *
            * \param system_matrix The system matrix.
            * \param right_hand_side The right hand side of the system.
            * \param konv_rad The parameter for convergence control.
            *
            */

            /// \{

            template <typename DT1_, typename DT2_>
            static DenseVector<DT1_> value(SparseMatrix<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, double konv_rad)
            {
                CONTEXT("When solving sparse linear system with PCG-Jacobi (with given convergence parameter):");


                DenseVector<DT1_> x(right_hand_side.size(), DT1_(0));
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
                return x;


            }


            template <typename DT_>
            static void value(SparseMatrixELL<DT_> & system_matrix,
                                           DenseVector<DT_> & right_hand_side,
                                           DenseVector<DT_> & x,
                                           DenseVector<DT_> & dd_inverted,
                                           unsigned long max_iters)
            {
                CONTEXT("When solving sparse ELL linear system with PCG-Jacobi: ");

                DenseVector<DT_> t_0(right_hand_side.size());

                DenseVector<DT_> r(right_hand_side.size());
                Defect<Tag_>::value(r, right_hand_side, system_matrix, x);
                DT_ initial_defect_norm(Norm<vnt_l_two, false, tags::GPU::CUDA>::value(r));
                std::cout << "Initial defect NORM: " << initial_defect_norm << std::endl;

                DenseVector<DT_> z(right_hand_side.size());
                ElementProduct<Tag_>::value(z, dd_inverted, r);
                DenseVector<DT_> d(z.size());
                copy<Tag_>(z , d);

                DenseVector<DT_> t_1(right_hand_side.size());
                for(unsigned long i(0) ; i < max_iters ; ++i)
                {
                    cg_kernel(system_matrix, r, z, d, x, dd_inverted, t_0);
                    Defect<Tag_>::value(t_1, right_hand_side, system_matrix, x);
                    DT_ current_defect_norm(Norm<vnt_l_two, false, tags::GPU::CUDA>::value(t_1));

                    if(current_defect_norm <= initial_defect_norm * 10e-8)
                    {
                        std::cout << "Converged after " << i + 1 << " iterations: NORM: " << current_defect_norm << std::endl;
                        break;
                    }
                    if(i == max_iters - 1)
                    {
                        std::cout << "NO convergence after " << i + 1 << " iterations: NORM: " << current_defect_norm << std::endl;
                    }
                }
                /*
                DenseVector<DT_> t_0(right_hand_side.size());
                DenseVector<DT_> t_1(right_hand_side.size());

                //r_0 = b - Ax_0
                DenseVector<DT_> r(Defect<Tag_>::value(right_hand_side, system_matrix, x));
                //DT_ initial_defect_norm(Norm<vnt_l_two, false, Tag_>::value(r));
                r.lock(lm_read_only);
                DT_ initial_defect_norm(Norm<vnt_l_two, false, tags::CPU>::value(r));
                r.unlock(lm_read_only);
                std::cout << "Initial defect NORM: " << initial_defect_norm << std::endl;

                //z_0 = C^-1 r_0
                DenseVector<DT_> z(right_hand_side.size());
                //copy<Tag_>(dd_inverted, z);
                z = dd_inverted.copy();
                ElementProduct<DT_>::value(z, r);

                //d_0 = z_0
                DenseVector<DT_> d(right_hand_side.size());
                copy<Tag_>(z, d);
                //d = z.copy();

                for(unsigned long i(0) ; i < max_iters ; ++i)
                {
                    cg_kernel(system_matrix, r, z, d, x, dd_inverted, t_0, t_1);
                    //DT_ current_defect_norm(Norm<vnt_l_two, false, Tag_>::value(Defect<Tag_>::value(right_hand_side, system_matrix, x)));
                    DenseVector<DT_> blub = Defect<Tag_>::value(right_hand_side, system_matrix, x);
                    blub.lock(lm_read_only);
                    DT_ current_defect_norm(Norm<vnt_l_two, false, tags::CPU>::value(blub));
                    blub.unlock(lm_read_only);
                    if(current_defect_norm <= initial_defect_norm * 10e-8)
                    {
                        std::cout << "Converged after " << i + 1 << " iterations: NORM: " << current_defect_norm << std::endl;
                        break;
                    }
                    if(i == max_iters - 1)
                    {
                        std::cout << "NO convergence after " << i + 1 << " iterations: NORM: " << current_defect_norm << std::endl;
                    }
                }*/

            }

        };
}

#endif
