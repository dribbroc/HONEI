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

#ifndef LIBMATH_GUARD_JACOBI_HH
#define LIBMATH_GUARD_JACOBI_HH 1

#include <honei/util/tags.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/product.hh>
#include <honei/la/sum.hh>
#include <honei/la/difference.hh>
#include <honei/la/dot_product.hh>
#include <honei/la/scale.hh>
#include <honei/la/norm.hh>
#include <honei/la/element_inverse.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/la/algorithm.hh>
#include <honei/math/defect.hh>

namespace honei
{
    /**
     * \brief Solution of LES with Jacobi method.
     *
     * The referenced containers are invariant under this operation.
     * In every case, a new object is created and returned.
     *
     * \ingroup grpmatrixoperations
     * \ingroup grpvectoroperations
     */

    template <typename Tag_=tags::CPU>
    struct Jacobi
    {
        private:
            template<typename DT1_, typename DT2_>
            static inline void jacobi_kernel(DenseMatrix<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, DenseVector<DT1_> & former_result, DenseVector<DT1_> & diag, DenseVector<DT1_> & diag_inverted, DenseMatrix<DT1_> & difference)
            {
                DenseVector<DT1_> temp = Product<Tag_>::value(difference, former_result.copy());

                DenseVector<DT1_> temp2(right_hand_side.copy());

                Difference<Tag_>::value(temp2, temp);
                ElementProduct<Tag_>::value(temp2, diag_inverted);
                former_result = temp2.copy();
            }

            template<typename DT1_, typename DT2_>
            static inline void jacobi_kernel(BandedMatrix<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, DenseVector<DT1_> & former_result, DenseVector<DT1_> & diag, DenseVector<DT1_> & diag_inverted, BandedMatrix<DT1_> & difference)
            {
                DenseVector<DT1_> temp = Product<Tag_>::value(difference, former_result.copy());

                DenseVector<DT1_> temp2(right_hand_side.copy());

                Difference<Tag_>::value(temp2, temp);
                ElementProduct<Tag_>::value(temp2, diag_inverted);
                former_result = temp2;
            }


            template<typename DT1_, typename DT2_>
            static inline void jacobi_kernel(BandedMatrixQ1<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, DenseVector<DT1_> & former_result, DenseVector<DT1_> & diag_inverted, BandedMatrixQ1<DT1_> & difference)
            {
                DenseVector<DT1_> temp = Product<Tag_>::value(difference, former_result);

                DenseVector<DT1_> temp2(right_hand_side.size());
                copy<Tag_>(right_hand_side, temp2);

                Difference<Tag_>::value(temp2, temp);
                ElementProduct<Tag_>::value(temp2, diag_inverted);
                former_result = temp2;
            }
            template<typename DT1_, typename DT2_>
            static inline void jacobi_kernel(DenseVector<DT2_> & right_hand_side, DenseVector<DT1_> & former_result, DenseVector<DT1_> & diag_inverted, SparseMatrixELL<DT1_> & difference)
            {
                DenseVector<DT1_> temp(right_hand_side.size());
                Product<Tag_>::value(temp, difference, former_result);

                DenseVector<DT1_> temp2(right_hand_side.size());
                copy<Tag_>(right_hand_side, temp2);

                Difference<Tag_>::value(temp2, temp);
                ElementProduct<Tag_>::value(temp2, diag_inverted);

                former_result = temp2;
            }
//MG types:
            template<typename DT1_, typename DT2_>
            static inline void jacobi_kernel(DenseVector<DT1_> to_smooth, BandedMatrixQ1<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, DenseVector<DT1_> & former_result, DenseVector<DT1_> & diag_inverted, BandedMatrixQ1<DT1_> & difference, DT1_ omega)
            {
                //OLD HONEI:
                /*DenseVector<DT1_> temp(Product<Tag_>::value(difference, to_smooth));

                DenseVector<DT1_> temp2(right_hand_side.size());
                copy<Tag_>(right_hand_side, temp2);

                Difference<Tag_>::value(temp2, temp);
                ElementProduct<Tag_>::value(temp2, diag_inverted);
                Sum<Tag_>::value(temp2, to_smooth);
                former_result = temp2;*/

                //NEW:

                DenseVector<DT1_> temp(Defect<Tag_>::value(right_hand_side, system_matrix, to_smooth));

                former_result = to_smooth;
                ScaledSum<Tag_>::value(former_result, temp, diag_inverted);

                //OLD NONHONEI
                /*unsigned long n = right_hand_side.size();
                unsigned long root_n = (unsigned long)sqrt(n);

                to_smooth.lock(lm_read_only);
                right_hand_side.lock(lm_read_only);
                difference.lock(lm_read_only);
                former_result.lock(lm_write_only);

                DT2_ * rhs = right_hand_side.elements();
                DT1_ * x_old = to_smooth.elements();
                DT1_ * x_new = former_result.elements();

                DT1_ * ll = difference.band(LL).elements();
                DT1_ * ld = difference.band(LD).elements();
                DT1_ * lu = difference.band(LU).elements();

                DT1_ * dl = difference.band(DL).elements();
                DT1_ * dd = difference.band(DD).elements();
                DT1_ * du = difference.band(DU).elements();

                DT1_ * ul = difference.band(UL).elements();
                DT1_ * ud = difference.band(UD).elements();
                DT1_ * uu = difference.band(UU).elements();

                unsigned long i(0);
                //index 0
                x_new[i] = x_old[i] + ((rhs [i] - (dd[i] * x_old[i] +
                           du[i] * x_old[1] +
                           ul[i] * x_old[root_n - 1] +
                           ud[i] * x_old[root_n] +
                           uu[i] * x_old[root_n + 1])) * (omega / dd[i]));

                //index in [1, root_n -1[
                i = 1;
                for(; i < root_n - 1 ; ++i)
                {
                    x_new[i] = x_old[i] + ((rhs[i] - (dl[i] * x_old[i - 1] +
                               dd[i] * x_old[i] +
                               du[i] * x_old[i + 1] +
                               ul[i] * x_old[i + root_n - 1] +
                               ud[i] * x_old[i + root_n] +
                               uu[i] * x_old[i + root_n + 1])) * (omega/dd[i]));
                }

                //index root_n -1
                i = root_n - 1;
                x_new[i] = x_old[i] + ((rhs[i] - (lu[i] * x_old[i - (root_n - 1)] +
                           dl[i] * x_old[i - 1] +
                           dd[i] * x_old[i] +
                           du[i] * x_old[i + 1] +
                           ul[i] * x_old[i + root_n - 1] +
                           ud[i] * x_old[i + root_n] +
                           uu[i] * x_old[i + root_n + 1])) * (omega/dd[i]));

                //index root_n
                i = root_n;
                x_new[i] = x_old[i] + ((rhs[i] - (ld[i] * x_old[i - root_n] +
                           lu[i] * x_old[i - (root_n - 1)] +
                           dl[i] * x_old[i - 1] +
                           dd[i] * x_old[i] +
                           du[i] * x_old[ i + 1] +
                           ul[i] * x_old[i + root_n - 1] +
                           ud[i] * x_old[i + root_n] +
                           uu[i] * x_old[i + root_n + 1])) * (omega/dd[i]));

                //index in [root_n + 1, n - (root_n + 1)[
                i = root_n + 1;
                for(; i < n - (root_n  + 1) ; ++i)
                {
                    x_new[i] = x_old[i] + ((rhs[i] - (ll[i] * x_old[i - root_n - 1] +
                               ld[i] * x_old[i - root_n] +
                               lu[i] * x_old[i - root_n + 1] +
                               dl[i] * x_old[i - 1] +
                               dd[i] * x_old[i] +
                               du[i] * x_old[ i + 1] +
                               ul[i] * x_old[i + root_n - 1] +
                               ud[i] * x_old[i + root_n] +
                               uu[i] * x_old[i + root_n + 1])) * (omega/dd[i]));
                }

                //index n - (root_n + 1)
                i = n - (root_n + 1);
                x_new[i] = x_old[i] + ((rhs[i] - (ll[i] * x_old[i - (n - (root_n + 1))] +
                           ld[i] * x_old[i - root_n] +
                           lu[i] * x_old[i - root_n + 1] +
                           dl[i] * x_old[i - 1] +
                           dd[i] * x_old[i] +
                           du[i] * x_old[ i + 1] +
                           ul[i] * x_old[i + root_n - 1] +
                           ud[i] * x_old[i + root_n])) * (omega/dd[i]));

                //index n - root_n
                i = n - root_n;
                x_new[i] = x_old[i] + ((rhs[i] - (ll[i] * x_old[i - (n - (root_n + 1))] +
                           ld[i] * x_old[i - root_n] +
                           lu[i] * x_old[i - root_n + 1] +
                           dl[i] * x_old[i - 1] +
                           dd[i] * x_old[i] +
                           du[i] * x_old[ i + 1] +
                           ul[i] * x_old[i + root_n - 1])) * (omega/dd[i]));

                //index in [n - root_n + 1, n -1[
                i = n - root_n + 1;
                for(; i < n - 1; ++i)
                {
                    x_new[i] = x_old[i] + ((rhs[i] - (ll[i] * x_old[i - (n - (root_n + 1))] +
                               ld[i] * x_old[i - root_n] +
                               lu[i] * x_old[i - root_n + 1] +
                               dl[i] * x_old[i - 1] +
                               dd[i] * x_old[i] +
                               du[i] * x_old[ i + 1])) * (omega/dd[i]));
                }

                //index n - 1
                i = n - 1;
                x_new[i] = x_old[i] + ((rhs[i] - (ll[i] * x_old[i - (n - (root_n + 1))] +
                    ld[i] * x_old[i - root_n] +
                    lu[i] * x_old[i - root_n + 1] +
                    dl[i] * x_old[i - 1] +
                    dd[i] * x_old[i])) * (omega/dd[i]));

                to_smooth.unlock(lm_read_only);
                right_hand_side.unlock(lm_read_only);
                difference.unlock(lm_read_only);
                former_result.unlock(lm_write_only);*/
            }
            template<typename DT1_, typename DT2_>
            static inline void jacobi_kernel(DenseVector<DT1_> to_smooth, SparseMatrixELL<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, DenseVector<DT1_> & former_result, DenseVector<DT1_> & diag_inverted, SparseMatrixELL<DT1_> & difference, DT1_ omega)
            {
                //todo DIRK
                //TODO use Defect including result vector to avoid temp vector usage
                DenseVector<DT1_> temp(Defect<Tag_>::value(right_hand_side, system_matrix, to_smooth));
                former_result = to_smooth;
                ScaledSum<Tag_>::value(former_result, temp, diag_inverted);
            }

            template<typename DT1_, typename DT2_>
            static inline void jacobi_kernel(BandedMatrixQ1<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, DenseVector<DT1_> & former_result, DenseVector<DT1_> & diag_inverted)
            {
                DenseVector<DT1_> temp2(right_hand_side.size());
                copy<Tag_>(right_hand_side, temp2);

                ElementProduct<Tag_>::value(temp2, diag_inverted);
                former_result = temp2;
            }
//end MG types
            template<typename DT1_, typename DT2_>
            static inline void jacobi_kernel(SparseMatrix<DT1_> & system_matrix, DenseVector<DT2_> & right_hand_side, DenseVector<DT1_> & former_result, DenseVector<DT1_> & diag, DenseVector<DT1_> & diag_inverted, SparseMatrix<DT1_> & difference)
            {
                DenseVector<DT1_> temp = Product<Tag_>::value(difference, former_result.copy());

                DenseVector<DT1_> temp2(right_hand_side.copy());

                Difference<Tag_>::value(temp2, temp);
                ElementProduct<Tag_>::value(temp2, diag_inverted);
                former_result = temp2.copy();
            }

        public:

            //---------------------------------------Q1---------------------------------------------------
            /**
            * \brief Returns solution of LES with the Jacobi method given by a BandedMatrixQ1 and a Vector.
            *
            * \param system_matrix The system matrix.
            * \param right_hand_side The right hand side of the system.
            * \param iter_number The fixed number of iterations.
            *
            */
            template <typename DT1_, typename DT2_>
            static void value(BandedMatrixQ1<DT1_> & system_matrix,
                    DenseVector<DT2_> & right_hand_side,
                    DenseVector<DT2_> & x,
                    unsigned long iter_number)
            {
                CONTEXT("When solving banded linear system (Q1) with Jacobi (fixed # iterations):");
#ifdef SOLVER_VERBOSE_L2
                std::cout << "Calling JACOBI solver, datalayout=Q1" << std::endl;
#endif
                DenseVector<DT1_> diag(right_hand_side.size(), DT1_(0));

                DenseVector<DT1_> diag_inverted(right_hand_side.size(), DT1_(0));

                BandedMatrixQ1<DT1_> difference(system_matrix.size(), DenseVector<DT1_>(system_matrix.band(LL).size()), DenseVector<DT1_>(system_matrix.band(LL).size()), DenseVector<DT1_>(system_matrix.band(LL).size()), DenseVector<DT1_>(system_matrix.band(LL).size()), DenseVector<DT1_>(system_matrix.band(LL).size()), DenseVector<DT1_>(system_matrix.band(LL).size()), DenseVector<DT1_>(system_matrix.band(LL).size()), DenseVector<DT1_>(system_matrix.band(LL).size()), DenseVector<DT1_>(system_matrix.band(LL).size()));
                copy<Tag_>(system_matrix, difference);

                DenseVector<DT1_> zeros(right_hand_side.size(), DT1_(0));
                difference.lock(lm_read_and_write);
                difference.unlock(lm_read_and_write);
                difference.band(DD) = zeros;
                ///Create Diagonal, invert, compute difference on the fly.
                system_matrix.lock(lm_read_only);
                for(unsigned long i =0; i < diag.size(); ++i)
                {
                    diag[i] = system_matrix.band(DD)[i];
                    if(fabs(diag[i]) >= std::numeric_limits<DT1_>::epsilon())
                    {
                        diag_inverted[i] = DT1_(1) / diag[i];
                    }
                    else
                    {
                        diag_inverted[i] = DT1_(1) / std::numeric_limits<DT1_>::epsilon();
                    }
                }
                system_matrix.unlock(lm_read_only);

                copy<Tag_>(right_hand_side, x);

                for(unsigned long i = 0; i<iter_number; ++i)
                {
                    jacobi_kernel(system_matrix, right_hand_side, x, diag_inverted, difference);
                }
            }

            template <typename DT1_, typename DT2_>
            static void value(BandedMatrixQ1<DT1_> & system_matrix,
                    DenseVector<DT2_> & right_hand_side,
                    DenseVector<DT2_> & x,
                    unsigned long iter_number,
                    DT1_ omega)
            {
                CONTEXT("When solving banded linear system (Q1) with Jacobi (fixed # iterations):");
#ifdef SOLVER_VERBOSE_L2
                std::cout << "Calling JACOBI solver, datalayout=Q1" << std::endl;
#endif
                DenseVector<DT1_> diag_inverted(right_hand_side.size());

                BandedMatrixQ1<DT1_> difference(system_matrix.size(), DenseVector<DT1_>(system_matrix.band(LL).size()), DenseVector<DT1_>(system_matrix.band(LL).size()), DenseVector<DT1_>(system_matrix.band(LL).size()), DenseVector<DT1_>(system_matrix.band(LL).size()), DenseVector<DT1_>(system_matrix.band(LL).size()), DenseVector<DT1_>(system_matrix.band(LL).size()), DenseVector<DT1_>(system_matrix.band(LL).size()), DenseVector<DT1_>(system_matrix.band(LL).size()), DenseVector<DT1_>(system_matrix.band(LL).size()));
                copy<Tag_>(system_matrix, difference);

                DenseVector<DT1_> zeros(right_hand_side.size(), DT1_(0));
                difference.lock(lm_read_and_write);
                difference.unlock(lm_read_and_write);
                difference.band(DD) = zeros;
                ///Create Diagonal, invert, compute difference on the fly.
                copy<Tag_>(system_matrix.band(DD), diag_inverted);
                ElementInverse<Tag_>::value(diag_inverted);
                Scale<Tag_>::value(diag_inverted, omega);

                copy<Tag_>(right_hand_side, x);

                for(unsigned long i = 0; i<iter_number; ++i)
                {
                    jacobi_kernel(system_matrix, right_hand_side, x, diag_inverted, difference);
                }
            }
//MG types:
            template <typename DT1_, typename DT2_>
            static inline void value(BandedMatrixQ1<DT1_> & system_matrix,
                    DenseVector<DT2_> & right_hand_side,
                    DenseVector<DT2_> & x,
                    DT1_ omega,
                    DenseVector<DT1_> & diag_inverted)
            {
                CONTEXT("When solving banded linear system (Q1) with Jacobi (fixed # iterations):");
#ifdef SOLVER_VERBOSE_L2
                std::cout << "Calling JACOBI smoother, datalayout=Q1, MULTIGRID (1)" << std::endl;
#endif
                /*DenseVector<DT1_> diag_inverted(right_hand_side.size());
                copy<Tag_>(system_matrix.band(DD), diag_inverted);
                ElementInverse<Tag_>::value(diag_inverted);
                Scale<Tag_>::value(diag_inverted, omega);
                DenseVector<DT1_> x(right_hand_side.size());
                x.lock(lm_write_only, Tag_::memory_value);
                x.unlock(lm_write_only);

                jacobi_kernel(system_matrix, right_hand_side, x, diag_inverted);
                return x;*/

                /*DenseVector<DT1_> x(right_hand_side.size());
                x.lock(lm_write_only);
                system_matrix.lock(lm_read_only);
                right_hand_side.lock(lm_read_only);

                for(unsigned long i(0) ; i < diag_inverted.size() ; ++i)
                {
                    x[i] = diag_inverted[i] * right_hand_side[i];
                }

                x.unlock(lm_write_only);
                system_matrix.unlock(lm_read_only);
                right_hand_side.unlock(lm_read_only);
                return x;*/
                //NEW:
                ElementProduct<Tag_>::value(x, diag_inverted, right_hand_side);
            }

            template <typename DT1_, typename DT2_>
            static inline void value(DenseVector<DT1_>& to_smooth,
                    BandedMatrixQ1<DT1_> & system_matrix,
                    DenseVector<DT2_> & right_hand_side,
                    DenseVector<DT2_> & x,
                    unsigned long iter_number,
                    DT1_ omega,
                    DenseVector<DT1_> & diag_inverted)
            {
                CONTEXT("When solving banded linear system (Q1) with Jacobi (fixed # iterations):");
#ifdef SOLVER_VERBOSE_L2
                std::cout << "Calling JACOBI smoother, datalayout=Q1, MULTIGRID (2)" << std::endl;
#endif

                    /*DenseVector<DT1_> diag_inverted(right_hand_side.size());
                    copy<Tag_>(system_matrix.band(DD), diag_inverted);
                    ElementInverse<Tag_>::value(diag_inverted);
                    Scale<Tag_>::value(diag_inverted, omega);
                    DenseVector<DT1_> x(right_hand_side.size());

                    x.lock(lm_write_only, Tag_::memory_value);
                    x.unlock(lm_write_only);

                    for(unsigned long i = 0; i<iter_number; ++i)
                    {
                        jacobi_kernel(to_smooth, system_matrix, right_hand_side, x, diag_inverted, system_matrix, omega);
                        DenseVector<DT1_> ts_c(to_smooth.size());
                        copy<Tag_>(to_smooth, ts_c);
                        ts_c = to_smooth;
                        to_smooth = x;
                        x = ts_c;
                    }
                    if(iter_number % 2 != 0)
                        return x;
                    else
                        return to_smooth;*/

                for(unsigned long i = 0; i<iter_number; ++i)
                {
                    jacobi_kernel(to_smooth, system_matrix, right_hand_side, x, diag_inverted, system_matrix, omega);
                }
            }

            template <typename DT1_, typename DT2_>
            static inline void value(SparseMatrixELL<DT1_> & system_matrix,
                    DenseVector<DT2_> & right_hand_side,
                    DenseVector<DT2_> & x,
                    DT1_ omega,
                    DenseVector<DT1_> & diag_inverted)
            {
                CONTEXT("When smoothing sparse linear system with Jacobi (fixed # iterations):");
#ifdef SOLVER_VERBOSE_L2
                std::cout << "Calling JACOBI smoother, datalayout=ELLPACK, MULTIGRID (1)" << std::endl;
#endif
                ElementProduct<Tag_>::value(x, diag_inverted, right_hand_side);
            }

            template <typename DT1_, typename DT2_>
            static inline void value(DenseVector<DT1_>& to_smooth,
                    SparseMatrixELL<DT1_> & system_matrix,
                    DenseVector<DT2_> & right_hand_side,
                    DenseVector<DT2_> & x,
                    unsigned long iter_number,
                    DT1_ omega,
                    DenseVector<DT1_> & diag_inverted)
            {
                CONTEXT("When solving sparse linear system (ELL) with Jacobi (fixed # iterations):");
#ifdef SOLVER_VERBOSE_L2
                std::cout << "Calling JACOBI smoother, datalayout=ELLPACK, MULTIGRID (2)" << std::endl;
#endif
                for(unsigned long i = 0; i<iter_number; ++i)
                {
                    jacobi_kernel(to_smooth, system_matrix, right_hand_side, x, diag_inverted, system_matrix, omega);
                }
            }
//end MG types

            /**
            * \brief Returns solution of LES with the Jacobi method given by a DenseMatrix and a Vector.
            *
            * \param system_matrix The system matrix.
            * \param right_hand_side The right hand side of the system.
            * \param x The result.
            * \param iter_number The fixed number of iterations.
            *
            */
            template <typename DT1_, typename DT2_>
            static void value(DenseMatrix<DT1_> & system_matrix,
                    DenseVector<DT2_> & right_hand_side,
                    DenseVector<DT2_> & x,
                    unsigned long iter_number)
            {
                CONTEXT("When solving dense linear system with Jacobi (fixed # iterations):");
#ifdef SOLVER_VERBOSE_L2
                std::cout << "Calling JACOBI solver, datalayout=DENSE" << std::endl;
#endif
                DenseVector<DT1_> diag(right_hand_side.size(), DT1_(0));

                DenseVector<DT1_> diag_inverted(right_hand_side.size(), DT1_(0));

                DenseMatrix<DT1_> difference(system_matrix.copy());
                ///Create Diagonal, invert, compute difference on the fly.
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
                    difference[i][i] = DT1_(0);
                }

                //DenseVector<DT1_> x(right_hand_side.copy());

                for(unsigned long i = 0; i<iter_number; ++i)
                {
                    jacobi_kernel(system_matrix, right_hand_side, x, diag, diag_inverted, difference);
                }
            }

            template <typename DT1_, typename DT2_>
            static void value(DenseMatrix<DT1_> & system_matrix,
                    DenseVector<DT2_> & right_hand_side,
                    DenseVector<DT2_> & x,
                    double konv_rad)
            {
                CONTEXT("When solving dense linear system with Jacobi (with given convergence parameter):");
#ifdef SOLVER_VERBOSE_L2
                std::cout << "Calling JACOBI solver, datalayout=DENSE" << std::endl;
#endif
                DenseVector<DT1_> x_last(x.copy());

                DT1_ norm_x_last = DT1_(0);
                DT1_ norm_x = DT1_(1);
                DenseVector<DT1_> diag(right_hand_side.size(), DT1_(0));

                DenseVector<DT1_> diag_inverted(right_hand_side.size(), DT1_(0));

                DenseMatrix<DT1_> difference(system_matrix.copy());
                ///Create Diagonal, invert, compute difference on the fly.
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
                    difference[i][i] = DT1_(0);
                }

                while(fabs(norm_x - norm_x_last) > konv_rad)
                {
                    jacobi_kernel(system_matrix, right_hand_side, x, diag, diag_inverted, difference);
                    norm_x = Norm<vnt_l_two, false, Tag_>::value(x);
                    norm_x_last = Norm<vnt_l_two, false, Tag_>::value(x_last);
                    x_last = x.copy();
                }
            }

            /**
             * \brief Returns solution of LES with the Jacobi method given by a BandedMatrix and a Vector.
             *
             * \param system_matrix The system matrix.
             * \param right_hand_side The right hand side of the system.
             * \param iter_number The fixed number of iterations.
             *
             */
            template <typename DT1_, typename DT2_>
            static void value(BandedMatrix<DT1_> & system_matrix,
                    DenseVector<DT2_> & right_hand_side,
                    DenseVector<DT2_> & x,
                    unsigned long iter_number)
            {
                CONTEXT("When solving banded linear system with Jacobi (fixed # iterations):");
#ifdef SOLVER_VERBOSE_L2
                std::cout << "Calling JACOBI solver, datalayout=BANDED" << std::endl;
#endif
                DenseVector<DT1_> diag(right_hand_side.size(), DT1_(0));

                DenseVector<DT1_> diag_inverted(right_hand_side.size(), DT1_(0));

                BandedMatrix<DT1_> difference(system_matrix.copy());

                DenseVector<DT1_> zeros(right_hand_side.size(), DT1_(0));
                difference.insert_band(0, zeros);
                ///Create Diagonal, invert, compute difference on the fly.
                for(unsigned long i =0; i < diag.size(); ++i)
                {
                    diag[i] = system_matrix.band(0)[i];
                    if(fabs(diag[i]) >= std::numeric_limits<DT1_>::epsilon())
                    {
                        diag_inverted[i] = DT1_(1) / diag[i];
                    }
                    else
                    {
                        diag_inverted[i] = DT1_(1) / std::numeric_limits<DT1_>::epsilon();
                    }
                }
                //DenseVector<DT1_> x(right_hand_side.copy());
                for(unsigned long i = 0; i<iter_number; ++i)
                {
                    jacobi_kernel(system_matrix, right_hand_side, x, diag, diag_inverted, difference);
                }
            }


            template <typename DT1_, typename DT2_>
            static void value(BandedMatrix<DT1_> & system_matrix,
                    DenseVector<DT2_> & right_hand_side,
                    DenseVector<DT2_> & x,
                    double konv_rad)
            {
                CONTEXT("When solving banded linear system with Jacobi (with given convergence parameter):");
#ifdef SOLVER_VERBOSE_L2
                std::cout << "Calling JACOBI solver, datalayout=BANDED" << std::endl;
#endif
                DenseVector<DT1_> x_last(x.copy());

                DT1_ norm_x_last = DT1_(0);
                DT1_ norm_x = DT1_(1);
                DenseVector<DT1_> diag(right_hand_side.size(), DT1_(0));

                DenseVector<DT1_> diag_inverted(right_hand_side.size(), DT1_(0));

                BandedMatrix<DT1_> difference(system_matrix.copy());
                ///Create Diagonal, invert, compute difference on the fly.

                DenseVector<DT1_> zeros(right_hand_side.size(), DT1_(0));
                difference.insert_band(0, zeros);
                for(unsigned long i =0; i < diag.size(); ++i)
                {
                    diag[i] = system_matrix.band(0)[i];
                    if(fabs(diag[i]) >= std::numeric_limits<DT1_>::epsilon())
                    {
                        diag_inverted[i] = DT1_(1) / diag[i];
                    }
                    else
                    {
                        diag_inverted[i] = DT1_(1) / std::numeric_limits<DT1_>::epsilon();
                    }
                }

                while(fabs(norm_x - norm_x_last) > konv_rad)
                {
                    jacobi_kernel(system_matrix, right_hand_side, x, diag, diag_inverted, difference);
                    norm_x = Norm<vnt_l_two, false, Tag_>::value(x);
                    norm_x_last = Norm<vnt_l_two, false, Tag_>::value(x_last);
                    x_last = x.copy();
                }
            }

            template <typename DT1_, typename DT2_>
            static void value(SparseMatrixELL<DT1_> & system_matrix,
                    SparseMatrixELL<DT1_> & difference,
                    DenseVector<DT2_> & right_hand_side,
                    DenseVector<DT1_> & x,
                    DenseVector<DT1_> & diag_inverted,
                    unsigned long max_iters,
                    unsigned long & used_iters,
                    DT1_ eps_relative = 1e-8)
            {
                CONTEXT("When solving sparse linear system (ELL) with Jacobi:");
#ifdef SOLVER_VERBOSE_L2
                std::cout << "Calling JACOBI solver, datalayout=ELLPACK" << std::endl;
#endif
                DenseVector<DT1_> defect(right_hand_side.size());
                Defect<Tag_>::value(defect, right_hand_side, system_matrix, x);
                DT1_ initial_defect_norm(Norm<vnt_l_two, false, Tag_>::value(defect));
                DT1_ current_defect_norm(1000);

                for(unsigned long i = 0; i<max_iters; ++i)
                {
                    jacobi_kernel(right_hand_side, x, diag_inverted, difference);
                    Defect<Tag_>::value(defect, right_hand_side, system_matrix, x);
                    current_defect_norm = Norm<vnt_l_two, false, Tag_>::value(defect);

                    if(current_defect_norm < initial_defect_norm * eps_relative)
                    {
                        std::cout << "Converged after " << i + 1 << " iterations: NORM: " << current_defect_norm << std::endl;
                        used_iters = i + 1;
                        break;
                    }
                    if(i == max_iters - 1)
                    {
                        std::cout << "NO convergence after " << i + 1 << " iterations: NORM: " << current_defect_norm << std::endl;
                        used_iters = i + 1;
                    }
                }
            }

            /**
             * \brief Returns solution of LES with the Jacobi method given by a SparseMatrix and a Vector.
             *
             * \param system_matrix The system matrix.
             * \param right_hand_side The right hand side of the system.
             * \param x The solution (and initial guess)
             * \param iter_number The fixed number of iterations.
             *
             */
            template <typename DT1_, typename DT2_>
                static void value(SparseMatrix<DT1_> & system_matrix,
                        DenseVector<DT2_> & right_hand_side,
                        DenseVector<DT2_> & x,
                        unsigned long iter_number)
                {
                    CONTEXT("When solving sparse linear system with Jacobi (fixed # iterations):");
#ifdef SOLVER_VERBOSE_L2
                    std::cout << "Calling JACOBI solver, datalayout=SPARSE" << std::endl;
#endif
                    DenseVector<DT1_> diag(right_hand_side.size(), DT1_(0));

                    DenseVector<DT1_> diag_inverted(right_hand_side.size(), DT1_(0));

                    SparseMatrix<DT1_> difference(system_matrix.copy());
                    ///Create Diagonal, invert, compute difference on the fly.
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
                        difference[i][i] = DT1_(0);
                    }

                    //DenseVector<DT1_> x(right_hand_side.copy());

                    for(unsigned long i = 0; i<iter_number; ++i)
                    {
                        jacobi_kernel(system_matrix, right_hand_side, x, diag, diag_inverted, difference);
                    }
                }

            template <typename DT1_, typename DT2_>
                static void value(SparseMatrix<DT1_> & system_matrix,
                        DenseVector<DT2_> & right_hand_side,
                        DenseVector<DT2_> & x,
                        double konv_rad)
                {
                    CONTEXT("When solving sparse linear system with Jacobi (with given convergence parameter):");
#ifdef SOLVER_VERBOSE_L2
                    std::cout << "Calling JACOBI solver, datalayout=SPARSE" << std::endl;
#endif
                    DenseVector<DT1_> x_last(x.copy());

                    DT1_ norm_x_last = DT1_(0);
                    DT1_ norm_x = DT1_(1);
                    DenseVector<DT1_> diag(right_hand_side.size(), DT1_(0));

                    DenseVector<DT1_> diag_inverted(right_hand_side.size(), DT1_(0));

                    SparseMatrix<DT1_> difference(system_matrix.copy());
                    ///Create Diagonal, invert, compute difference on the fly.
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
                        difference[i][i] = DT1_(0);
                    }

                    while(fabs(norm_x - norm_x_last) > konv_rad)
                    {
                        jacobi_kernel(system_matrix, right_hand_side, x, diag, diag_inverted, difference);
                        norm_x = Norm<vnt_l_two, false, Tag_>::value(x);
                        norm_x_last = Norm<vnt_l_two, false, Tag_>::value(x_last);
                        x_last = x.copy();
                    }
                }
    };
}
#endif
