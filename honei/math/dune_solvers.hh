/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
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
#ifndef LIBMATH_GUARD_DUNE_SOLVERS_HH
#define LIBMATH_GUARD_DUNE_SOLVERS_HH 1

#include <honei/util/tags.hh>
#include <honei/math/methods.hh>
#include <honei/la/algorithm.hh>
#include <iostream>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#ifdef HAVE_SUPERLU
#include <dune/istl/superlu.hh>
#endif
#include <dune/istl/io.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

using namespace honei;
using namespace Dune;
using namespace methods;


namespace honei
{

    template<typename PreCon_>
    struct DuneConjugateGradients
    {
    };

    template <>
    struct DuneConjugateGradients<methods::JAC>
    {
        template <typename DT_>
            static void value(SparseMatrixELL<DT_> & system_matrix,
                    DenseVector<DT_> & right_hand_side,
                    DenseVector<DT_> & x,
                    DenseVector<DT_> & dd_inverted,
                    unsigned long max_iters,
                    DT_ eps_relative = 1e-8)
            {
#ifdef SOLVER_VERBOSE_L2
                std::cout << "Calling Dune CG solver, preconditioning=JACOBI, datalayout=ELL" << std::endl;
#endif

                typedef FieldMatrix<DT_, 1, 1> M;
                typedef FieldVector<DT_, 1> V;
                typedef BCRSMatrix<M> SM;
                typedef BlockVector<V> BV;

                SparseMatrix<DT_> smatrix2(system_matrix);
                SM smatrix_dune(smatrix2.rows(), smatrix2.columns(), SM::random);
                for (unsigned long i(0) ; i < smatrix2.rows() ; i++)
                {
                    smatrix_dune.setrowsize(i, smatrix2[i].used_elements());
                }
                smatrix_dune.endrowsizes();
                for (typename SparseMatrix<DT_>::NonZeroElementIterator i(smatrix2.begin_non_zero_elements()), i_end(smatrix2.end_non_zero_elements()) ; i != i_end ; ++i)
                {
                    smatrix_dune.addindex(i.row(), i.column());
                }
                smatrix_dune.endindices();
                for (typename SparseMatrix<DT_>::NonZeroElementIterator i(smatrix2.begin_non_zero_elements()), i_end(smatrix2.end_non_zero_elements()) ; i != i_end ; ++i)
                {
                    smatrix_dune[i.row()][i.column()] = *i;
                }
                BV result_dune(x.size(), 0);
                BV rhs_dune(right_hand_side.size(),0);
                for (unsigned long i(0) ; i < right_hand_side.size() ; ++i)
                {
                    rhs_dune[i] = right_hand_side[i];
                    result_dune[i] = x[i];
                }
                InverseOperatorResult irs;

                typedef SeqJac<SM, BV, BV> PREC;
                PREC prec(smatrix_dune, 1, 1);
                MatrixAdapter<SM,BV,BV> op(smatrix_dune);
#ifdef SOLVER_VERBOSE_L2
                int cgverbose = 2;
#else
                int cgverbose = 0;
#endif
                CGSolver<BV> cg(op, prec, eps_relative, max_iters, cgverbose);
                cg.apply(result_dune, rhs_dune, irs);

                for (unsigned long i(0) ; i < x.size() ; ++i)
                {
                    x[i] = result_dune[i];
                }

            }

    };

    struct DuneLU
    {
        template <typename DT_>
            static void value(SparseMatrixELL<DT_> & system_matrix,
                    DenseVector<DT_> & right_hand_side,
                    DenseVector<DT_> & x)
            {
#ifdef SOLVER_VERBOSE_L2
                std::cout << "Calling Dune LU solver, datalayout=ELL" << std::endl;
#endif

                typedef FieldMatrix<DT_, 1, 1> M;
                typedef FieldVector<DT_, 1> V;
                typedef BCRSMatrix<M> SM;
                typedef BlockVector<V> BV;

                SparseMatrix<DT_> smatrix2(system_matrix);
                SM smatrix_dune(smatrix2.rows(), smatrix2.columns(), SM::random);
                for (unsigned long i(0) ; i < smatrix2.rows() ; i++)
                {
                    smatrix_dune.setrowsize(i, smatrix2[i].used_elements());
                }
                smatrix_dune.endrowsizes();
                for (typename SparseMatrix<DT_>::NonZeroElementIterator i(smatrix2.begin_non_zero_elements()), i_end(smatrix2.end_non_zero_elements()) ; i != i_end ; ++i)
                {
                    smatrix_dune.addindex(i.row(), i.column());
                }
                smatrix_dune.endindices();
                for (typename SparseMatrix<DT_>::NonZeroElementIterator i(smatrix2.begin_non_zero_elements()), i_end(smatrix2.end_non_zero_elements()) ; i != i_end ; ++i)
                {
                    smatrix_dune[i.row()][i.column()] = *i;
                }
                BV result_dune(x.size(), 0);
                BV rhs_dune(right_hand_side.size(),0);
                for (unsigned long i(0) ; i < right_hand_side.size() ; ++i)
                {
                    rhs_dune[i] = right_hand_side[i];
                    result_dune[i] = x[i];
                }
                InverseOperatorResult irs;

#ifdef SOLVER_VERBOSE_L2
                bool luverbose = true;
#else
                bool luverbose = false;
#endif
                SuperLU<SM> slu(smatrix_dune, luverbose);
                slu.apply(result_dune, rhs_dune, irs);

                for (unsigned long i(0) ; i < x.size() ; ++i)
                {
                    x[i] = result_dune[i];
                }

            }

    };

}

#endif
