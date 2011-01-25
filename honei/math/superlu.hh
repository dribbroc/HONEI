/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de>
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
#ifndef LIBMATH_GUARD_SUPERLU_HH
#define LIBMATH_GUARD_SUPERLU_HH 1

#include <honei/util/tags.hh>
#include <honei/la/algorithm.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/math/matrix_io.hh>
#include "honei/math/SuperLU_4.1/SRC/slu_ddefs.h"

extern "C"
{
#include <honei/math/spai/src/spai.h>
#include <honei/math/spai/src/read_mm_matrix.h>
}

#include <iostream>
#include <cstdio>


namespace honei
{

    struct SuperLU
    {
        template <typename DT_>
        static void value(const SparseMatrixELL<DT_> & in_matrix, const DenseVector<DT_> & in_rhs, DenseVector<DT_> & out_result)
        {
            SuperMatrix A, L, U, B;
            double   *a, *rhs, *ta;
            int      *asub, *xa, *tasub, *txa;
            int      *perm_r; /* row permutations from partial pivoting */
            int      *perm_c; /* column permutation vector */
            int      nrhs, info, i, m, n, nnz;//, permc_spec;
            superlu_options_t options;
            SuperLUStat_t stat;

            /* Initialize matrix A. */
            m = in_matrix.rows();
            n = in_matrix.columns();
            nnz = in_matrix.used_elements();
            if ( !(a = doubleMalloc(nnz)) ) ABORT("Malloc fails for a[].");
            if ( !(asub = intMalloc(nnz)) ) ABORT("Malloc fails for asub[].");
            if ( !(xa = intMalloc(n+1)) ) ABORT("Malloc fails for xa[].");
            if ( !(ta = doubleMalloc(nnz)) ) ABORT("Malloc fails for ta[].");
            if ( !(tasub = intMalloc(nnz)) ) ABORT("Malloc fails for tasub[].");
            if ( !(txa = intMalloc(m+1)) ) ABORT("Malloc fails for txa[].");

            unsigned long ti(0);
            /*for (unsigned long ii(0) ; ii < in_matrix.rows() ; ++ii)
            {
                xa[ii] = ti;
                for (unsigned long ij(0) ; ij < in_matrix.columns() ; ++ij)
                {
                    if (in_matrix(ii, ij) != DT_(0))
                    {
                        a[ti] = in_matrix(ii, ij);
                        asub[ti] = ij;
                        ++ti;
                    }
                }
            }*/
            for (unsigned long srow(0) ; srow < in_matrix.rows() ; ++srow)
            {
                xa[srow] = ti;
                for (unsigned long ii(srow) ; ii < in_matrix.Ax().size() && ii/in_matrix.stride() < in_matrix.Arl()[srow] ; ii+=in_matrix.stride())
                {
                    a[ti] = in_matrix.Ax()[ii];
                    asub[ti] = in_matrix.Aj()[ii];
                    ++ti;
                }
            }


            xa[n] = nnz;

            dCompRow_to_CompCol(m, n, nnz, a, asub, xa, &ta, &tasub, &txa);

            /* Create matrix A in the format expected by SuperLU. */
            dCreate_CompRow_Matrix(&A, m, n, nnz, ta, tasub, txa, SLU_NC, SLU_D, SLU_GE);

            /* Create right-hand side matrix B. */
            nrhs = 1;
            if ( !(rhs = doubleMalloc(m * nrhs)) ) ABORT("Malloc fails for rhs[].");
            for (i = 0; i < m; ++i) rhs[i] = in_rhs[i];
            dCreate_Dense_Matrix(&B, m, nrhs, rhs, m, SLU_DN, SLU_D, SLU_GE);

            if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");
            if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");

            /* Set the default input options. */
            set_default_options(&options);
            options.ColPerm = NATURAL;

            /* Initialize the statistics variables. */
            StatInit(&stat);

            /* Solve the linear system. */
            dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);


            double *sol = (double*) ((DNformat*) B.Store)->nzval;
            for (i = 0; i < m; ++i) out_result[i] = sol[i];

            //dPrint_CompCol_Matrix("A", &A);
            //dPrint_CompCol_Matrix("U", &U);
            //dPrint_SuperNode_Matrix("L", &L);
            //print_int_vec("\nperm_r", m, perm_r);

            /* De-allocate storage */
            SUPERLU_FREE (rhs);
            SUPERLU_FREE (perm_r);
            SUPERLU_FREE (perm_c);
            SUPERLU_FREE (a);
            SUPERLU_FREE (asub);
            SUPERLU_FREE (xa);
            Destroy_CompCol_Matrix(&A);
            Destroy_SuperMatrix_Store(&B);
            Destroy_SuperNode_Matrix(&L);
            Destroy_CompCol_Matrix(&U);
            StatFree(&stat);
        }
    };
}
#endif
