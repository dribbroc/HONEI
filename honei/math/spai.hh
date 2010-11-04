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
#ifndef LIBMATH_GUARD_SPAI_HH
#define LIBMATH_GUARD_SPAI_HH 1

#include <honei/util/tags.hh>
#include <honei/la/algorithm.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/math/matrix_io.hh>

extern "C"
{
#include <honei/math/spai/src/spai.h>
#include <honei/math/spai/src/read_mm_matrix.h>
}

#include <iostream>
#include <cstdio>


namespace honei
{

    struct SPAI
    {
        template <typename DT_>
        static SparseMatrix<DT_> value(const SparseMatrix<DT_> & src, double epsilon_param = 0.6, int nbsteps_param = 5, int maxnew_param = 5)
        {
            /* SPAI parameters */
            int    max_param=100000;
            int    cache_size_param=5;
            int    block_size_param=0;
            int    symmetric_pattern_param=0;
            int    left_precon_param=0;
            int    verbose_param=0;
            int    spar_param =0;
            int    lower_diag_param=0;
            int    upper_diag_param=0;
            double tau_param=0;

            matrix *A = NULL;
            matrix *M = NULL;

            unsigned long used_elements(src.used_elements());
            mm_data * rows_array = new mm_data[used_elements];
            mm_data * cols_array = new mm_data[used_elements];
            unsigned long tmp_i(0);

            for (typename SparseMatrix<DT_>::NonZeroConstElementIterator iter(src.begin_non_zero_elements()) ; iter != src.end_non_zero_elements() ; ++iter)
            {
                rows_array[tmp_i].i = iter.row();
                rows_array[tmp_i].j = iter.column();
                rows_array[tmp_i].val = *iter;
                cols_array[tmp_i].i = iter.row();
                cols_array[tmp_i].j = iter.column();
                cols_array[tmp_i].val = *iter;
                tmp_i++;
            }

            A = mm_to_matrix
                (left_precon_param,
                 symmetric_pattern_param,
                 src.rows(),
                 1,
                 src.rows(),
                 rows_array, used_elements,
                 cols_array, used_elements,
                 NULL);

            order_pointers(A);

            bspai
                (A, &M,
                 NULL,
                 epsilon_param,
                 nbsteps_param,
                 max_param,
                 maxnew_param,
                 block_size_param,
                 cache_size_param,
                 verbose_param,
                 spar_param,
                 lower_diag_param,
                 upper_diag_param,
                 tau_param);


            SparseMatrix<DT_> tsmatrix(src.rows(), src.columns());
            int ptr;
            for (long j=0; j<M->mnl; j++) {

                for (long i=0; i<M->lines->len[j]; i++) {
                    unsigned long row, col;

                    ptr = M->lines->ptrs[j][i];

                    if (! M->transposed) {
                        row = ptr;
                        col = j+M->my_start_index;
                    }
                    else {
                        col = ptr;
                        row = j+M->my_start_index;
                    }

                    if (M->bs == 1) {
                        tsmatrix(row, col) = M->lines->A[j][i];
                    }
                    else {
                        std::cout<<"Error: Dirk hat blocksize != 1 vergessen!"<<std::endl;
                        exit(1);
                    }
                }
            }
            delete [] rows_array;
            delete [] cols_array;
            return tsmatrix;
        }
    };
}
#endif
