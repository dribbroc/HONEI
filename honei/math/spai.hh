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

extern "C"
{
#include <honei/math/spai/src/spai.h>
#include <honei/math/spai/src/read_mm_matrix.h>
}

#include <iostream>
#include <cstdio>
#include <vector>


namespace honei
{

    struct SPAI
    {
        template <typename DT_>
        static SparseMatrix<DT_> value(const SparseMatrix<DT_> & src)
        {
            /* SPAI parameters */
            int    max_param=1000000;
            int    cache_size_param=5;
            int    block_size_param=1;
            int    symmetric_pattern_param=0;
            int    left_precon_param=0;
            int    verbose_param=0;
            int    spar_param =1; // 0: eps, 1: tau
            int    lower_diag_param=0;
            int    upper_diag_param=0;
            double tau_param=1;
            double epsilon_param = 0.6;
            int nbsteps_param = 5;
            int maxnew_param = 5;

            matrix *A = NULL;
            matrix *M = NULL;

            SparseMatrix<DT_> sm(src.copy());
            unsigned long root(sqrt(src.rows()));

            // 0
            for (unsigned long i(0) ; i < sm.rows() ; ++i)
            {
                sm(i, i, src(i, i));
                sm(i,i) += std::numeric_limits<double>::epsilon();
            }

            for (unsigned long i(0) ; i < sm.rows() - 1 ; ++i)
            {
                sm(i, i+1, src(i, i+1));
                sm(i,i+1) += std::numeric_limits<double>::epsilon();
            }

            for (unsigned long i(1) ; i < sm.rows() ; ++i)
            {
                sm(i, i-1, src(i, i-1));
                sm(i,i-1) += std::numeric_limits<double>::epsilon();
            }

            for (unsigned long i(0) ; i < sm.rows() - 2 ; ++i)
            {
                sm(i, i+2, src(i, i+2));
                sm(i,i+2) += std::numeric_limits<double>::epsilon();
            }

            for (unsigned long i(2) ; i < sm.rows() ; ++i)
            {
                sm(i, i-2, src(i, i-2));
                sm(i,i-2) += std::numeric_limits<double>::epsilon();
            }

            // root
            for (unsigned long i(0) ; i < sm.rows() - root ; ++i)
            {
                sm(i, i+root, src(i, i+root));
                sm(i,i+root) += std::numeric_limits<double>::epsilon();
            }

            for (unsigned long i(0) ; i < sm.rows() - root -1; ++i)
            {
                sm(i, i+(root-1), src(i, i+(root-1)));
                sm(i,i+(root-1)) += std::numeric_limits<double>::epsilon();
            }

            for (unsigned long i(0) ; i < sm.rows() - (root + 1); ++i)
            {
                sm(i, i+root+1, src(i, i+root+1));
                sm(i,i+root+1) += std::numeric_limits<double>::epsilon();
            }

            for (unsigned long i(0) ; i < sm.rows() - root -2; ++i)
            {
                sm(i, i+(root-2), src(i, i+(root-2)));
                sm(i,i+(root-2)) += std::numeric_limits<double>::epsilon();
            }

            for (unsigned long i(0) ; i < sm.rows() - (root + 2); ++i)
            {
                sm(i, i+root+2, src(i, i+root+2));
                sm(i,i+root+2) += std::numeric_limits<double>::epsilon();
            }

            // -root
            for (unsigned long i(root) ; i < sm.rows() ; ++i)
            {
                sm(i, i-root, src(i, i-root));
                sm(i,i-root) += std::numeric_limits<double>::epsilon();
            }

            for (unsigned long i(root-1) ; i < sm.rows() ; ++i)
            {
                sm(i, i-(root-1), src(i, i-(root-1)));
                sm(i,i-(root-1)) += std::numeric_limits<double>::epsilon();
            }

            for (unsigned long i(root+1) ; i < sm.rows() ; ++i)
            {
                sm(i, i-(root+1), src(i, i-(root+1)));
                sm(i,i-(root+1)) += std::numeric_limits<double>::epsilon();
            }

            for (unsigned long i(root-2) ; i < sm.rows() ; ++i)
            {
                sm(i, i-(root-2), src(i, i-(root-2)));
                sm(i,i-(root-2)) += std::numeric_limits<double>::epsilon();
            }

            for (unsigned long i(root+2) ; i < sm.rows() ; ++i)
            {
                sm(i, i-(root+2), src(i, i-(root+2)));
                sm(i,i-(root+2)) += std::numeric_limits<double>::epsilon();
            }

            // (root*2)
            for (unsigned long i(0) ; i < sm.rows() - (root*2) ; ++i)
            {
                sm(i, i+(root*2), src(i, i+(root*2)));
                sm(i,i+(root*2)) += std::numeric_limits<double>::epsilon();
            }

            for (unsigned long i(0) ; i < sm.rows() - (root*2) -1; ++i)
            {
                sm(i, i+((root*2)-1), src(i, i+((root*2)-1)));
                sm(i,i+((root*2)-1)) += std::numeric_limits<double>::epsilon();
            }

            for (unsigned long i(0) ; i < sm.rows() - ((root*2) + 1); ++i)
            {
                sm(i, i+(root*2)+1, src(i, i+(root*2)+1));
                sm(i,i+(root*2)+1) += std::numeric_limits<double>::epsilon();
            }

            for (unsigned long i(0) ; i < sm.rows() - (root*2) -2; ++i)
            {
                sm(i, i+((root*2)-2), src(i, i+((root*2)-2)));
                sm(i,i+((root*2)-2)) += std::numeric_limits<double>::epsilon();
            }

            for (unsigned long i(0) ; i < sm.rows() - ((root*2) + 2); ++i)
            {
                sm(i, i+(root*2)+2, src(i, i+(root*2)+2));
                sm(i,i+(root*2)+2) += std::numeric_limits<double>::epsilon();
            }

            // -(root*2)
            for (unsigned long i((root*2)) ; i < sm.rows() ; ++i)
            {
                sm(i, i-(root*2), src(i, i-(root*2)));
                sm(i,i-(root*2)) += std::numeric_limits<double>::epsilon();
            }

            for (unsigned long i((root*2)-1) ; i < sm.rows() ; ++i)
            {
                sm(i, i-((root*2)-1), src(i, i-((root*2)-1)));
                sm(i,i-((root*2)-1)) += std::numeric_limits<double>::epsilon();
            }

            for (unsigned long i((root*2)+1) ; i < sm.rows() ; ++i)
            {
                sm(i, i-((root*2)+1), src(i, i-((root*2)+1)));
                sm(i,i-((root*2)+1)) += std::numeric_limits<double>::epsilon();
            }

            for (unsigned long i((root*2)-2) ; i < sm.rows() ; ++i)
            {
                sm(i, i-((root*2)-2), src(i, i-((root*2)-2)));
                sm(i,i-((root*2)-2)) += std::numeric_limits<double>::epsilon();
            }

            for (unsigned long i((root*2)+2) ; i < sm.rows() ; ++i)
            {
                sm(i, i-((root*2)+2), src(i, i-((root*2)+2)));
                sm(i,i-((root*2)+2)) += std::numeric_limits<double>::epsilon();
            }

            std::vector<mm_data> elements;
            for (typename SparseMatrix<DT_>::NonZeroConstElementIterator iter(sm.begin_non_zero_elements()) ; iter != sm.end_non_zero_elements() ; ++iter)
            {
                mm_data temp;
                temp.i = iter.row();
                temp.j = iter.column();
                temp.val = *iter;
                elements.push_back(temp);
            }

            unsigned long used_elements(elements.size());
            mm_data * rows_array = new mm_data[used_elements];
            mm_data * cols_array = new mm_data[used_elements];
            for (unsigned long i(0) ; i < elements.size() ; ++i)
            {
                rows_array[i].i = elements.at(i).i;
                rows_array[i].j = elements.at(i).j;
                rows_array[i].val = elements.at(i).val;
                cols_array[i].i = elements.at(i).i;
                cols_array[i].j = elements.at(i).j;
                cols_array[i].val = elements.at(i).val;
            }

            A = mm_to_matrix
                (left_precon_param,
                 symmetric_pattern_param,
                 sm.rows(),
                 1,
                 sm.rows(),
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


            SparseMatrix<DT_> tsmatrix(sm.rows(), sm.columns());
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
