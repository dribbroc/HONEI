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
        static SparseMatrix<DT_> value(SparseMatrix<DT_> & src)
        {
            /* SPAI parameters */
            double epsilon_param=0.6;
            int    nbsteps_param=5;
            int    max_param=100000;
            int    maxnew_param=5;
            int    cache_size_param=5;
            int    block_size_param=0;
            int    symmetric_pattern_param=0;
            int    left_precon_param=0;
            int    verbose_param=0;
            int    spar_param =0;
            int    lower_diag_param=0;
            int    upper_diag_param=0;
            double tau_param=0;
            int binary_param=0;

            std::string temp_file("spai_temp.mtx");
            // src schreiben in mm file
            MatrixIO<io_formats::MTX>::write_matrix(temp_file.c_str(), src);

            matrix *A = NULL;
            matrix *M = NULL;
            A = read_mm_matrix(temp_file.c_str(),
                    1,
                    1,
                    symmetric_pattern_param,
                    left_precon_param,
                    binary_param,
                    verbose_param,
                    NULL);

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

            write_matrix_mm(M,temp_file.c_str(),left_precon_param);

            // m einlesen in sparsematrix
            unsigned long non_zeros(MatrixIO<io_formats::MTX>::get_non_zeros(temp_file.c_str()));
            unsigned long rows, columns, ax, bx;
            DenseVector<unsigned long> r(non_zeros);
            DenseVector<unsigned long> c(non_zeros);
            DenseVector<DT_> data(non_zeros);

            MatrixIO<io_formats::MTX>::read_matrix(temp_file.c_str(), r, c, data);
            MatrixIO<io_formats::MTX>::get_sizes(temp_file.c_str(), rows, columns, ax, bx);
            SparseMatrix<DT_> tsmatrix(rows, columns, r, c, data);

            remove(temp_file.c_str());
            return tsmatrix;
        }
    };

}
#endif
