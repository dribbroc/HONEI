/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2011 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the LA C++ library. LibLa is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibLa is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#pragma once
#ifndef MPI_GUARD_OPERATIONS_HH
#define MPI_GUARD_OPERATIONS_HH 1

#include <honei/mpi/dense_vector_mpi.hh>
#include <honei/mpi/sparse_matrix_ell_mpi.hh>

namespace honei
{
    template <typename Tag_>
        struct MPIOps
        {
            template <typename DT_>
                static void difference(DenseVectorMPI<DT_> & r, const DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y);

            template <typename DT_>
                static DT_ dot_product(const DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y);

            template <typename DT_>
                static void element_product(DenseVectorMPI<DT_> & r, const DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y);

            template <typename DT_>
                static DT_ norm_l2_false(const DenseVectorMPI<DT_> & x);

            template <typename DT_>
                static void product(DenseVectorMPI<DT_> & r, const SparseMatrixELLMPI<DT_> & a, const DenseVectorMPI<DT_> & b);

            template <typename DT_>
                static void scale(DenseVectorMPI<DT_> & x, DT_ a);

            template <typename DT_>
                static void scaled_sum(DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y, DT_ a);

            template <typename DT_>
                static void scaled_sum(DenseVectorMPI<DT_> & r, const DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y, DT_ a);

            template <typename DT_>
                static void sum(DenseVectorMPI<DT_> & x, const DenseVectorMPI<DT_> & y);
        };
}

#endif
