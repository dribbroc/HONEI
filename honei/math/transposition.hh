/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the HONEI MATH C++ library. Math is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#pragma once
#ifndef MATH_GUARD_TRANSPOSITION_HH
#define MATH_GUARD_TRANSPOSITION_HH 1

#include <honei/la/sparse_matrix.hh>
#include <honei/mpi/operations.hh>
#include <honei/mpi/sparse_matrix_ell_mpi-fwd.hh>
#include <honei/mpi/sparse_matrix_csr_mpi-fwd.hh>

namespace honei
{
    template <typename Tag_>
    struct Transposition
    {
        public:
        template <typename MT_>
        static inline MT_ value(const MT_ & src)
        {
            CONTEXT("When transposing sparse matrix: ");
            SparseMatrix<typename MT_::iDT_> source(src);
            SparseMatrix<typename MT_::iDT_> target(src.columns(), src.rows());

            for(typename SparseMatrix<typename MT_::iDT_>::NonZeroConstElementIterator i(source.begin_non_zero_elements()) ; i != source.end_non_zero_elements() ; ++i)
            {
                target(i.column(), i.row(), *i);
            }

            MT_ result(target);

            return result;
        }

        template <typename DT_>
        static inline SparseMatrixELLMPI<DT_> value(const SparseMatrixELLMPI<DT_> & src)
        {
            return MPIOps<tags::CPU>::transposition(src);
        }

        template <typename DT_>
        static inline SparseMatrixCSRMPI<DT_> value(const SparseMatrixCSRMPI<DT_> & src)
        {
            return MPIOps<tags::CPU>::transposition(src);
        }
    };
}

#endif
