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
#ifndef LIBMATH_GUARD_SPAI2_HH
#define LIBMATH_GUARD_SPAI2_HH 1

#include <honei/util/tags.hh>
#include <honei/la/sparse_matrix.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/util/operation_wrapper.hh>
#include <honei/util/configuration.hh>
#include <honei/backends/multicore/operation.hh>
#include <honei/backends/multicore/thread_pool.hh>
#include <vector>
#include <algorithm>
#include <iostream>
#include <mkl_lapacke.h>

//based on "Parallel Preconditioning with Sparse Approximate Inverses" by Grote et al.
namespace honei
{
    template <typename Tag_>
    struct SPAI2
    {
    };

    template <>
    struct SPAI2<tags::CPU>
    {
        template <typename DT_>
        static SparseMatrix<DT_> & value(SparseMatrix<DT_> & M, const SparseMatrix<DT_> & A, unsigned long col_start = 0, unsigned long col_end = 0)
        {
            if (col_end == 0)
                col_end = A.columns();
            for (unsigned long idx(col_start) ; idx < col_end ; ++idx)
            {
                unsigned long n2(A.column(idx).used_elements());
                if (n2 == 0)
                    continue;
                unsigned long J[n2];
                for (unsigned long i(0) ; i < n2 ; ++i)
                {
                    J[i] = A.column(idx).indices()[i];
                }

                unsigned long n1(0);
                for (unsigned long i(0) ; i < n2 ; ++i)
                {
                    n1 += A.column(J[i]).used_elements();
                }
                if (n1 == 0)
                    continue;
                unsigned long I[n1];
                DenseVector<DT_> et(n1, DT_(0));
                unsigned long tmp(0);
                for (unsigned long i(0) ; i < n2 ; ++i)
                {
                    for (unsigned long j(0) ; j < A.column(J[i]).used_elements() ; ++j)
                    {
                        I[tmp] = A.column(J[i]).indices()[j];
                        et[tmp] = (I[tmp] == idx);
                        ++tmp;
                    }
                }

                DenseMatrix<DT_> At(n1, n2, DT_(0));
                for (unsigned long i(0) ; i < n1 ; ++i)
                {
                    SparseVector<DT_> row = A[I[i]];
                    unsigned long it(0);
                    const unsigned long * indices(row.indices());
                    const DT_ * elements(row.elements());
                    const unsigned long used_elements(row.used_elements());
                    for (unsigned long j(0) ; j < n2 ; ++j)
                    {
                        const unsigned long index(J[j]);
                        for ( ; (indices[it] < index) && (it < used_elements) ; ++it)
                            ;
                        if (indices[it] == index)
                            At(i, j) = elements[it];
                    }
                }

                LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', At.rows(), At.columns(), 1, At.elements(), At.columns(), et.elements(), 1);

                for (unsigned long i(0) ; i < n2 ; ++i)
                {
                    M(A.column(idx).indices()[i], idx, et[i]);
                }
            }

            return M;
        }

        template <typename DT_>
        static DenseVector<DT_> & value_0(DenseVector<DT_> & result, const SparseMatrix<DT_> & A, unsigned long col_start = 0, unsigned long col_end = 0)
        {
            if (col_end == 0)
                col_end = A.columns();
            for (unsigned long idx(col_start) ; idx < col_end ; ++idx)
            {
                unsigned long n2(1);
                unsigned long J[n2];
                J[0] = idx;

                unsigned long n1(0);
                for (unsigned long i(0) ; i < n2 ; ++i)
                {
                    n1 += A.column(J[i]).used_elements();
                }
                if (n1 == 0)
                    continue;
                unsigned long I[n1];
                DenseVector<DT_> et(n1, DT_(0));
                unsigned long tmp(0);
                for (unsigned long i(0) ; i < n2 ; ++i)
                {
                    for (unsigned long j(0) ; j < A.column(J[i]).used_elements() ; ++j)
                    {
                        I[tmp] = A.column(J[i]).indices()[j];
                        et[tmp] = (I[tmp] == idx);
                        ++tmp;
                    }
                }

                DenseMatrix<DT_> At(n1, n2, DT_(0));
                for (unsigned long i(0) ; i < n1 ; ++i)
                {
                    SparseVector<DT_> row = A[I[i]];
                    unsigned long it(0);
                    const unsigned long * indices(row.indices());
                    const DT_ * elements(row.elements());
                    const unsigned long used_elements(row.used_elements());
                    for (unsigned long j(0) ; j < n2 ; ++j)
                    {
                        const unsigned long index(J[j]);
                        for ( ; (indices[it] < index) && (it < used_elements) ; ++it)
                            ;
                        if (indices[it] == index)
                            At(i, j) = elements[it];
                    }
                }

                LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', At.rows(), At.columns(), 1, At.elements(), At.columns(), et.elements(), 1);

                result[idx] = et[0];
            }

            return result;
        }
    };

    namespace mc
    {
        template <typename Tag_> struct SPAI2
        {
            template <typename DT_>
            static SparseMatrix<DT_> & value(SparseMatrix<DT_> & M, const SparseMatrix<DT_> & A)
            {
                unsigned long max_count(Configuration::instance()->get_value("mc::Product(DV,SMELL,DV)::max_count",
                            mc::ThreadPool::instance()->num_threads()));

                TicketVector tickets;

                unsigned long limits[max_count + 1];
                limits[0] = 0;
                for (unsigned long i(1) ; i < max_count; ++i)
                {
                    limits[i] = limits[i-1] + A.columns() / max_count;
                }
                limits[max_count] = A.columns();

                for (unsigned long i(0) ; i < max_count ; ++i)
                {
                    OperationWrapper<honei::SPAI2<typename Tag_::DelegateTo>, SparseMatrix<DT_>,
                        SparseMatrix<DT_>, SparseMatrix<DT_>, unsigned long, unsigned long > wrapper(M);
                    tickets.push_back(mc::ThreadPool::instance()->enqueue(bind(wrapper, M, A, limits[i], limits[i+1])));
                }

                tickets.wait();

                return M;
            }
        };
    }

    template <> struct SPAI2<tags::CPU::MultiCore> :
        public mc::SPAI2<tags::CPU::MultiCore>
    {
    };
}
#endif
