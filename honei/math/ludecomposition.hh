/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c)  2010 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
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
#ifndef MATH_GUARD_LU_DECOMPOSITION_HH
#define MATH_GUARD_LU_DECOMPOSITION_HH 1

#include <honei/la/dense_vector.hh>
#include <honei/la/dense_matrix.hh>

#include <iostream>

namespace honei
{
    template <typename Tag_> struct LUDecomposition;

    template <> struct LUDecomposition<tags::CPU>
    {
        public:
            template <typename DT_>
            static void value(DenseMatrix<DT_> & a, DenseVector<DT_> & b, DenseVector<DT_> & x)
            {
                if (a.rows() != a.columns())
                {
                    throw VectorSizeDoesNotMatch(a.rows(), a.columns());
                }
                if (a.rows() != b.size())
                {
                    throw VectorSizeDoesNotMatch(a.rows(), b.size());
                }
                if (a.rows() != x.size())
                {
                    throw VectorSizeDoesNotMatch(a.rows(), x.size());
                }

                DenseMatrix<DT_> u(a.copy());
                DenseMatrix<DT_> l(a.rows(), a.columns(), 0);
                for (unsigned long i(0) ; i < a.rows() ; ++i)
                {
                    l(i, i) = 1;
                }

                for (unsigned long k(0) ; k < u.rows() - 1 ; ++k)
                {
                    //search maximum pivot in column k
                    unsigned long pivot = k;
                    for (unsigned long t(k + 1) ; t < u.rows() ; ++t)
                    {
                        if (abs(u(t, k)) > abs(u(pivot, k)))
                                pivot = t;
                    }
                    //switch row k and row pivot
                    if (pivot != k)
                    {
                        for (unsigned long i(k) ; i < a.columns() ; ++i)
                        {
                            DT_ temp(u(k, i));
                            u(k, i) = u(pivot, i);
                            u(pivot, i) = temp;
                        }
                        for (unsigned long i(0) ; i < k  ; ++i)
                        {
                            DT_ temp(l(k, i));
                            l(k, i) = l(pivot, i);
                            l(pivot, i) = temp;
                        }
                        DT_ temp(b[k]);
                        b[k] = b[pivot];
                        b[pivot] = temp;
                    }

                    //todo calc and store LU insitu in A
                    for (unsigned long j(k + 1) ; j < u.rows() ; ++j)
                    {
                        l(j, k) = u(j, k) / u(k, k);
                        for (unsigned long i(k) ; i < u.rows() ; ++i)
                        {
                            u(j, i) = u(j, i) - l(j, k) * u(k, i);
                        }
                    }
                }

                for (unsigned long i(0) ; i < x.size() ; ++i)
                {
                    DT_ sum(0);
                    for (unsigned long j(0) ; j < i ; ++j)
                    {
                        sum += l(i,j) * x[j];
                    }
                    x[i] = DT_(1) / l(i, i) * (b[i] - sum);
                }

                for (long i(x.size() - 1) ; i >= 0 ; --i)
                {
                    DT_ sum(0);
                    for (unsigned long j(i+1) ; j < x.size() ; ++j)
                    {
                        sum += u(i,j) * x[j];
                    }
                    x[i] = DT_(1) / u(i, i) * (x[i] - sum);
                }
            }

            template <typename DT_>
            static void value(SparseMatrix<DT_> & a, DenseVector<DT_> & b, DenseVector<DT_> & x)
            {
                if (a.rows() != a.columns())
                {
                    throw VectorSizeDoesNotMatch(a.rows(), a.columns());
                }
                if (a.rows() != b.size())
                {
                    throw VectorSizeDoesNotMatch(a.rows(), b.size());
                }
                if (a.rows() != x.size())
                {
                    throw VectorSizeDoesNotMatch(a.rows(), x.size());
                }

                //todo use non zero elemment iterators

                SparseMatrix<DT_> u(a.copy());
                SparseMatrix<DT_> l(a.rows(), a.columns());
                for (unsigned long i(0) ; i < a.rows() ; ++i)
                {
                    l(i, i, 1);
                }

                for (unsigned long k(0) ; k < u.rows() - 1 ; ++k)
                {
                    //search maximum pivot in column k
                    unsigned long pivot = k;
                    DT_ piv_value = ((const SparseMatrix<DT_>)u)(pivot, k);
                    SparseVector<DT_> act_col(u.column(k));
                    for (typename SparseVector<DT_>::NonZeroConstElementIterator ci(act_col.begin_non_zero_elements()) ;
                            ci != act_col.end_non_zero_elements() ; ++ci)
                    {
                        if (ci.index() < k + 1)
                            continue;

                        if (abs(*ci) > abs(piv_value))
                        {
                                pivot = ci.index();
                                piv_value = *ci;
                        }
                    }

                    //switch row k and row pivot
                    if (pivot != k)
                    {
                        SparseVector<DT_> tempu = u[k];
                        u[k] = u[pivot];
                        u[pivot] = tempu;

                        SparseVector<DT_> templ = l[k];
                        l[k] = l[pivot];
                        l[pivot] = templ;

                        DT_ temp(b[k]);
                        b[k] = b[pivot];
                        b[pivot] = temp;
                    }
                }
                u._synch_column_vectors();
                l._synch_column_vectors();

                DT_ nzt;
                for (unsigned long k(0) ; k < u.rows() - 1 ; ++k)
                {
                    //todo calc and store LU insitu in A
                    DT_ ukk(((const SparseMatrix<DT_>)u)(k, k));
                    SparseVector<DT_> act_col(u.column(k));
                    for (typename SparseVector<DT_>::NonZeroConstElementIterator ci(act_col.begin_non_zero_elements()) ;
                            ci != act_col.end_non_zero_elements() ; ++ci)
                    {
                        if (ci.index() < k + 1)
                            continue;

                        l(ci.index(), k, *ci / ukk);
                    }

                    act_col = l.column(k);
                    for (typename SparseVector<DT_>::NonZeroConstElementIterator ci(act_col.begin_non_zero_elements()) ;
                            ci != act_col.end_non_zero_elements() ; ++ci)
                    {
                        if (ci.index() < k + 1)
                            continue;

                        SparseVector<DT_> act_row(u[k]);
                        for (typename SparseVector<DT_>::NonZeroConstElementIterator ri(act_row.begin_non_zero_elements()) ;
                                ri != act_row.end_non_zero_elements() ; ++ri)
                        {
                            if (ri.index() < k)
                                continue;

                            nzt = ((const SparseMatrix<DT_>)u)(ci.index(), ri.index()) - *ci * *ri;
                            if (nzt != DT_(0))
                                u(ci.index(), ri.index(), nzt);
                        }
                    }
                }

                for (unsigned long i(0) ; i < x.size() ; ++i)
                {
                    DT_ sum(0);
                    for (unsigned long j(0) ; j < i ; ++j)
                    {
                        sum += ((const SparseMatrix<DT_>)l)(i,j) * x[j];
                    }
                    x[i] = DT_(1) / ((const SparseMatrix<DT_>)l)(i, i) * (b[i] - sum);
                }

                for (long i(x.size() - 1) ; i >= 0 ; --i)
                {
                    DT_ sum(0);
                    for (unsigned long j(i+1) ; j < x.size() ; ++j)
                    {
                        sum += ((const SparseMatrix<DT_>)u)(i,j) * x[j];
                    }
                    x[i] = DT_(1) / ((const SparseMatrix<DT_>)u)(i, i) * (x[i] - sum);
                }

                unsigned long nza(0);
                unsigned long nzl(0);
                unsigned long nzu(0);
                for (unsigned long i(0) ; i < a.rows() ; ++i)
                {
                    nza += a[i].used_elements();
                    nzl += l[i].used_elements();
                    nzu += u[i].used_elements();
                }
                std::cout<<"NZ A: "<<nza<<" L: "<<nzl<<" U: "<<nzu<<std::endl;

            }
    };

    template <> struct LUDecomposition<tags::CPU::SSE>
    {
        public:
            template <typename DT_>
            static void value(DenseMatrix<DT_> & a, DenseVector<DT_> & b, DenseVector<DT_> & x)
            {
                if (a.rows() != a.columns())
                {
                    throw VectorSizeDoesNotMatch(a.rows(), a.columns());
                }
                if (a.rows() != b.size())
                {
                    throw VectorSizeDoesNotMatch(a.rows(), b.size());
                }
                if (a.rows() != x.size())
                {
                    throw VectorSizeDoesNotMatch(a.rows(), x.size());
                }

                DenseMatrix<DT_> u(a.copy());
                DenseMatrix<DT_> l(a.rows(), a.columns(), 0);
                DT_ * ue(u.elements());
                DT_ * le(l.elements());
                DT_ * be(b.elements());
                DT_ * xe(x.elements());
                const unsigned long rows(a.rows());
                const unsigned long columns(a.columns());
                const unsigned long size(x.size());
                for (unsigned long i(0) ; i < a.rows() ; ++i)
                {
                    le[columns * i + i] = 1;
                }

                for (unsigned long k(0) ; k < rows - 1 ; ++k)
                {
                    //todo calc and store LU insitu in A
                    for (unsigned long j(k + 1) ; j < rows ; ++j)
                    {
                        le[columns * j + k] = ue[columns * j +  k] / ue[columns * k + k];
                        for (unsigned long i(k) ; i < rows ; ++i)
                        {
                            ue[columns * j +  i] = ue[columns * j + i] - le[columns * j + k] * ue[columns * k + i];
                        }
                    }
                }

                for (unsigned long i(0) ; i < size ; ++i)
                {
                    DT_ sum(0);
                    for (unsigned long j(0) ; j < i ; ++j)
                    {
                        sum += le[columns * i +j] * xe[j];
                    }
                    xe[i] = DT_(1) / le[columns * i + i] * (be[i] - sum);
                }

                for (long i(size - 1) ; i >= 0 ; --i)
                {
                    DT_ sum(0);
                    for (unsigned long j(i+1) ; j < size ; ++j)
                    {
                        sum += ue[columns * i + j] * xe[j];
                    }
                    xe[i] = DT_(1) / ue[columns * i + i] * (xe[i] - sum);
                }
            }
    };
}
#endif
