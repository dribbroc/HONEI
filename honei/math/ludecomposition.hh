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


#ifndef MATH_GUARD_LU_DECOMPOSITION_HH
#define MATH_GUARD_LU_DECOMPOSITION_HH 1

#include <honei/la/dense_vector.hh>
#include <honei/la/dense_matrix.hh>

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
                        if (fabs(u(t, k)) > fabs(u(pivot, k)))
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

                SparseMatrix<DT_> u(a.copy());
                SparseMatrix<DT_> l(a.rows(), a.columns());
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
                        if (fabs(((const SparseMatrix<DT_>)u)(t, k)) > fabs(((const SparseMatrix<DT_>)u)(pivot, k)))
                                pivot = t;
                    }
                    //switch row k and row pivot
                    DT_ nzt;
                    if (pivot != k)
                    {
                        for (unsigned long i(k) ; i < a.columns() ; ++i)
                        {
                            const DT_ temp(((const SparseMatrix<DT_>)u)(k, i));
                            nzt = ((const SparseMatrix<DT_>)u)(pivot, i);
                            if (nzt != DT_(0) || temp != DT_(0))
                            {
                                nzt = u(k, i);
                                u(pivot, i) = temp;
                            }
                        }
                        for (unsigned long i(0) ; i < k  ; ++i)
                        {
                            const DT_ temp(((const SparseMatrix<DT_>)l)(k, i));
                            nzt = ((const SparseMatrix<DT_>)l)(pivot, i);
                            if (nzt != DT_(0) || temp != DT_(0))
                            {
                                l(k, i) = nzt;
                                l(pivot, i) = temp;
                            }
                        }
                        DT_ temp(b[k]);
                        b[k] = b[pivot];
                        b[pivot] = temp;
                    }

                    //todo calc and store LU insitu in A
                    //todo use non zero elemment iterators
                    for (unsigned long j(k + 1) ; j < u.rows() ; ++j)
                    {
                        nzt = ((const SparseMatrix<DT_>)u)(j, k) / ((const SparseMatrix<DT_>)u)(k, k);
                        if (nzt != DT_(0))
                            l(j, k) = nzt;
                        for (unsigned long i(k) ; i < u.rows() ; ++i)
                        {
                            nzt = ((const SparseMatrix<DT_>)u)(j, i) - ((const SparseMatrix<DT_>)l)(j, k) * ((const SparseMatrix<DT_>)u)(k, i);
                            if (nzt != DT_(0))
                                u(j, i) = nzt;
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
            }
    };
}
#endif
