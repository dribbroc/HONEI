/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007, 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/la/reduction.hh>
#include <honei/backends/sse/operations.hh>



using namespace honei;

float Reduction<rt_sum, tags::CPU::SSE>::value(const DenseVectorContinuousBase<float> & a)
{
    CONTEXT("When reducing DenseVectorContinuousBase<float> to sum (SSE):");

    return sse::reduction_sum(a.elements(), a.size());
}

double Reduction<rt_sum, tags::CPU::SSE>::value(const DenseVectorContinuousBase<double> & a)
{
    CONTEXT("When reducing DenseVectorContinuousBase<double> to sum (SSE):");

    return sse::reduction_sum(a.elements(), a.size());
}

DenseVector<float> Reduction<rt_sum, tags::CPU::SSE>::value(const DenseMatrix<float> & a)
{
    CONTEXT("When reducing DenseMatrix<float> to sum (SSE):");

    DenseVector<float> result(a.rows());

    for (unsigned long i(0) ; i < a.rows() ; ++i)
    {
        result[i] = sse::reduction_sum(a[i].elements(), a[i].size());
    }

    return result;
}

DenseVector<double> Reduction<rt_sum, tags::CPU::SSE>::value(const DenseMatrix<double> & a)
{
    CONTEXT("When reducing DenseMatrix<double> to sum (SSE):");

    DenseVector<double> result(a.rows());

    for (unsigned long i(0) ; i < a.rows() ; ++i)
    {
        result[i] = sse::reduction_sum(a[i].elements(), a[i].size());
    }

    return result;
}

float Reduction<rt_sum, tags::CPU::SSE>::value(const SparseVector<float> & a)
{
    CONTEXT("When reducing SparseVector<float> to sum (SSE):");

    return sse::reduction_sum(a.elements(), a.used_elements());
}

double Reduction<rt_sum, tags::CPU::SSE>::value(const SparseVector<double> & a)
{
    CONTEXT("When reducing SparseVector<double> to sum (SSE):");

    return sse::reduction_sum(a.elements(), a.used_elements());
}

DenseVector<float> Reduction<rt_sum, tags::CPU::SSE>::value(const SparseMatrix<float> & a)
{
    CONTEXT("When reducing SparseMatrix<float> to sum (SSE):");

    DenseVector<float> result(a.rows());

    for (unsigned long i(0) ; i < a.rows() ; ++i)
    {
        result[i] = sse::reduction_sum(a[i].elements(), a[i].used_elements());
    }

    return result;
}

DenseVector<double> Reduction<rt_sum, tags::CPU::SSE>::value(const SparseMatrix<double> & a)
{
    CONTEXT("When reducing SparseMatrix<double> to sum (SSE):");

    DenseVector<double> result(a.rows());

    for (unsigned long i(0) ; i < a.rows() ; ++i)
    {
        result[i] = sse::reduction_sum(a[i].elements(), a[i].used_elements());
    }

    return result;
}

