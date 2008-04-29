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

#include <honei/la/element_product.hh>
#include <honei/backends/sse/operations.hh>


using namespace honei;

DenseVectorContinuousBase<float> & ElementProduct<tags::CPU::SSE>::value(DenseVectorContinuousBase<float> & a,
        const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When multiplying DenseVectorContinuousBase<float> and DenseVectorContinuousBase<float> elementwise "
            "(SSE):");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    sse::element_product(a.elements(), b.elements(), a.size());

    return a;
}

DenseVectorContinuousBase<double> & ElementProduct<tags::CPU::SSE>::value(DenseVectorContinuousBase<double> & a,
        const DenseVectorContinuousBase<double> & b)
{
    CONTEXT("When multiplying DenseVectorContinuousBase<double> and DenseVectorContinuousBase<double> elementwise "
            "(SSE):");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    sse::element_product(a.elements(), b.elements(), a.size());

    return a;
}

DenseMatrix<float> & ElementProduct<tags::CPU::SSE>::value(DenseMatrix<float> & a, const DenseMatrix<float> & b)
{
    CONTEXT("When multiplying DenseMatrix<float> and DenseMatrix<float> elementwise (SSE):");

    if (a.columns() != b.columns())
    {
        throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
    }

    if (a.rows() != b.rows())
    {
        throw MatrixRowsDoNotMatch(b.rows(), a.rows());
    }

    sse::element_product(a.elements(), b.elements(), a.rows() * a.columns());

    return a;
}

DenseMatrix<double> & ElementProduct<tags::CPU::SSE>::value(DenseMatrix<double> & a, const DenseMatrix<double> & b)
{
    CONTEXT("When multiplying DenseMatrix<double> and DenseMatrix<double> elementwise (SSE):");

    if (a.columns() != b.columns())
    {
        throw MatrixColumnsDoNotMatch(b.columns(), a.columns());
    }

    if (a.rows() != b.rows())
    {
        throw MatrixRowsDoNotMatch(b.rows(), a.rows());
    }

    sse::element_product(a.elements(), b.elements(), a.rows() * a.columns());

    return a;
}

