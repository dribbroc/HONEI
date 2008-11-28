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

#include <honei/la/element_inverse.hh>
#include <honei/backends/sse/operations.hh>


using namespace honei;

DenseVectorContinuousBase<float> & ElementInverse<tags::CPU::SSE>::value(DenseVectorContinuousBase<float> & x)
{
    CONTEXT("When inverting DenseVectorContinuousBase<float> (SSE):");

    sse::element_inverse(x.elements(), x.size());

    return x;
}

DenseVectorContinuousBase<double> & ElementInverse<tags::CPU::SSE>::value(DenseVectorContinuousBase<double> & x)
{
    CONTEXT("When inverting DenseVectorContinuousBase<double> (SSE):");

    sse::element_inverse(x.elements(), x.size());

    return x;
}

DenseMatrix<float> & ElementInverse<tags::CPU::SSE>::value(DenseMatrix<float> & x)
{
    CONTEXT("When inverting DenseMatrix<float> (SSE):");

    sse::element_inverse(x.elements(), x.rows() * x.columns());

    return x;
}

DenseMatrix<double> & ElementInverse<tags::CPU::SSE>::value(DenseMatrix<double> & x)
{
    CONTEXT("When inverting DenseMatrix<double> (SSE):");

    sse::element_inverse(x.elements(), x.rows() * x.columns());

    return x;
}

SparseVector<float> & ElementInverse<tags::CPU::SSE>::value(SparseVector<float> & x)
{
    CONTEXT("When inverting SparseVector<float> (SSE):");

    sse::element_inverse(x.elements(), x.used_elements());

    return x;
}

SparseVector<double> & ElementInverse<tags::CPU::SSE>::value(SparseVector<double> & x)
{
    CONTEXT("When inverting SparseVector<double> (SSE):");

    sse::element_inverse(x.elements(), x.used_elements());

    return x;
}

/*SparseMatrix<float> & ElementInverse<tags::CPU::SSE>::value(SparseMatrix<float> & x)
{
    CONTEXT("When invertingSparseMatrix<float> (SSE):");

    for (SparseMatrix<float>::RowIterator l(x.begin_non_zero_rows()),
            l_end(x.end_non_zero_rows()) ; l != l_end ; ++l)
    {
        sse::element_inverse(l->elements(), l->used_elements());
    }

    return x;
}

SparseMatrix<double> & ElementInverse<tags::CPU::SSE>::value(SparseMatrix<double> & x)
{
    CONTEXT("When inverting SparseMatrix<double> (SSE):");

    for (SparseMatrix<double>::RowIterator l(x.begin_non_zero_rows()),
            l_end(x.end_non_zero_rows()) ; l != l_end ; ++l)
    {
        sse::element_inverse(l->elements(), l->used_elements());
    }

    return x;
}*/

