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

#include <honei/la/scale.hh>
#include <honei/backends/sse/operations.hh>


using namespace honei;

DenseVectorContinuousBase<float> & Scale<tags::CPU::SSE>::value(DenseVectorContinuousBase<float> & x, const float a)
{
    CONTEXT("When scaling DenseVectorContinuousBase<float> by float (SSE):");

    sse::scale(a, x.elements(), x.size());

    return x;
}

DenseVectorContinuousBase<double> & Scale<tags::CPU::SSE>::value(DenseVectorContinuousBase<double> & x, const double a)
{
    CONTEXT("When scaling DenseVectorContinuousBase<double> by double (SSE):");

    sse::scale(a, x.elements(), x.size());

    return x;
}

DenseMatrix<float> & Scale<tags::CPU::SSE>::value(DenseMatrix<float> & x, const float a)
{
    CONTEXT("When scaling DenseMatrix<float> by float (SSE):");

    sse::scale(a, x.elements(), x.rows() * x.columns());

    return x;
}

DenseMatrix<double> & Scale<tags::CPU::SSE>::value(DenseMatrix<double> & x, const double a)
{
    CONTEXT("When scaling DenseMatrix<double> by double (SSE):");

    sse::scale(a, x.elements(), x.rows() * x.columns());

    return x;
}

SparseVector<float> & Scale<tags::CPU::SSE>::value(SparseVector<float> & x, const float a)
{
    CONTEXT("When scaling SparseVector<float> by float (SSE):");

    sse::scale(a, x.elements(), x.used_elements());

    return x;
}

SparseVector<double> & Scale<tags::CPU::SSE>::value(SparseVector<double> & x, const double a)
{
    CONTEXT("When scaling SparseVector<double> by double (SSE):");

    sse::scale(a, x.elements(), x.used_elements());

    return x;
}

/*SparseMatrix<float> & Scale<tags::CPU::SSE>::value(SparseMatrix<float> & x, const float a)
{
    CONTEXT("When scaling SparseMatrix<float> by float (SSE):");

    for (SparseMatrix<float>::RowIterator l(x.begin_non_zero_rows()),
            l_end(x.end_non_zero_rows()) ; l != l_end ; ++l)
    {
        sse::scale(a, l->elements(), l->used_elements());
    }

    return x;
}

SparseMatrix<double> & Scale<tags::CPU::SSE>::value(SparseMatrix<double> & x, const double a)
{
    CONTEXT("When scaling SparseMatrix<double> by double (SSE):");

    for (SparseMatrix<double>::RowIterator l(x.begin_non_zero_rows()),
            l_end(x.end_non_zero_rows()) ; l != l_end ; ++l)
    {
        sse::scale(a, l->elements(), l->used_elements());
    }

    return x;
}*/

