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

#include <honei/la/scaled_sum.hh>
#include <honei/backends/sse/operations.hh>


using namespace honei;

DenseVectorContinuousBase<float> & ScaledSum<tags::CPU::SSE>::value(DenseVectorContinuousBase<float> & x,
        const DenseVectorContinuousBase<float> & y, float b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<float> (SSE):");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());

    sse::scaled_sum(x.elements(), y.elements(), b, x.size());

    return x;
}

DenseVectorContinuousBase<double> & ScaledSum<tags::CPU::SSE>::value(DenseVectorContinuousBase<double> & x,
        const DenseVectorContinuousBase<double> & y, double b)
{
    CONTEXT("When calculating ScaledSum form DenseVectoriContinuousBase<double> (SSE):");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());

    sse::scaled_sum(x.elements(), y.elements(), b, x.size());

    return x;
}

DenseVectorContinuousBase<float> & ScaledSum<tags::CPU::SSE>::value(DenseVectorContinuousBase<float> & a,
        const DenseVectorContinuousBase<float> & b, const DenseVectorContinuousBase<float> & c)
{
    CONTEXT("When calculating ScaledSum (DenseVectorContinuousBase<float>, DenseVectorContinuousBase<float>, "
            "DenseVectorContinuousBase<float>) (SSE):");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    if (a.size() != c.size())
        throw VectorSizeDoesNotMatch(c.size(), a.size());

    sse::scaled_sum(a.elements(), b.elements(), c.elements(), a.size());

    return a;
}

DenseVectorContinuousBase<double> & ScaledSum<tags::CPU::SSE>::value(DenseVectorContinuousBase<double> & a,
        const DenseVectorContinuousBase<double> & b, const DenseVectorContinuousBase<double> & c)
{
    CONTEXT("When calculating ScaledSum (DenseVectorContinuousBase<doule>, DenseVectorContinuousBase<double>, "
            "DenseVectorContinuousBase<double>) (SSE):");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    if (a.size() != c.size())
        throw VectorSizeDoesNotMatch(c.size(), a.size());

    sse::scaled_sum(a.elements(), b.elements(), c.elements(), a.size());

    return a;
}

DenseVectorContinuousBase<float> & ScaledSum<tags::CPU::SSE>::value(DenseVectorContinuousBase<float> & x,
        const DenseVectorContinuousBase<float> & y, const DenseVectorContinuousBase<float> & z, float b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<float> (SSE):");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());
    if (x.size() != z.size())
        throw VectorSizeDoesNotMatch(x.size(), z.size());

    sse::scaled_sum(x.elements(), y.elements(), z.elements(), b, x.size());

    return x;
}

DenseVectorContinuousBase<double> & ScaledSum<tags::CPU::SSE>::value(DenseVectorContinuousBase<double> & x,
        const DenseVectorContinuousBase<double> & y, const DenseVectorContinuousBase<double> & z, double b)
{
    CONTEXT("When calculating ScaledSum form DenseVectorContinuousBase<double> (SSE):");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());
    if (x.size() != z.size())
        throw VectorSizeDoesNotMatch(x.size(), z.size());

    sse::scaled_sum(x.elements(), y.elements(), z.elements(), b, x.size());

    return x;
}
