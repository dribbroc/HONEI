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

#include <honei/la/dot_product.hh>
#include <honei/backends/sse/operations.hh>


using namespace honei;

float DotProduct<tags::CPU::SSE>::value(const DenseVectorContinuousBase<float> & a,
        const DenseVectorContinuousBase<float> & b)
{
    CONTEXT("When calculating dot-product of DenseVectorContinuousBase<float> with DenseVectorContinuousBase<float> "
            "(SSE):");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    a.lock(lm_read_only);
    b.lock(lm_read_only);
    a.unlock(lm_read_only);
    b.unlock(lm_read_only);
    return sse::dot_product(a.elements(), b.elements(), a.size());
}

double DotProduct<tags::CPU::SSE>::value(const DenseVectorContinuousBase<double> & a,
        const DenseVectorContinuousBase<double> & b)
{
    CONTEXT("When calculating dot-product of DenseVectorContinuousBase<double> with DenseVectorContinuousBase<double> "
            "(SSE):");

    if (a.size() != b.size())
        throw VectorSizeDoesNotMatch(b.size(), a.size());

    a.lock(lm_read_only);
    b.lock(lm_read_only);
    a.unlock(lm_read_only);
    b.unlock(lm_read_only);
    return sse::dot_product(a.elements(), b.elements(), a.size());
}

