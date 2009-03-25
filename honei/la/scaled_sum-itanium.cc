/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2009 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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
#include <honei/backends/itanium/operations.hh>


using namespace honei;

template <typename DT1_, typename DT2_>
DenseVectorContinuousBase<DT1_> & ScaledSum<tags::CPU::Itanium>::value(DenseVectorContinuousBase<DT1_> & x,
        const DenseVectorContinuousBase<DT2_> & y, DT2_ b)
{
    CONTEXT("When calculating ScaledSum from DenseVectorContinuousBase (Itanium):");

    if (x.size() != y.size())
        throw VectorSizeDoesNotMatch(x.size(), y.size());

    itanium::scaled_sum(x.elements(), y.elements(), b, x.size());

    return x;
}

template DenseVectorContinuousBase<float> & ScaledSum<tags::CPU::Itanium>::value<float, float>(DenseVectorContinuousBase<float> &, const DenseVectorContinuousBase<float> &, float);
template DenseVectorContinuousBase<double> & ScaledSum<tags::CPU::Itanium>::value<double, double>(DenseVectorContinuousBase<double> &, const DenseVectorContinuousBase<double> &, double);
