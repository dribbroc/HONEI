/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007, 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#include <honei/la/absolute.hh>
#include <cmath>

using namespace honei;

template <>
DenseVectorBase<float> &
Absolute<tags::CPU>::value(DenseVectorBase<float> & x)
{
    CONTEXT("When calculating the absolute value of DenseVectorBase<float> elements:");

    for (DenseVectorBase<float>::ElementIterator i(x.begin_elements()), i_end(x.end_elements()) ;
            i != i_end ; ++i)
    {
        *i = fabs(*i);
    }

    return x;
}

template <>
SparseVector<float> &
Absolute<tags::CPU>::value(SparseVector<float> & x)
{
    CONTEXT("When calculating the absolute value of SparseVector<float> elements:");

    for (SparseVector<float>::NonZeroElementIterator i(x.begin_non_zero_elements()), i_end(x.end_non_zero_elements()) ;
            i != i_end ; ++i)
    {
        *i = fabs(*i);
    }

    return x;
}

template <>
DenseVectorBase<double> &
Absolute<tags::CPU>::value(DenseVectorBase<double> & x)
{
    CONTEXT("When calculating the absolute value of DenseVectorBase<double> elements:");

    for (DenseVectorBase<double>::ElementIterator i(x.begin_elements()), i_end(x.end_elements()) ;
            i != i_end ; ++i)
    {
        *i = fabs(*i);
    }

    return x;
}

template <>
SparseVector<double> &
Absolute<tags::CPU>::value(SparseVector<double> & x)
{
    CONTEXT("When calculating the absolute value of SparseVector<double> elements:");

    for (SparseVector<double>::NonZeroElementIterator i(x.begin_non_zero_elements()), i_end(x.end_non_zero_elements()) ;
            i != i_end ; ++i)
    {
        *i = fabs(*i);
    }

    return x;
}

