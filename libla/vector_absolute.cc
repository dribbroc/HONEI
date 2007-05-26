/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#include <libla/vector_absolute.hh>

using namespace pg512;

DenseVector<float> & VectorAbsolute<float, tags::CPU>::value(DenseVector<float> & vector)
{
    for (Vector<float>::ElementIterator i(vector.begin_elements()), i_end(vector.end_elements()) ;
            i != i_end ; ++i)
    {
        *i = fabs(*i);
    }

    return vector;
}

SparseVector<float> & VectorAbsolute<float, tags::CPU>::value(SparseVector<float> & vector)
{
    for (Vector<float>::ElementIterator i(vector.begin_non_zero_elements()), i_end(vector.end_non_zero_elements()) ;
            i != i_end ; ++i)
    {
        *i = fabs(*i);
    }

    return vector;
}

DenseVector<double> & VectorAbsolute<double, tags::CPU>::value(DenseVector<double> & vector)
{
    for (Vector<double>::ElementIterator i(vector.begin_elements()), i_end(vector.end_elements()) ;
            i != i_end ; ++i)
    {
        *i = fabs(*i);
    }

    return vector;
}

SparseVector<double> & VectorAbsolute<double, tags::CPU>::value(SparseVector<double> & vector)
{
    for (Vector<double>::ElementIterator i(vector.begin_non_zero_elements()), i_end(vector.end_non_zero_elements()) ;
            i != i_end ; ++i)
    {
        *i = fabs(*i);
    }

    return vector;
}

