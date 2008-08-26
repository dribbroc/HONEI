/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Volker Jung <volker.jung@uni-dortmund.de>
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

#include <honei/la/dense_vector_slice.hh>
#include <honei/la/dense_vector_slice-impl.hh>

namespace honei
{
    template class DenseVectorSlice<float>;

    template bool operator== (const DenseVectorSlice<float> & a, const DenseVectorSlice<float> & b);

    template std::ostream & operator<< (std::ostream & lhs, const DenseVectorSlice<float> & vector);

    template class DenseVectorSlice<double>;

    template bool operator== (const DenseVectorSlice<double> & a, const DenseVectorSlice<double> & b);

    template std::ostream & operator<< (std::ostream & lhs, const DenseVectorSlice<double> & vector);

    template class DenseVectorSlice<int>;

    template bool operator== (const DenseVectorSlice<int> & a, const DenseVectorSlice<int> & b);

    template std::ostream & operator<< (std::ostream & lhs, const DenseVectorSlice<int> & vector);

    template class DenseVectorSlice<long>;

    template bool operator== (const DenseVectorSlice<long> & a, const DenseVectorSlice<long> & b);

    template std::ostream & operator<< (std::ostream & lhs, const DenseVectorSlice<long> & vector);

    template class DenseVectorSlice<unsigned long>;

    template bool operator== (const DenseVectorSlice<unsigned long> & a, const DenseVectorSlice<unsigned long> & b);

    template std::ostream & operator<< (std::ostream & lhs, const DenseVectorSlice<unsigned long> & vector);
}

