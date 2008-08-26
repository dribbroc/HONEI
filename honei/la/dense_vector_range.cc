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

#include <honei/la/dense_vector_range.hh>
#include <honei/la/dense_vector_range-impl.hh>

namespace honei
{
    template class DenseVectorRange<float>;

    template bool operator== (const DenseVectorRange<float> & a, const DenseVectorRange<float> & b);

    template std::ostream & operator<< (std::ostream & lhs, const DenseVectorRange<float> & vector);

    template class DenseVectorRange<double>;

    template bool operator== (const DenseVectorRange<double> & a, const DenseVectorRange<double> & b);

    template std::ostream & operator<< (std::ostream & lhs, const DenseVectorRange<double> & vector);

    template class DenseVectorRange<int>;

    template bool operator== (const DenseVectorRange<int> & a, const DenseVectorRange<int> & b);

    template std::ostream & operator<< (std::ostream & lhs, const DenseVectorRange<int> & vector);

    template class DenseVectorRange<long>;

    template bool operator== (const DenseVectorRange<long> & a, const DenseVectorRange<long> & b);

    template std::ostream & operator<< (std::ostream & lhs, const DenseVectorRange<long> & vector);

    template class DenseVectorRange<unsigned long>;

    template bool operator== (const DenseVectorRange<unsigned long> & a, const DenseVectorRange<unsigned long> & b);

    template std::ostream & operator<< (std::ostream & lhs, const DenseVectorRange<unsigned long> & vector);
}

