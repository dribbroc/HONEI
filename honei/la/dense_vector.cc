/* vim: set sw=4 sts=4 et nofoldenable : */

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

#include <honei/la/dense_vector.hh>
#include <honei/la/dense_vector-impl.hh>

namespace honei
{
    template class ConstElementIterator<storage::Dense, container::Vector, float>;

    template class DenseVector<float>;

    template class ElementIterator<storage::Dense, container::Vector, float>;

    template bool operator== (const DenseVectorBase<float> & a, const DenseVectorBase<float> & b);

    template std::ostream & operator<< (std::ostream & lhs, const DenseVector<float> & vector);

    template class ConstElementIterator<storage::Dense, container::Vector, double>;

    template class DenseVector<double>;

    template class ElementIterator<storage::Dense, container::Vector, double>;

    template bool operator== (const DenseVectorBase<double> & a, const DenseVectorBase<double> & b);

    template std::ostream & operator<< (std::ostream & lhs, const DenseVector<double> & vector);

    template class ConstElementIterator<storage::Dense, container::Vector, int>;

    template class DenseVector<int>;

    template class ElementIterator<storage::Dense, container::Vector, int>;

    template bool operator== (const DenseVectorBase<int> & a, const DenseVectorBase<int> & b);

    template std::ostream & operator<< (std::ostream & lhs, const DenseVector<int> & vector);

    template class ConstElementIterator<storage::Dense, container::Vector, long>;

    template class DenseVector<long>;

    template class ElementIterator<storage::Dense, container::Vector, long>;

    template bool operator== (const DenseVectorBase<long> & a, const DenseVectorBase<long> & b);

    template std::ostream & operator<< (std::ostream & lhs, const DenseVector<long> & vector);

    template class ConstElementIterator<storage::Dense, container::Vector, unsigned long>;

    template class DenseVector<unsigned long>;

    template class ElementIterator<storage::Dense, container::Vector, unsigned long>;

    template bool operator== (const DenseVectorBase<unsigned long> & a, const DenseVectorBase<unsigned long> & b);

    template std::ostream & operator<< (std::ostream & lhs, const DenseVector<unsigned long> & vector);
}

