/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#include <honei/la/sparse_vector.hh>
#include <honei/la/sparse_vector-impl.hh>

namespace honei
{
    template class ConstElementIterator<storage::Sparse, container::Vector, float>;

    template class ConstElementIterator<storage::SparseNonZero, container::Vector, float>;

    template class SparseVector<float>;

    template class ElementIterator<storage::Sparse, container::Vector, float>;

    template class ElementIterator<storage::SparseNonZero, container::Vector, float>;

    template bool operator== (const SparseVector<float> & a, const SparseVector<float> & b);

    template std::ostream & operator<< (std::ostream & lhs, const SparseVector<float> & vector);

    template class ConstElementIterator<storage::Sparse, container::Vector, double>;

    template class ConstElementIterator<storage::SparseNonZero, container::Vector, double>;

    template class SparseVector<double>;

    template class ElementIterator<storage::Sparse, container::Vector, double>;

    template class ElementIterator<storage::SparseNonZero, container::Vector, double>;

    template bool operator== (const SparseVector<double> & a, const SparseVector<double> & b);

    template std::ostream & operator<< (std::ostream & lhs, const SparseVector<double> & vector);

    template class ConstElementIterator<storage::Sparse, container::Vector, int>;

    template class ConstElementIterator<storage::SparseNonZero, container::Vector, int>;

    template class SparseVector<int>;

    template class ElementIterator<storage::Sparse, container::Vector, int>;

    template class ElementIterator<storage::SparseNonZero, container::Vector, int>;

    template bool operator== (const SparseVector<int> & a, const SparseVector<int> & b);

    template std::ostream & operator<< (std::ostream & lhs, const SparseVector<int> & vector);

    template class ConstElementIterator<storage::Sparse, container::Vector, long>;

    template class ConstElementIterator<storage::SparseNonZero, container::Vector, long>;

    template class SparseVector<long>;

    template class ElementIterator<storage::Sparse, container::Vector, long>;

    template class ElementIterator<storage::SparseNonZero, container::Vector, long>;

    template bool operator== (const SparseVector<long> & a, const SparseVector<long> & b);

    template std::ostream & operator<< (std::ostream & lhs, const SparseVector<long> & vector);

    template class ConstElementIterator<storage::Sparse, container::Vector, unsigned long>;

    template class ConstElementIterator<storage::SparseNonZero, container::Vector, unsigned long>;

    template class SparseVector<unsigned long>;

    template class ElementIterator<storage::Sparse, container::Vector, unsigned long>;

    template class ElementIterator<storage::SparseNonZero, container::Vector, unsigned long>;

    template bool operator== (const SparseVector<unsigned long> & a, const SparseVector<unsigned long> & b);

    template std::ostream & operator<< (std::ostream & lhs, const SparseVector<unsigned long> & vector);
}

