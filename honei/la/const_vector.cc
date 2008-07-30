/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/la/const_vector.hh>
#include <honei/la/const_vector-impl.hh>

namespace honei
{
    template class ConstElementIterator<storage::Const, container::Vector, float>;

    template class ConstVector<float>;

    template bool operator== (const ConstVector<float> & a, const ConstVector<float> & b);

    template std::ostream & operator<< (std::ostream & lhs, const ConstVector<float> & vector);

    template class ConstElementIterator<storage::Const, container::Vector, double>;

    template class ConstVector<double>;

    template bool operator== (const ConstVector<double> & a, const ConstVector<double> & b);

    template std::ostream & operator<< (std::ostream & lhs, const ConstVector<double> & vector);

    template class ConstElementIterator<storage::Const, container::Vector, int>;

    template class ConstVector<int>;

    template bool operator== (const ConstVector<int> & a, const ConstVector<int> & b);

    template std::ostream & operator<< (std::ostream & lhs, const ConstVector<int> & vector);

    template class ConstElementIterator<storage::Const, container::Vector, long>;

    template class ConstVector<long>;

    template bool operator== (const ConstVector<long> & a, const ConstVector<long> & b);

    template std::ostream & operator<< (std::ostream & lhs, const ConstVector<long> & vector);

    template class ConstElementIterator<storage::Const, container::Vector, unsigned long>;

    template class ConstVector<unsigned long>;

    template bool operator== (const ConstVector<unsigned long> & a, const ConstVector<unsigned long> & b);

    template std::ostream & operator<< (std::ostream & lhs, const ConstVector<unsigned long> & vector);
}

