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
    template class ConstVector<float>;

    template bool operator== (const ConstVector<float> & a, const ConstVector<float> & b);

    template std::ostream & operator<< (std::ostream & lhs, const ConstVector<float> & vector);

    template class ConstVector<double>;

    template bool operator== (const ConstVector<double> & a, const ConstVector<double> & b);

    template std::ostream & operator<< (std::ostream & lhs, const ConstVector<double> & vector);
}

