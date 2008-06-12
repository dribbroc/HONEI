/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the HONEI C++ library. HOENI is free software;
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

#include <honei/la/banded_matrix_q1.hh>
#include <honei/la/banded_matrix_q1-impl.hh>

namespace honei
{
    template <> const float Implementation<BandedMatrixQ1<float> >::zero_element(0.0f);

    template class BandedMatrixQ1<float>;

    template bool operator== (const BandedMatrixQ1<float> & a, const BandedMatrixQ1<float> & b);

    template std::ostream & operator<< (std::ostream & lhs, const BandedMatrixQ1<float> & matrix);

    template <> const double Implementation<BandedMatrixQ1<double> >::zero_element(0.0);

    template class BandedMatrixQ1<double>;

    template bool operator== (const BandedMatrixQ1<double> & a, const BandedMatrixQ1<double> & b);

    template std::ostream & operator<< (std::ostream & lhs, const BandedMatrixQ1<double> & matrix);
}

