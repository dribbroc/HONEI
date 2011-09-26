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

#include <honei/la/banded_matrix_qx.hh>
#include <honei/la/banded_matrix_qx-impl.hh>

namespace honei
{
    template <> const float Implementation<BandedMatrixQx<Q1Type, float> >::zero_element(0.0f);

    template class BandedMatrixQx<Q1Type, float>;

    template bool operator== (const BandedMatrixQx<Q1Type, float> & a, const BandedMatrixQx<Q1Type, float> & b);

    template std::ostream & operator<< (std::ostream & lhs, const BandedMatrixQx<Q1Type, float> & matrix);

    template <> const double Implementation<BandedMatrixQx<Q1Type, double> >::zero_element(0.0);

    template class BandedMatrixQx<Q1Type, double>;

    template bool operator== (const BandedMatrixQx<Q1Type, double> & a, const BandedMatrixQx<Q1Type, double> & b);

    template std::ostream & operator<< (std::ostream & lhs, const BandedMatrixQx<Q1Type, double> & matrix);
}

