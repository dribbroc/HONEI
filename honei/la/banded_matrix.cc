/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#include <honei/la/banded_matrix.hh>
#include <honei/la/banded_matrix-impl.hh>

namespace honei
{
    template <> const float Implementation<BandedMatrix<float> >::zero_element(0.0f);

    template class BandedMatrix<float>;

    template class BandIterator<type::Banded, float>;

    template class ConstBandIterator<type::Banded, float>;

    template class ConstElementIterator<storage::Banded, container::Matrix, float>;

    template class ElementIterator<storage::Banded, container::Matrix, float>;

    template bool operator== (const BandedMatrix<float> & a, const BandedMatrix<float> & b);

    template std::ostream & operator<< (std::ostream & lhs, const BandedMatrix<float> & matrix);

    template <> const double Implementation<BandedMatrix<double> >::zero_element(0.0);

    template class BandedMatrix<double>;

    template class BandIterator<type::Banded, double>;

    template class ConstBandIterator<type::Banded, double>;

    template class ConstElementIterator<storage::Banded, container::Matrix, double>;

    template class ElementIterator<storage::Banded, container::Matrix, double>;

    template bool operator== (const BandedMatrix<double> & a, const BandedMatrix<double> & b);

    template std::ostream & operator<< (std::ostream & lhs, const BandedMatrix<double> & matrix);
}

