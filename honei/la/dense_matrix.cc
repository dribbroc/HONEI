/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@uni-dortmund.de>
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

#include <honei/la/dense_matrix.hh>
#include <honei/la/dense_matrix-impl.hh>

namespace honei
{
    template class DenseMatrix<float>;

    template class ConstElementIterator<storage::Dense, container::Matrix, float>;

    template class ElementIterator<storage::Dense, container::Matrix, float>;

    template bool operator== (const DenseMatrix<float> & a, const DenseMatrix<float> & b);

    template std::ostream & operator<< (std::ostream & lhs, const DenseMatrix<float> & matrix);

    template class DenseMatrix<double>;

    template class ConstElementIterator<storage::Dense, container::Matrix, double>;

    template class ElementIterator<storage::Dense, container::Matrix, double>;

    template bool operator== (const DenseMatrix<double> & a, const DenseMatrix<double> & b);

    template std::ostream & operator<< (std::ostream & lhs, const DenseMatrix<double> & matrix);

    template class DenseMatrix<long>;

    template class ConstElementIterator<storage::Dense, container::Matrix, long>;

    template class ElementIterator<storage::Dense, container::Matrix, long>;

    template bool operator== (const DenseMatrix<long> & a, const DenseMatrix<long> & b);

    template std::ostream & operator<< (std::ostream & lhs, const DenseMatrix<long> & matrix);

    template class DenseMatrix<bool>;

    template class ConstElementIterator<storage::Dense, container::Matrix, bool>;

    template class ElementIterator<storage::Dense, container::Matrix, bool>;

    template bool operator== (const DenseMatrix<bool> & a, const DenseMatrix<bool> & b);

    template std::ostream & operator<< (std::ostream & lhs, const DenseMatrix<bool> & matrix);
}

