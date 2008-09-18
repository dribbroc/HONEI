/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the Math C++ library. Math is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * Math is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */


#ifndef LIBMATH_GUARD_DEFECT_HH
#define LIBMATH_GUARD_DEFECT_HH 1

#include<honei/la/banded_matrix_q1.hh>
#include<honei/la/dense_vector.hh>

using namespace honei;
namespace honei
{
    template<typename Tag_>
    struct Defect
    {
        public:
            template<typename DT_>
            static DenseVector<DT_>  value(DenseVector<DT_> & rhs, BandedMatrixQ1<DT_> & system, DenseVector<DT_> & x)
            {
                DenseVector<DT_> result(rhs.copy());

                system.lock(lm_read_only);
                x.lock(lm_read_only);
                DenseVector<DT_> prod(Product<Tag_>::value(system, x));
                system.unlock(lm_read_only);
                x.unlock(lm_read_only);

                Difference<Tag_>::value(result, prod);

                return result;
            }
    };
}

#endif
