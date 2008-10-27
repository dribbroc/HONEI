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
#include<honei/la/algorithm.hh>
#include<honei/la/product.hh>
#include<honei/la/difference.hh>

using namespace honei;
namespace honei
{
    template<typename Tag_>
    struct Defect
    {
        public:
            template<typename DT_>
            static DenseVector<DT_> value(const DenseVector<DT_> & rhs, const BandedMatrixQ1<DT_> & system,
                    const DenseVector<DT_> & x)
            {
                if (x.size() != system.columns())
                {
                    throw VectorSizeDoesNotMatch(x.size(), system.columns());
                }
                if (rhs.size() != system.columns())
                {
                    throw VectorSizeDoesNotMatch(rhs.size(), system.columns());
                }

                DenseVector<DT_> result(rhs.size());
                copy<Tag_>(rhs, result);

                DenseVector<DT_> prod(Product<Tag_>::value(system, x));

                Difference<Tag_>::value(result, prod);

                return result;
            }
    };

    template<>
    struct Defect<tags::GPU::CUDA>
    {
        public:
            static DenseVector<float> value(const DenseVectorContinuousBase<float> & right_hand_side,
                    const BandedMatrixQ1<float> & system, const DenseVectorContinuousBase<float> & x);

    };
}

#endif
