/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2009 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the MATH C++ library. Math is free software;
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

#ifndef LIBMATH_GUARD_SCALED_PRODUCT_SUM_NORM_HH
#define LIBMATH_GUARD_SCALED_PRODUCT_SUM_NORM_HH 1

#include <honei/util/tags.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/product.hh>
#include <honei/la/scaled_sum.hh>
#include <honei/la/norm.hh>
#include <honei/la/element_inverse.hh>
#include <honei/util/memory_arbiter.hh>
#include <honei/la/algorithm.hh>

namespace honei
{
    //Tutorial version
    template<typename Tag_>
        struct ScaledProductSumNorm_TUTORIAL
        {
            template<typename DataType_>
                static inline DataType_ value(DataType_ alpha,
                                              DenseVector<DataType_> & y,
                                              DataType_ beta,
                                              BandedMatrixQ1<DataType_> & A,
                                              DenseVector<DataType_> & x)
                {
                    y.lock(lm_read_only);
                    A.lock(lm_read_only);
                    x.lock(lm_read_only);

                    DenseVector<DataType_> result(y.copy());
                    Scale<Tag_>::value(result, alpha);
                    result = ScaledSum<Tag_>::value(result, Product<Tag_>::value(A , x), beta);

                    y.unlock(lm_read_only);
                    A.unlock(lm_read_only);
                    x.unlock(lm_read_only);

                    return Norm<vnt_l_two, false, Tag_>::value(result);

                }
        };

    template<>
        struct ScaledProductSumNorm_TUTORIAL<tags::CPU::SSE>
        {
            static float value(float a, DenseVector<float> & y, float b, BandedMatrixQ1<float> & A, DenseVector<float> & x);
            static double value(double a, DenseVector<double> & y, double b, BandedMatrixQ1<double> & A, DenseVector<double> & x);
        };


    //Fully optimized version
    template<typename Tag_>
        struct ScaledProductSumNorm
        {
            template<typename DataType_>
                static inline DataType_ value(DataType_ alpha,
                                              DenseVector<DataType_> & y,
                                              DataType_ beta,
                                              BandedMatrixQ1<DataType_> & A,
                                              DenseVector<DataType_> & x)
                {
                    return ScaledProductSumNorm_TUTORIAL<Tag_>::value(alpha, y, beta, A, x);
                }
        };

    template<>
        struct ScaledProductSumNorm<tags::CPU::SSE>
        {
            static float value(float a, DenseVector<float> & y, float b, BandedMatrixQ1<float> & A, DenseVector<float> & x);
            static double value(double a, DenseVector<double> & y, double b, BandedMatrixQ1<double> & A, DenseVector<double> & x);
        };

}

#endif
