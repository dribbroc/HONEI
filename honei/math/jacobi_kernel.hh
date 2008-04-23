/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
* Copyright (c) 2007 , 2008 Markus Geveler <apryde@gmx.de>
*
* This file is part of the MATH C++ library. LibMath is free software;
* you can redistribute it and/or modify it under the terms of the GNU General
* Public License version 2, as published by the Free Software Foundation.
*
* LibMath is distributed in the hope that it will be useful, but WITHOUT ANY
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
* FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
* details.
*
* You should have received a copy of the GNU General Public License along with
* this program; if not, write to the Free Software Foundation, Inc., 59 Temple
* Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef LIBMATH_GUARD_JACOBI_KERNEL_HH
#define LIBMATH_GUARD_JACOBI_KERNEL_HH 1

#include <honei/util/tags.hh>
#include <honei/la/dense_matrix.hh>
#include <honei/la/dense_vector.hh>
#include <honei/la/product.hh>
#include <honei/la/sum.hh>
#include <honei/la/difference.hh>
#include <honei/la/dot_product.hh>
#include <honei/la/scale.hh>
#include <honei/la/norm.hh>
#include <iostream>

namespace honei
{
    /**
     * \brief Solution of LES with Jacobi method.
     *
     * The referenced containers are invariant under this operation.
     * In every case, a new object is created and returned.
     *
     * \ingroup grpmatrixoperations
     * \ingroup grpvectoroperations
     */
    template <typename Tag_=tags::CPU>
    struct JacobiKernel
    {
    };

    template<>
    struct JacobiKernel<tags::CPU>
    {
        public:
            template<typename DT1_, typename DT2_>
            static inline  DenseVector<DT1_> value(DenseVector<DT2_> & right_hand_side, DenseVector<DT1_> & former_result, DenseVector<DT1_> & diag_inverted, BandedMatrix<DT1_> & difference)
            {
                DenseVector<DT1_> temp = Product<tags::CPU>::value(difference, former_result.copy());

                DenseVector<DT1_> temp2(right_hand_side.copy());

                Difference<tags::CPU>::value(temp2, temp);
                ElementProduct<tags::CPU>::value(temp2, diag_inverted);
                former_result = temp2;

                return former_result;
            }
    };

    template<>
    struct JacobiKernel<tags::CPU::SSE>
    {
        public:
            static /*inline*/ DenseVector<float> value(DenseVector<float> & b, DenseVector<float> & x, DenseVector<float> & d, BandedMatrix<float> & a);

            static /*inline*/ DenseVector<double> value(DenseVector<double> & b, DenseVector<double> & x, DenseVector<double> & d, BandedMatrix<double> & a);

    };

    template<>
    struct JacobiKernel<tags::Cell>
    {
        public:
            template<typename DT1_, typename DT2_>
            static inline  DenseVector<DT1_> value(DenseVector<DT2_> & right_hand_side, DenseVector<DT1_> & former_result, DenseVector<DT1_> & diag_inverted, BandedMatrix<DT1_> & difference);
    };
}
#endif
