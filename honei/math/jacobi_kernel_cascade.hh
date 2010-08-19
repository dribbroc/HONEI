/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
* Copyright (c) 2007 , 2008 Markus Geveler <apryde@gmx.de>
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

#ifndef MATH_GUARD_JACOBI_KERNEL_CASCADE_HH
#define MATH_GUARD_JACOBI_KERNEL_CASCADE_HH 1

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
#include <cmath>

namespace honei
{
    namespace cascade_quantities
    {
        class ONCE;
        class TWICE;
    }

    namespace cascade_specializations
    {
        class NONE;
        class Q1;
    }

    /**
     * \brief Solution of linear system with Jacobi method, coupling several
     * steps.
     *
     * The referenced containers are invariant under this operation.
     * In every case, a new object is created and returned.
     *
     * \ingroup grpmatrixoperations
     * \ingroup grpvectoroperations
     */

    using namespace cascade_quantities;
    using namespace cascade_specializations;

    template <typename Tag_=tags::CPU, typename Quantity_=ONCE, typename Spec_=NONE>
    struct JacobiKernelCascade
    {
    };

    template<>
    struct JacobiKernelCascade<tags::CPU, ONCE>
    {
        public:
            template<typename DT1_, typename DT2_>
            static inline  DenseVector<DT1_> value(DenseVector<DT2_> & right_hand_side, DenseVector<DT1_> & former_result, DenseVector<DT1_> & diag_inverted, BandedMatrix<DT1_> & difference, HONEI_UNUSED DenseVector<DT2_> & scaled_d)
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
    struct JacobiKernelCascade<tags::CPU, ONCE, Q1>
    {
        public:
            template<typename DT1_, typename DT2_>
            static inline  DenseVector<DT1_> value(DenseVector<DT2_> & right_hand_side, DenseVector<DT1_> & former_result, DenseVector<DT1_> & diag_inverted, BandedMatrix<DT1_> & difference, DenseVector<DT2_> & scaled_d)
            {
                unsigned long n = right_hand_side.size();
                unsigned long root_n = (unsigned long)sqrt(n);
                DenseVector<DT1_> result(n);

                DT1_ * x_old = former_result.elements();
                DT1_ * d_i = diag_inverted.elements();
                DT2_ * d_i_s = scaled_d.elements();
                DT1_ * x_new = result.elements();

                DT1_ * ll = difference.band((signed long)(-root_n - 1)).elements();
                DT1_ * ld = difference.band((signed long)(-root_n)).elements();
                DT1_ * lu = difference.band((signed long)(-root_n + 1)).elements();

                DT1_ * dl = difference.band((signed long)(-1)).elements();
                //DT1_ * dd = difference.band((signed long)(0)).elements();
                DT1_ * du = difference.band((signed long)(1)).elements();

                DT1_ * ul = difference.band((signed long)(root_n - 1)).elements();
                DT1_ * ud = difference.band((signed long)(root_n)).elements();
                DT1_ * uu = difference.band((signed long)(root_n - 1)).elements();

                unsigned long i(0);
                //index 0
                x_new[i] = //dd[i] * x_old[i] +
                           du[i] * x_old[1] +
                           ul[i] * x_old[root_n - 1] +
                           ud[i] * x_old[root_n] +
                           uu[i] * x_old[root_n + 1];
                x_new[i] *= d_i[i];
                x_new[i] += d_i_s[i];

                //index in [1, root_n -1[
                i = 1;
                for(; i < root_n - 1 ; ++i)
                {
                    x_new[i] = dl[i] * x_old[i - 1] +
                               //dd[i] * x_old[i] +
                               du[i] * x_old[i + 1] +
                               ul[i] * x_old[i + root_n - 1] +
                               ud[i] * x_old[i + root_n] +
                               uu[i] * x_old[i + root_n + 1];
                    x_new[i] *= d_i[i];
                    x_new[i] += d_i_s[i];
                }

                //index root_n -1
                i = root_n - 1;
                x_new[i] = lu[i] * x_old[i - (root_n - 1)] +
                           dl[i] * x_old[i - 1] +
                           //dd[i] * x_old[i] +
                           du[i] * x_old[i + 1] +
                           ul[i] * x_old[i + root_n - 1] +
                           ud[i] * x_old[i + root_n] +
                           uu[i] * x_old[i + root_n + 1];
                x_new[i] *= d_i[i];
                x_new[i] += d_i_s[i];

                //index root_n
                i = root_n;
                x_new[i] = ld[i] * x_old[i - root_n] +
                           lu[i] * x_old[i - (root_n - 1)] +
                           dl[i] * x_old[i - 1] +
                           //dd[i] * x_old[i] +
                           du[i] * x_old[ i + 1] +
                           ul[i] * x_old[i + root_n - 1] +
                           ud[i] * x_old[i + root_n] +
                           uu[i] * x_old[i + root_n + 1];
                x_new[i] *= d_i[i];
                x_new[i] += d_i_s[i];

                //index in [root_n + 1, n - (root_n + 1)[
                i = root_n + 1;
                for(; i < n - (root_n  + 1) ; ++i)
                {
                    x_new[i] = ll[i] * x_old[i - root_n - 1] +
                               ld[i] * x_old[i - root_n] +
                               lu[i] * x_old[i - root_n + 1] +
                               dl[i] * x_old[i - 1] +
                               //dd[i] * x_old[i] +
                               du[i] * x_old[ i + 1] +
                               ul[i] * x_old[i + root_n - 1] +
                               ud[i] * x_old[i + root_n] +
                               uu[i] * x_old[i + root_n + 1];
                    x_new[i] *= d_i[i];
                    x_new[i] += d_i_s[i];
                }

                //index n - (root_n + 1)
                i = n - (root_n + 1);
                x_new[i] = ll[i] * x_old[i - (n - (root_n + 1))] +
                           ld[i] * x_old[i - root_n] +
                           lu[i] * x_old[i - root_n + 1] +
                           dl[i] * x_old[i - 1] +
                           //dd[i] * x_old[i] +
                           du[i] * x_old[ i + 1] +
                           ul[i] * x_old[i + root_n - 1] +
                           ud[i] * x_old[i + root_n];
                x_new[i] *= d_i[i];
                x_new[i] += d_i_s[i];

                //index n - root_n
                i = n - root_n;
                x_new[i] = ll[i] * x_old[i - (n - (root_n + 1))] +
                           ld[i] * x_old[i - root_n] +
                           lu[i] * x_old[i - root_n + 1] +
                           dl[i] * x_old[i - 1] +
                           //dd[i] * x_old[i] +
                           du[i] * x_old[ i + 1] +
                           ul[i] * x_old[i + root_n - 1];
                x_new[i] *= d_i[i];
                x_new[i] += d_i_s[i];

                //index in [n - root_n + 1, n -1[
                i = n - root_n + 1;
                for(; i < n - 1; ++i)
                {
                    x_new[i] = ll[i] * x_old[i - (n - (root_n + 1))] +
                               ld[i] * x_old[i - root_n] +
                               lu[i] * x_old[i - root_n + 1] +
                               dl[i] * x_old[i - 1] +
                               //dd[i] * x_old[i] +
                               du[i] * x_old[ i + 1];
                    x_new[i] *= d_i[i];
                    x_new[i] += d_i_s[i];

                }

                //index n - 1
                i = n - 1;
                x_new[i] = ll[i] * x_old[i - (n - (root_n + 1))] +
                    ld[i] * x_old[i - root_n] +
                    lu[i] * x_old[i - root_n + 1] +
                    dl[i] * x_old[i - 1] ;//+
                    //dd[i] * x_old[i];
                x_new[i] *= d_i[i];
                x_new[i] += d_i_s[i];
                return result;
            }
    };
//----------------------------------------------------------------------------------------------------
    template<>
    struct JacobiKernelCascade<tags::CPU, TWICE, Q1>
    {
        public:
            template<typename DT1_, typename DT2_>
            static inline  DenseVector<DT1_> value(DenseVector<DT2_> & right_hand_side, DenseVector<DT1_> & former_result, DenseVector<DT1_> & diag_inverted, BandedMatrix<DT1_> & difference, DenseVector<DT2_> & scaled_d)
            {
                unsigned long n = right_hand_side.size();
                unsigned long root_n = (unsigned long)sqrt(n);
                DenseVector<DT1_> result(n);
                //Prepare result for SECOND iteration
                DenseVector<DT1_> result_next(n);

                DT1_ * x_old = former_result.elements();
                DT1_ * d_i = diag_inverted.elements();
                DT2_ * d_i_s = scaled_d.elements();
                DT1_ * x_new = result.elements();
                //SECOND iteration data
                DT1_ * x_new_2 = result_next.elements();

                DT1_ * ll = difference.band((signed long)(-root_n - 1)).elements();
                DT1_ * ld = difference.band((signed long)(-root_n)).elements();
                DT1_ * lu = difference.band((signed long)(-root_n + 1)).elements();

                DT1_ * dl = difference.band((signed long)(-1)).elements();
                //DT1_ * dd = difference.band((signed long)(0)).elements();
                DT1_ * du = difference.band((signed long)(1)).elements();

                DT1_ * ul = difference.band((signed long)(root_n - 1)).elements();
                DT1_ * ud = difference.band((signed long)(root_n)).elements();
                DT1_ * uu = difference.band((signed long)(root_n - 1)).elements();

                unsigned long i(0);
                //FIRST index 0
                x_new[i] = //dd[i] * x_old[i] +
                           du[i] * x_old[1] +
                           ul[i] * x_old[root_n - 1] +
                           ud[i] * x_old[root_n] +
                           uu[i] * x_old[root_n + 1];
                x_new[i] *= d_i[i];
                x_new[i] += d_i_s[i];

                //FIRST index in [1, root_n -1[
                i = 1;
                for(; i < root_n - 1 ; ++i)
                {
                    x_new[i] = dl[i] * x_old[i - 1] +
                               //dd[i] * x_old[i] +
                               du[i] * x_old[i + 1] +
                               ul[i] * x_old[i + root_n - 1] +
                               ud[i] * x_old[i + root_n] +
                               uu[i] * x_old[i + root_n + 1];
                    x_new[i] *= d_i[i];
                    x_new[i] += d_i_s[i];
                }

                //FIRST index root_n -1
                i = root_n - 1;
                x_new[i] = lu[i] * x_old[i - (root_n - 1)] +
                           dl[i] * x_old[i - 1] +
                           //dd[i] * x_old[i] +
                           du[i] * x_old[i + 1] +
                           ul[i] * x_old[i + root_n - 1] +
                           ud[i] * x_old[i + root_n] +
                           uu[i] * x_old[i + root_n + 1];
                x_new[i] *= d_i[i];
                x_new[i] += d_i_s[i];

                //FIRST index root_n
                i = root_n;
                x_new[i] = ld[i] * x_old[i - root_n] +
                           lu[i] * x_old[i - (root_n - 1)] +
                           dl[i] * x_old[i - 1] +
                           //dd[i] * x_old[i] +
                           du[i] * x_old[ i + 1] +
                           ul[i] * x_old[i + root_n - 1] +
                           ud[i] * x_old[i + root_n] +
                           uu[i] * x_old[i + root_n + 1];
                x_new[i] *= d_i[i];
                x_new[i] += d_i_s[i];

                //FIRST index in [root_n + 1, n - (root_n + 1)[
                //The loop is splitted in order to let the second iteration catch up
                i = root_n + 1;
                for(; i < 2*(root_n  + 1) ; ++i)
                {
                    x_new[i] = ll[i] * x_old[i - root_n - 1] +
                               ld[i] * x_old[i - root_n] +
                               lu[i] * x_old[i - root_n + 1] +
                               dl[i] * x_old[i - 1] +
                               du[i] * x_old[ i + 1] +
                               ul[i] * x_old[i + root_n - 1] +
                               ud[i] * x_old[i + root_n] +
                               uu[i] * x_old[i + root_n + 1];
                    x_new[i] *= d_i[i];
                    x_new[i] += d_i_s[i];
                }

                //SECOND index in [0, root_n[
                //SECOND index 0
                unsigned long j(0);
                x_new_2[j] = du[j] * x_new[1] +
                             ul[j] * x_new[root_n - 1] +
                             ud[j] * x_new[root_n] +
                             uu[j] * x_new[root_n + 1];
                x_new_2[j] *= d_i[j];
                x_new_2[j] += d_i_s[j];

                //SECOND index in [1, root_n -1[
                j = 1;
                for(; j < root_n - 1 ; ++j)
                {
                    x_new_2[j] = dl[j] * x_new[j - 1] +
                                 du[j] * x_new[j + 1] +
                                 ul[j] * x_new[j + root_n - 1] +
                                 ud[j] * x_new[j + root_n] +
                                 uu[j] * x_new[j + root_n + 1];
                    x_new_2[j] *= d_i[j];
                    x_new_2[j] += d_i_s[j];
                }

                //SECOND index root_n -1
                j = root_n - 1;
                x_new_2[j] = lu[j] * x_new[j - (root_n - 1)] +
                             dl[j] * x_new[j - 1] +
                             du[j] * x_new[j + 1] +
                             ul[j] * x_new[j + root_n - 1] +
                             ud[j] * x_new[j + root_n] +
                             uu[j] * x_new[j + root_n + 1];
                x_new_2[j] *= d_i[j];
                x_new_2[j] += d_i_s[j];

                //SECOND index root_n
                j = root_n;
                x_new_2[j] = ld[j] * x_new[j - root_n] +
                             lu[j] * x_new[j - (root_n - 1)] +
                             dl[j] * x_new[j - 1] +
                             du[j] * x_new[j + 1] +
                             ul[j] * x_new[j + root_n - 1] +
                             ud[j] * x_new[j + root_n] +
                             uu[j] * x_new[j + root_n + 1];
                x_new_2[j] *= d_i[j];
                x_new_2[j] += d_i_s[j];
                //END SECOND forecast

                //CONTINUE FIRST index
                for(; i < n - (root_n  + 1) ; ++i)
                {
                    x_new[i] = ll[i] * x_old[i - root_n - 1] +
                               ld[i] * x_old[i - root_n] +
                               lu[i] * x_old[i - root_n + 1] +
                               dl[i] * x_old[i - 1] +
                               //dd[i] * x_old[i] +
                               du[i] * x_old[ i + 1] +
                               ul[i] * x_old[i + root_n - 1] +
                               ud[i] * x_old[i + root_n] +
                               uu[i] * x_old[i + root_n + 1];
                    x_new[i] *= d_i[i];
                    x_new[i] += d_i_s[i];
                }

                //SECOND index in [root_n + 1, 2*(root_n +1)[
                j = root_n + 1;
                for(; j < 2*(root_n  + 1) ; ++j)
                {
                    x_new_2[j] = ll[j] * x_old[j - root_n - 1] +
                               ld[j] * x_old[j - root_n] +
                               lu[j] * x_old[j - root_n + 1] +
                               dl[j] * x_old[j - 1] +
                               du[j] * x_old[j + 1] +
                               ul[j] * x_old[j + root_n - 1] +
                               ud[j] * x_old[j + root_n] +
                               uu[j] * x_old[j + root_n + 1];
                    x_new_2[j] *= d_i[j];
                    x_new_2[j] += d_i_s[j];
                }
                //END SECOND forecast

                //FIRST index n - (root_n + 1)
                i = n - (root_n + 1);
                x_new[i] = ll[i] * x_old[i - (n - (root_n + 1))] +
                           ld[i] * x_old[i - root_n] +
                           lu[i] * x_old[i - root_n + 1] +
                           dl[i] * x_old[i - 1] +
                           //dd[i] * x_old[i] +
                           du[i] * x_old[ i + 1] +
                           ul[i] * x_old[i + root_n - 1] +
                           ud[i] * x_old[i + root_n];
                x_new[i] *= d_i[i];
                x_new[i] += d_i_s[i];

                //FIRST index n - root_n
                i = n - root_n;
                x_new[i] = ll[i] * x_old[i - (n - (root_n + 1))] +
                           ld[i] * x_old[i - root_n] +
                           lu[i] * x_old[i - root_n + 1] +
                           dl[i] * x_old[i - 1] +
                           //dd[i] * x_old[i] +
                           du[i] * x_old[ i + 1] +
                           ul[i] * x_old[i + root_n - 1];
                x_new[i] *= d_i[i];
                x_new[i] += d_i_s[i];

                //FIRST index in [n - root_n + 1, n -1[
                i = n - root_n + 1;
                for(; i < n - 1; ++i)
                {
                    x_new[i] = ll[i] * x_old[i - (n - (root_n + 1))] +
                               ld[i] * x_old[i - root_n] +
                               lu[i] * x_old[i - root_n + 1] +
                               dl[i] * x_old[i - 1] +
                               //dd[i] * x_old[i] +
                               du[i] * x_old[ i + 1];
                    x_new[i] *= d_i[i];
                    x_new[i] += d_i_s[i];

                }

                //FIRST index n - 1
                i = n - 1;
                x_new[i] = ll[i] * x_old[i - (n - (root_n + 1))] +
                    ld[i] * x_old[i - root_n] +
                    lu[i] * x_old[i - root_n + 1] +
                    dl[i] * x_old[i - 1] ;//+
                x_new[i] *= d_i[i];
                x_new[i] += d_i_s[i];

                //SECOND index rest
                for(; j < n - (root_n  + 1) ; ++j)
                {
                    x_new_2[j] = ll[j] * x_new[j - root_n - 1] +
                                 ld[j] * x_new[j - root_n] +
                                 lu[j] * x_new[j - root_n + 1] +
                                 dl[j] * x_new[j - 1] +
                                 du[j] * x_new[j + 1] +
                                 ul[j] * x_new[j + root_n - 1] +
                                 ud[j] * x_new[j + root_n] +
                                 uu[j] * x_new[j + root_n + 1];
                    x_new_2[j] *= d_i[j];
                    x_new_2[j] += d_i_s[j];
                }
                //SECOND index n - (root_n + 1)
                j = n - (root_n + 1);
                x_new_2[j] = ll[j] * x_new[j - (n - (root_n + 1))] +
                             ld[j] * x_new[j - root_n] +
                             lu[j] * x_new[j - root_n + 1] +
                             dl[j] * x_new[j - 1] +
                             du[j] * x_new[j + 1] +
                             ul[j] * x_new[j + root_n - 1] +
                             ud[j] * x_new[j + root_n];
                x_new_2[j] *= d_i[j];
                x_new_2[j] += d_i_s[j];

                //SECOND index n - root_n
                j = n - root_n;
                x_new_2[j] = ll[j] * x_old[j - (n - (root_n + 1))] +
                             ld[j] * x_old[j - root_n] +
                             lu[j] * x_old[j - root_n + 1] +
                             dl[j] * x_old[j - 1] +
                             du[j] * x_old[j + 1] +
                             ul[j] * x_old[j + root_n - 1];
                x_new_2[j] *= d_i[j];
                x_new_2[j] += d_i_s[j];

                //SECOND index in [n - root_n + 1, n -1[
                j = n - root_n + 1;
                for(; j < n - 1; ++j)
                {
                    x_new_2[j] = ll[j] * x_new[j - (n - (root_n + 1))] +
                                 ld[j] * x_new[j - root_n] +
                                 lu[j] * x_new[j - root_n + 1] +
                                 dl[j] * x_new[j - 1] +
                                 du[j] * x_new[j + 1];
                    x_new_2[j] *= d_i[j];
                    x_new_2[j] += d_i_s[j];
                }

                //SECOND index n - 1
                j = n - 1;
                x_new_2[j] = ll[j] * x_new[j - (n - (root_n + 1))] +
                             ld[j] * x_new[j - root_n] +
                             lu[j] * x_new[j - root_n + 1] +
                             dl[j] * x_new[j - 1] ;//+
                x_new_2[j] *= d_i[j];
                x_new_2[j] += d_i_s[j];
                return result_next;
            }
    };

    template<>
        struct JacobiKernelCascade<tags::CPU::SSE, ONCE, NONE>
        {
            public:
                static DenseVector<float> value(DenseVector<float> & b, DenseVector<float> & x, DenseVector<float> & d, BandedMatrix<float> & a, DenseVector<float> & scaled_d);
        };

    template<>
        struct JacobiKernelCascade<tags::CPU::SSE, ONCE, Q1>
        {
            public:
                static DenseVector<float> value(DenseVector<float> & b, DenseVector<float> & x, DenseVector<float> & d, BandedMatrix<float> & a, DenseVector<float> & scaled_d);
        };

}
#endif
