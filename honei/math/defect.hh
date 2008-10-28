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
            static DenseVector<DT_> value(DenseVector<DT_> & right_hand_side, BandedMatrixQ1<DT_> & system, DenseVector<DT_> & x)
            {
                if (x.size() != system.columns())
                {
                    throw VectorSizeDoesNotMatch(x.size(), system.columns());
                }
                if (right_hand_side.size() != system.columns())
                {
                    throw VectorSizeDoesNotMatch(right_hand_side.size(), system.columns());
                }

                right_hand_side.lock(lm_read_only);
                system.lock(lm_read_only);
                x.lock(lm_read_only);
                DenseVector<DT_> result(right_hand_side.size());
                result.lock(lm_write_only);
                unsigned long n = right_hand_side.size();
                unsigned long root_n = (unsigned long)sqrt(n);

                DT_ * rhs = right_hand_side.elements();
                DT_ * x_old = x.elements();
                DT_ * x_new = result.elements();

                DT_ * ll = system.band(LL).elements();
                DT_ * ld = system.band(LD).elements();
                DT_ * lu = system.band(LU).elements();

                DT_ * dl = system.band(DL).elements();
                DT_ * dd = system.band(DD).elements();
                DT_ * du = system.band(DU).elements();

                DT_ * ul = system.band(UL).elements();
                DT_ * ud = system.band(UD).elements();
                DT_ * uu = system.band(UU).elements();

                unsigned long i(0);
                //index 0
                x_new[i] = ((rhs [i] - (dd[i] * x_old[i] +
                           du[i] * x_old[1] +
                           ul[i] * x_old[root_n - 1] +
                           ud[i] * x_old[root_n] +
                           uu[i] * x_old[root_n + 1])));

                //index in [1, root_n -1[
                i = 1;
                for(; i < root_n - 1 ; ++i)
                {
                    x_new[i] = ((rhs[i] - (dl[i] * x_old[i - 1] +
                               dd[i] * x_old[i] +
                               du[i] * x_old[i + 1] +
                               ul[i] * x_old[i + root_n - 1] +
                               ud[i] * x_old[i + root_n] +
                               uu[i] * x_old[i + root_n + 1])));
                }

                //index root_n -1
                i = root_n - 1;
                x_new[i] = ((rhs[i] - (lu[i] * x_old[i - (root_n - 1)] +
                           dl[i] * x_old[i - 1] +
                           dd[i] * x_old[i] +
                           du[i] * x_old[i + 1] +
                           ul[i] * x_old[i + root_n - 1] +
                           ud[i] * x_old[i + root_n] +
                           uu[i] * x_old[i + root_n + 1])));

                //index root_n
                i = root_n;
                x_new[i] = ((rhs[i] - (ld[i] * x_old[i - root_n] +
                           lu[i] * x_old[i - (root_n - 1)] +
                           dl[i] * x_old[i - 1] +
                           dd[i] * x_old[i] +
                           du[i] * x_old[ i + 1] +
                           ul[i] * x_old[i + root_n - 1] +
                           ud[i] * x_old[i + root_n] +
                           uu[i] * x_old[i + root_n + 1])));

                //index in [root_n + 1, n - (root_n + 1)[
                i = root_n + 1;
                for(; i < n - (root_n  + 1) ; ++i)
                {
                    x_new[i] = ((rhs[i] - (ll[i] * x_old[i - root_n - 1] +
                               ld[i] * x_old[i - root_n] +
                               lu[i] * x_old[i - root_n + 1] +
                               dl[i] * x_old[i - 1] +
                               dd[i] * x_old[i] +
                               du[i] * x_old[ i + 1] +
                               ul[i] * x_old[i + root_n - 1] +
                               ud[i] * x_old[i + root_n] +
                               uu[i] * x_old[i + root_n + 1])));
                }

                //index n - (root_n + 1)
                i = n - (root_n + 1);
                x_new[i] = ((rhs[i] - (ll[i] * x_old[i - (n - (root_n + 1))] +
                           ld[i] * x_old[i - root_n] +
                           lu[i] * x_old[i - root_n + 1] +
                           dl[i] * x_old[i - 1] +
                           dd[i] * x_old[i] +
                           du[i] * x_old[ i + 1] +
                           ul[i] * x_old[i + root_n - 1] +
                           ud[i] * x_old[i + root_n])));

                //index n - root_n
                i = n - root_n;
                x_new[i] = ((rhs[i] - (ll[i] * x_old[i - (n - (root_n + 1))] +
                           ld[i] * x_old[i - root_n] +
                           lu[i] * x_old[i - root_n + 1] +
                           dl[i] * x_old[i - 1] +
                           dd[i] * x_old[i] +
                           du[i] * x_old[ i + 1] +
                           ul[i] * x_old[i + root_n - 1])));

                //index in [n - root_n + 1, n -1[
                i = n - root_n + 1;
                for(; i < n - 1; ++i)
                {
                    x_new[i] = ((rhs[i] - (ll[i] * x_old[i - (n - (root_n + 1))] +
                               ld[i] * x_old[i - root_n] +
                               lu[i] * x_old[i - root_n + 1] +
                               dl[i] * x_old[i - 1] +
                               dd[i] * x_old[i] +
                               du[i] * x_old[ i + 1])));
                }

                //index n - 1
                i = n - 1;
                x_new[i] = ((rhs[i] - (ll[i] * x_old[i - (n - (root_n + 1))] +
                    ld[i] * x_old[i - root_n] +
                    lu[i] * x_old[i - root_n + 1] +
                    dl[i] * x_old[i - 1] +
                    dd[i] * x_old[i])));

                right_hand_side.unlock(lm_read_only);
                system.unlock(lm_read_only);
                x.unlock(lm_read_only);
                result.unlock(lm_write_only);
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
    template<>
        struct Defect<tags::CPU::SSE>
        {

            public:
                static DenseVector<double> value(DenseVector<double> & right_hand_side, BandedMatrixQ1<double> & system, DenseVector<double> & x);
                static DenseVector<float> value(DenseVector<float> & right_hand_side, BandedMatrixQ1<float> & system, DenseVector<float> & x);
        };
}

#endif
