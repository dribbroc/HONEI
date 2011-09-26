/* vim: set sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
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

#include <honei/math/defect.hh>
#include <honei/backends/sse/operations.hh>
#include <honei/util/profiler.hh>

#include <xmmintrin.h>
#include <emmintrin.h>
#include <iostream>

using namespace std;
using namespace honei;

namespace honei
{
    DenseVector<float> & Defect<tags::CPU::SSE>::value(DenseVector<float> & result, const DenseVector<float> & right_hand_side, const BandedMatrixQx<Q1Type, float> & system, const DenseVector<float> & x)
    {
        unsigned long n = right_hand_side.size();
        unsigned long root_n = (unsigned long)sqrt(n);

        float * b = right_hand_side.elements();
        float * x_old = x.elements();
        float * x_new = result.elements();

        float * ll = system.band(LL).elements();
        float * ld = system.band(LD).elements();
        float * lu = system.band(LU).elements();

        float * dl = system.band(DL).elements();
        float * dd = system.band(DD).elements();
        float * du = system.band(DU).elements();

        float * ul = system.band(UL).elements();
        float * ud = system.band(UD).elements();
        float * uu = system.band(UU).elements();

        __m128 m1, m2, m3, m4;

        unsigned long index(0);

        //index 0
        x_new[index] = b[index] - (dd[index] * x_old[index]
                + du[index] * x_old[index + 1]
                + ul[index] * x_old[index + root_n - 1]
                + ud[index] * x_old[index + root_n]
                + uu[index] * x_old[index + root_n + 1]);

        //index in [1, root_n -1[
        index = 1;
        unsigned long quad_start(index + 4 - (index % 4));
        unsigned long quad_end(root_n - 1 - ((root_n - 1 - quad_start) % 4));

        for (unsigned long si(index) ; si < quad_start ; ++si)
        {
            x_new[si] = b[si] - (dd[si] * x_old[si]
                    + dl[si] * x_old[si - 1]
                    + du[si] * x_old[si + 1]
                    + ul[si] * x_old[si + root_n - 1]
                    + ud[si] * x_old[si + root_n]
                    + uu[si] * x_old[si + root_n + 1]);
        }

        for (unsigned long si(quad_start) ; si < quad_end ; si+=4)
        {
            m4 = _mm_load_ps(b + si);

            //x_new[si] = dd[si] * x_old[si]
            m2 = _mm_loadu_ps( x_old + si);
            m3 = _mm_load_ps(dd + si);
            m1 = _mm_mul_ps(m3, m2);

            //+ dl[si] * x_old[si - 1]
            m2 = _mm_loadu_ps( x_old + si - 1);
            m3 = _mm_load_ps(dl + si);
            m3 = _mm_mul_ps(m3, m2);
            m1 = _mm_add_ps(m1, m3);

            //+ du[si] * x_old[si + 1]
            m2 = _mm_loadu_ps( x_old + si + 1);
            m3 = _mm_load_ps(du + si);
            m3 = _mm_mul_ps(m3, m2);
            m1 = _mm_add_ps(m1, m3);

            //+ ul[si] * x_old[si + root_n - 1]
            m2 = _mm_loadu_ps( x_old + si + root_n - 1);
            m3 = _mm_load_ps(ul + si);
            m3 = _mm_mul_ps(m3, m2);
            m1 = _mm_add_ps(m1, m3);

            //+ ud[si] * x_old[si + root_n]
            m2 = _mm_loadu_ps( x_old + si + root_n);
            m3 = _mm_load_ps(ud + si);
            m3 = _mm_mul_ps(m3, m2);
            m1 = _mm_add_ps(m1, m3);

            //+ uu[si] * x_old[si + root_n + 1];
            m2 = _mm_loadu_ps( x_old + si + root_n + 1);
            m3 = _mm_load_ps(uu + si);
            m3 = _mm_mul_ps(m3, m2);
            m1 = _mm_add_ps(m1, m3);

            m1 = _mm_sub_ps(m4 , m1);

            _mm_store_ps(x_new + si, m1);
        }

        for (unsigned long si(quad_end) ; si < root_n - 1 ; ++si)
        {
            x_new[si] = b[si] - (dd[si] * x_old[si]
                    + dl[si] * x_old[si - 1]
                    + du[si] * x_old[si + 1]
                    + ul[si] * x_old[si + root_n - 1]
                    + ud[si] * x_old[si + root_n]
                    + uu[si] * x_old[si + root_n + 1]);
        }


        //index root_n -1
        index = root_n - 1;
        x_new[index] = b[index] - (dd[index] * x_old[index]
                + lu[index] * x_old[index - root_n + 1]
                + dl[index] * x_old[index - 1]
                + du[index] * x_old[index + 1]
                + ul[index] * x_old[index + root_n - 1]
                + ud[index] * x_old[index + root_n]
                + uu[index] * x_old[index + root_n + 1]);

        //index root_n
        index = root_n;
        x_new[index] = b[index] - (dd[index] * x_old[index]
                + lu[index] * x_old[index - root_n + 1]
                + ld[index] * x_old[index - root_n]
                + dl[index] * x_old[index - 1]
                + du[index] * x_old[index + 1]
                + ul[index] * x_old[index + root_n - 1]
                + ud[index] * x_old[index + root_n]
                + uu[index] * x_old[index + root_n + 1]);

        //index in [root_n + 1, n - (root_n + 1)[
        index = root_n + 1;
        quad_start = index + 4 - (index % 4);
        quad_end = n - root_n - 1 - ((n - root_n - 1 - quad_start) % 4);

        for (unsigned long si(index) ; si < quad_start ; ++si)
        {
            x_new[si] = b[si] - (dd[si] * x_old[si]
                    + ll[si] * x_old[si - root_n - 1]
                    + lu[si] * x_old[si - root_n + 1]
                    + ld[si] * x_old[si - root_n]
                    + dl[si] * x_old[si - 1]
                    + du[si] * x_old[si + 1]
                    + ul[si] * x_old[si + root_n - 1]
                    + ud[si] * x_old[si + root_n]
                    + uu[si] * x_old[si + root_n + 1]);
        }

        for (unsigned long si(quad_start) ; si < quad_end ; si+=4)
        {
            m4 = _mm_load_ps(b + si);

            //x_new[si] = dd[si] * x_old[si]
            m2 = _mm_loadu_ps( x_old + si);
            m3 = _mm_load_ps(dd + si);
            m1 = _mm_mul_ps(m3, m2);

            //+ ll[si] * x_old[si - root_n - 1]
            m2 = _mm_loadu_ps( x_old + si - root_n - 1);
            m3 = _mm_load_ps(ll + si);
            m3 = _mm_mul_ps(m3, m2);
            m1 = _mm_add_ps(m1, m3);

            //+ lu[si] * x_old[si - root_n + 1]
            m2 = _mm_loadu_ps( x_old + si - root_n + 1);
            m3 = _mm_load_ps(lu + si);
            m3 = _mm_mul_ps(m3, m2);
            m1 = _mm_add_ps(m1, m3);

            //+ ld[si] * x_old[si - root_n]
            m2 = _mm_loadu_ps( x_old + si - root_n);
            m3 = _mm_load_ps(ld + si);
            m3 = _mm_mul_ps(m3, m2);
            m1 = _mm_add_ps(m1, m3);

            //+ dl[si] * x_old[si - 1]
            m2 = _mm_loadu_ps( x_old + si - 1);
            m3 = _mm_load_ps(dl + si);
            m3 = _mm_mul_ps(m3, m2);
            m1 = _mm_add_ps(m1, m3);

            //+ du[si] * x_old[si + 1]
            m2 = _mm_loadu_ps( x_old + si + 1);
            m3 = _mm_load_ps(du + si);
            m3 = _mm_mul_ps(m3, m2);
            m1 = _mm_add_ps(m1, m3);

            //+ ul[si] * x_old[si + root_n - 1]
            m2 = _mm_loadu_ps( x_old + si + root_n - 1);
            m3 = _mm_load_ps(ul + si);
            m3 = _mm_mul_ps(m3, m2);
            m1 = _mm_add_ps(m1, m3);

            //+ ud[si] * x_old[si + root_n]
            m2 = _mm_loadu_ps( x_old + si + root_n);
            m3 = _mm_load_ps(ud + si);
            m3 = _mm_mul_ps(m3, m2);
            m1 = _mm_add_ps(m1, m3);

            //+ uu[si] * x_old[si + root_n + 1];
            m2 = _mm_loadu_ps( x_old + si + root_n + 1);
            m3 = _mm_load_ps(uu + si);
            m3 = _mm_mul_ps(m3, m2);
            m1 = _mm_add_ps(m1, m3);

            m1 = _mm_sub_ps(m4 , m1);

            _mm_store_ps(x_new + si, m1);
        }

        for (unsigned long si(quad_end) ; si < n - root_n - 1 ; ++si)
        {
            x_new[si] = b[si] - (dd[si] * x_old[si]
                    + ll[si] * x_old[si - root_n - 1]
                    + lu[si] * x_old[si - root_n + 1]
                    + ld[si] * x_old[si - root_n]
                    + dl[si] * x_old[si - 1]
                    + du[si] * x_old[si + 1]
                    + ul[si] * x_old[si + root_n - 1]
                    + ud[si] * x_old[si + root_n]
                    + uu[si] * x_old[si + root_n + 1]);
        }

        //index n - (root_n + 1)
        index = n - root_n - 1;
        x_new[index] = b[index] - (dd[index] * x_old[index]
                + ll[index] * x_old[index - root_n - 1]
                + lu[index] * x_old[index - root_n + 1]
                + ld[index] * x_old[index - root_n]
                + dl[index] * x_old[index - 1]
                + du[index] * x_old[index + 1]
                + ul[index] * x_old[index + root_n - 1]
                + ud[index] * x_old[index + root_n]);

        //index n - root_n
        index = n - root_n;
        x_new[index] = b[index] - (dd[index] * x_old[index]
                + ll[index] * x_old[index - root_n - 1]
                + lu[index] * x_old[index - root_n + 1]
                + ld[index] * x_old[index - root_n]
                + dl[index] * x_old[index - 1]
                + du[index] * x_old[index + 1]
                + ul[index] * x_old[index + root_n - 1]);

        //index in [n - root_n + 1, n -1[
        index = n - root_n + 1;
        quad_start = index + 4 - (index % 4);
        quad_end = n - 1 - ((n - 1 - quad_start) % 4);

        for (unsigned long si(index) ; si < quad_start ; ++si)
        {
            x_new[si] = b[si] - (dd[si] * x_old[si]
                    + ll[si] * x_old[si - root_n - 1]
                    + lu[si] * x_old[si - root_n + 1]
                    + ld[si] * x_old[si - root_n]
                    + dl[si] * x_old[si - 1]
                    + du[si] * x_old[si + 1]);
        }

        for (unsigned long si(quad_start) ; si < quad_end ; si+=4)
        {
            m4 = _mm_load_ps(b + si);

            //x_new[si] = dd[si] * x_old[si]
            m2 = _mm_loadu_ps( x_old + si);
            m3 = _mm_load_ps(dd + si);
            m1 = _mm_mul_ps(m3, m2);

            //+ ll[si] * x_old[si - root_n - 1]
            m2 = _mm_loadu_ps( x_old + si - root_n - 1);
            m3 = _mm_load_ps(ll + si);
            m3 = _mm_mul_ps(m3, m2);
            m1 = _mm_add_ps(m1, m3);

            //+ lu[si] * x_old[si - root_n + 1]
            m2 = _mm_loadu_ps( x_old + si - root_n + 1);
            m3 = _mm_load_ps(lu + si);
            m3 = _mm_mul_ps(m3, m2);
            m1 = _mm_add_ps(m1, m3);

            //+ ld[si] * x_old[si - root_n]
            m2 = _mm_loadu_ps( x_old + si - root_n);
            m3 = _mm_load_ps(ld + si);
            m3 = _mm_mul_ps(m3, m2);
            m1 = _mm_add_ps(m1, m3);

            //+ dl[si] * x_old[si - 1]
            m2 = _mm_loadu_ps( x_old + si - 1);
            m3 = _mm_load_ps(dl + si);
            m3 = _mm_mul_ps(m3, m2);
            m1 = _mm_add_ps(m1, m3);

            //+ du[si] * x_old[si + 1]
            m2 = _mm_loadu_ps( x_old + si + 1);
            m3 = _mm_load_ps(du + si);
            m3 = _mm_mul_ps(m3, m2);
            m1 = _mm_add_ps(m1, m3);

                m1 = _mm_sub_ps(m4 , m1);

            _mm_store_ps(x_new + si, m1);
        }

        for (unsigned long si(quad_end) ; si < n - root_n - 1 ; ++si)
        {
            x_new[si] = b[si] - (dd[si] * x_old[si]
                    + ll[si] * x_old[si - root_n - 1]
                    + lu[si] * x_old[si - root_n + 1]
                    + ld[si] * x_old[si - root_n]
                    + dl[si] * x_old[si - 1]
                    + du[si] * x_old[si + 1]);
        }

        //index n - 1
        index = n - 1;
        x_new[index] = b[index] - (dd[index] * x_old[index]
                + ll[index] * x_old[index - root_n - 1]
                + lu[index] * x_old[index - root_n + 1]
                + ld[index] * x_old[index - root_n]
                + dl[index] * x_old[index - 1]);

        return result;
    }

    DenseVector<float> Defect<tags::CPU::SSE>::value(const DenseVector<float> & right_hand_side, const BandedMatrixQx<Q1Type, float> & system, const DenseVector<float> & x)
    {
        DenseVector<float> result(right_hand_side.size());
        Defect<tags::CPU::SSE>::value(result, right_hand_side, system, x);
        return result;
    }

    DenseVector<double> & Defect<tags::CPU::SSE>::value(DenseVector<double> & result, const DenseVector<double> & right_hand_side, const BandedMatrixQx<Q1Type, double> & system, const DenseVector<double> & x)
    {
        unsigned long n = right_hand_side.size();
        unsigned long root_n = (unsigned long)sqrt(n);

        double * b = right_hand_side.elements();
        double * x_old = x.elements();
        double * x_new = result.elements();

        double * ll = system.band(LL).elements();
        double * ld = system.band(LD).elements();
        double * lu = system.band(LU).elements();

        double * dl = system.band(DL).elements();
        double * dd = system.band(DD).elements();
        double * du = system.band(DU).elements();

        double * ul = system.band(UL).elements();
        double * ud = system.band(UD).elements();
        double * uu = system.band(UU).elements();

        __m128d m1, m2, m3, m4;

        unsigned long index(0);

        //index 0
        x_new[index] = b[index] - (dd[index] * x_old[index]
                + du[index] * x_old[index + 1]
                + ul[index] * x_old[index + root_n - 1]
                + ud[index] * x_old[index + root_n]
                + uu[index] * x_old[index + root_n + 1]);

        //index in [1, root_n -1[
        index = 1;
        unsigned long quad_start(index + (index % 2));
        unsigned long quad_end(root_n - 1 - ((root_n - 1 - quad_start) % 2));

        for (unsigned long si(index) ; si < quad_start ; ++si)
        {
            x_new[si] = b[si] - (dd[si] * x_old[si]
                    + dl[si] * x_old[si - 1]
                    + du[si] * x_old[si + 1]
                    + ul[si] * x_old[si + root_n - 1]
                    + ud[si] * x_old[si + root_n]
                    + uu[si] * x_old[si + root_n + 1]);
        }

        for (unsigned long si(quad_start) ; si < quad_end ; si+=2)
        {
            m4 = _mm_load_pd(b + si);

            //x_new[si] = dd[si] * x_old[si]
            m2 = _mm_loadu_pd( x_old + si);
            m3 = _mm_load_pd(dd + si);
            m1 = _mm_mul_pd(m3, m2);

            //+ dl[si] * x_old[si - 1]
            m2 = _mm_loadu_pd( x_old + si - 1);
            m3 = _mm_load_pd(dl + si);
            m3 = _mm_mul_pd(m3, m2);
            m1 = _mm_add_pd(m1, m3);

            //+ du[si] * x_old[si + 1]
            m2 = _mm_loadu_pd( x_old + si + 1);
            m3 = _mm_load_pd(du + si);
            m3 = _mm_mul_pd(m3, m2);
            m1 = _mm_add_pd(m1, m3);

            //+ ul[si] * x_old[si + root_n - 1]
            m2 = _mm_loadu_pd( x_old + si + root_n - 1);
            m3 = _mm_load_pd(ul + si);
            m3 = _mm_mul_pd(m3, m2);
            m1 = _mm_add_pd(m1, m3);

            //+ ud[si] * x_old[si + root_n]
            m2 = _mm_loadu_pd( x_old + si + root_n);
            m3 = _mm_load_pd(ud + si);
            m3 = _mm_mul_pd(m3, m2);
            m1 = _mm_add_pd(m1, m3);

            //+ uu[si] * x_old[si + root_n + 1];
            m2 = _mm_loadu_pd( x_old + si + root_n + 1);
            m3 = _mm_load_pd(uu + si);
            m3 = _mm_mul_pd(m3, m2);
            m1 = _mm_add_pd(m1, m3);

            m1 = _mm_sub_pd(m4 , m1);

            _mm_store_pd(x_new + si, m1);
        }

        for (unsigned long si(quad_end) ; si < root_n - 1 ; ++si)
        {
            x_new[si] = b[si] - (dd[si] * x_old[si]
                    + dl[si] * x_old[si - 1]
                    + du[si] * x_old[si + 1]
                    + ul[si] * x_old[si + root_n - 1]
                    + ud[si] * x_old[si + root_n]
                    + uu[si] * x_old[si + root_n + 1]);
        }


        //index root_n -1
        index = root_n - 1;
        x_new[index] = b[index] - (dd[index] * x_old[index]
                + lu[index] * x_old[index - root_n + 1]
                + dl[index] * x_old[index - 1]
                + du[index] * x_old[index + 1]
                + ul[index] * x_old[index + root_n - 1]
                + ud[index] * x_old[index + root_n]
                + uu[index] * x_old[index + root_n + 1]);

        //index root_n
        index = root_n;
        x_new[index] = b[index] - (dd[index] * x_old[index]
                + lu[index] * x_old[index - root_n + 1]
                + ld[index] * x_old[index - root_n]
                + dl[index] * x_old[index - 1]
                + du[index] * x_old[index + 1]
                + ul[index] * x_old[index + root_n - 1]
                + ud[index] * x_old[index + root_n]
                + uu[index] * x_old[index + root_n + 1]);

        //index in [root_n + 1, n - (root_n + 1)[
        index = root_n + 1;
        quad_start = index + (index % 2);
        quad_end = n - root_n - 1 - ((n - root_n - 1 - quad_start) % 2);

        for (unsigned long si(index) ; si < quad_start ; ++si)
        {
            x_new[si] = b[si] - (dd[si] * x_old[si]
                    + ll[si] * x_old[si - root_n - 1]
                    + lu[si] * x_old[si - root_n + 1]
                    + ld[si] * x_old[si - root_n]
                    + dl[si] * x_old[si - 1]
                    + du[si] * x_old[si + 1]
                    + ul[si] * x_old[si + root_n - 1]
                    + ud[si] * x_old[si + root_n]
                    + uu[si] * x_old[si + root_n + 1]);
        }

        for (unsigned long si(quad_start) ; si < quad_end ; si+=2)
        {
            m4 = _mm_load_pd(b + si);

            //x_new[si] = dd[si] * x_old[si]
            m2 = _mm_loadu_pd( x_old + si);
            m3 = _mm_load_pd(dd + si);
            m1 = _mm_mul_pd(m3, m2);

            //+ ll[si] * x_old[si - root_n - 1]
            m2 = _mm_loadu_pd( x_old + si - root_n - 1);
            m3 = _mm_load_pd(ll + si);
            m3 = _mm_mul_pd(m3, m2);
            m1 = _mm_add_pd(m1, m3);

            //+ lu[si] * x_old[si - root_n + 1]
            m2 = _mm_loadu_pd( x_old + si - root_n + 1);
            m3 = _mm_load_pd(lu + si);
            m3 = _mm_mul_pd(m3, m2);
            m1 = _mm_add_pd(m1, m3);

            //+ ld[si] * x_old[si - root_n]
            m2 = _mm_loadu_pd( x_old + si - root_n);
            m3 = _mm_load_pd(ld + si);
            m3 = _mm_mul_pd(m3, m2);
            m1 = _mm_add_pd(m1, m3);

            //+ dl[si] * x_old[si - 1]
            m2 = _mm_loadu_pd( x_old + si - 1);
            m3 = _mm_load_pd(dl + si);
            m3 = _mm_mul_pd(m3, m2);
            m1 = _mm_add_pd(m1, m3);

            //+ du[si] * x_old[si + 1]
            m2 = _mm_loadu_pd( x_old + si + 1);
            m3 = _mm_load_pd(du + si);
            m3 = _mm_mul_pd(m3, m2);
            m1 = _mm_add_pd(m1, m3);

            //+ ul[si] * x_old[si + root_n - 1]
            m2 = _mm_loadu_pd( x_old + si + root_n - 1);
            m3 = _mm_load_pd(ul + si);
            m3 = _mm_mul_pd(m3, m2);
            m1 = _mm_add_pd(m1, m3);

            //+ ud[si] * x_old[si + root_n]
            m2 = _mm_loadu_pd( x_old + si + root_n);
            m3 = _mm_load_pd(ud + si);
            m3 = _mm_mul_pd(m3, m2);
            m1 = _mm_add_pd(m1, m3);

            //+ uu[si] * x_old[si + root_n + 1];
            m2 = _mm_loadu_pd( x_old + si + root_n + 1);
            m3 = _mm_load_pd(uu + si);
            m3 = _mm_mul_pd(m3, m2);
            m1 = _mm_add_pd(m1, m3);

            m1 = _mm_sub_pd(m4 , m1);

            _mm_store_pd(x_new + si, m1);
        }

        for (unsigned long si(quad_end) ; si < n - root_n - 1 ; ++si)
        {
            x_new[si] = b[si] - (dd[si] * x_old[si]
                    + ll[si] * x_old[si - root_n - 1]
                    + lu[si] * x_old[si - root_n + 1]
                    + ld[si] * x_old[si - root_n]
                    + dl[si] * x_old[si - 1]
                    + du[si] * x_old[si + 1]
                    + ul[si] * x_old[si + root_n - 1]
                    + ud[si] * x_old[si + root_n]
                    + uu[si] * x_old[si + root_n + 1]);
        }

        //index n - (root_n + 1)
        index = n - root_n - 1;
        x_new[index] = b[index] - (dd[index] * x_old[index]
                + ll[index] * x_old[index - root_n - 1]
                + lu[index] * x_old[index - root_n + 1]
                + ld[index] * x_old[index - root_n]
                + dl[index] * x_old[index - 1]
                + du[index] * x_old[index + 1]
                + ul[index] * x_old[index + root_n - 1]
                + ud[index] * x_old[index + root_n]);

        //index n - root_n
        index = n - root_n;
        x_new[index] = b[index] - (dd[index] * x_old[index]
                + ll[index] * x_old[index - root_n - 1]
                + lu[index] * x_old[index - root_n + 1]
                + ld[index] * x_old[index - root_n]
                + dl[index] * x_old[index - 1]
                + du[index] * x_old[index + 1]
                + ul[index] * x_old[index + root_n - 1]);

        //index in [n - root_n + 1, n -1[
        index = n - root_n + 1;
        quad_start = index + (index % 2);
        quad_end = n - 1 - ((n - 1 - quad_start) % 2);

        for (unsigned long si(index) ; si < quad_start ; ++si)
        {
            x_new[si] = b[si] - (dd[si] * x_old[si]
                    + ll[si] * x_old[si - root_n - 1]
                    + lu[si] * x_old[si - root_n + 1]
                    + ld[si] * x_old[si - root_n]
                    + dl[si] * x_old[si - 1]
                    + du[si] * x_old[si + 1]);
        }

        for (unsigned long si(quad_start) ; si < quad_end ; si+=2)
        {
            m4 = _mm_load_pd(b + si);

            //x_new[si] = dd[si] * x_old[si]
            m2 = _mm_loadu_pd( x_old + si);
            m3 = _mm_load_pd(dd + si);
            m1 = _mm_mul_pd(m3, m2);

            //+ ll[si] * x_old[si - root_n - 1]
            m2 = _mm_loadu_pd( x_old + si - root_n - 1);
            m3 = _mm_load_pd(ll + si);
            m3 = _mm_mul_pd(m3, m2);
            m1 = _mm_add_pd(m1, m3);

            //+ lu[si] * x_old[si - root_n + 1]
            m2 = _mm_loadu_pd( x_old + si - root_n + 1);
            m3 = _mm_load_pd(lu + si);
            m3 = _mm_mul_pd(m3, m2);
            m1 = _mm_add_pd(m1, m3);

            //+ ld[si] * x_old[si - root_n]
            m2 = _mm_loadu_pd( x_old + si - root_n);
            m3 = _mm_load_pd(ld + si);
            m3 = _mm_mul_pd(m3, m2);
            m1 = _mm_add_pd(m1, m3);

            //+ dl[si] * x_old[si - 1]
            m2 = _mm_loadu_pd( x_old + si - 1);
            m3 = _mm_load_pd(dl + si);
            m3 = _mm_mul_pd(m3, m2);
            m1 = _mm_add_pd(m1, m3);

            //+ du[si] * x_old[si + 1]
            m2 = _mm_loadu_pd( x_old + si + 1);
            m3 = _mm_load_pd(du + si);
            m3 = _mm_mul_pd(m3, m2);
            m1 = _mm_add_pd(m1, m3);

            m1 = _mm_sub_pd(m4 , m1);

            _mm_store_pd(x_new + si, m1);
        }

        for (unsigned long si(quad_end) ; si < n - root_n - 1 ; ++si)
        {
            x_new[si] = b[si] - (dd[si] * x_old[si]
                    + ll[si] * x_old[si - root_n - 1]
                    + lu[si] * x_old[si - root_n + 1]
                    + ld[si] * x_old[si - root_n]
                    + dl[si] * x_old[si - 1]
                    + du[si] * x_old[si + 1]);
        }

        //index n - 1
        index = n - 1;
        x_new[index] = b[index] - (dd[index] * x_old[index]
                + ll[index] * x_old[index - root_n - 1]
                + lu[index] * x_old[index - root_n + 1]
                + ld[index] * x_old[index - root_n]
                + dl[index] * x_old[index - 1]);

        return result;
    }

    DenseVector<double> Defect<tags::CPU::SSE>::value(const DenseVector<double> & right_hand_side, const BandedMatrixQx<Q1Type, double> & system, const DenseVector<double> & x)
    {
        DenseVector<double> result(right_hand_side.size());
        Defect<tags::CPU::SSE>::value(result, right_hand_side, system, x);
        return result;
    }

    DenseVector<float> & Defect<tags::CPU::SSE>::value(DenseVector<float> & result, const DenseVector<float> & right_hand_side, const SparseMatrixELL<float> & a, const DenseVector<float> & b,
            unsigned long row_start, unsigned long row_end)
    {
        CONTEXT("When calculating defect of SparseMatrixELL<float> with DenseVector<float> (SSE):");
        PROFILER_START("Defect SMELL float tags::CPU::SSE");

        if (b.size() != a.columns())
        {
            throw VectorSizeDoesNotMatch(b.size(), a.columns());
        }
        if (right_hand_side.size() != result.size())
        {
            throw VectorSizeDoesNotMatch(result.size(), right_hand_side.size());
        }


        if (row_end == 0)
        {
            fill<tags::CPU::SSE>(result, float(0));
            honei::sse::defect_smell_dv(result.elements(), right_hand_side.elements(), a.Aj().elements(), a.Ax().elements(), b.elements(),
                    a.stride(), a.rows(), a.num_cols_per_row(), a.threads());
        }

        else
        {
            honei::sse::defect_smell_dv(result.elements(), right_hand_side.elements(), a.Aj().elements(), a.Ax().elements(), a.Arl().elements(), b.elements(),
                    a.stride(), a.rows(), a.num_cols_per_row(), row_start, row_end, a.threads());
        }

        PROFILER_STOP("Defect SMELL float tags::CPU::SSE");
        return result;
    }

    DenseVector<double> & Defect<tags::CPU::SSE>::value(DenseVector<double> & result, const DenseVector<double> & right_hand_side, const SparseMatrixELL<double> & a, const DenseVector<double> & b,
            unsigned long row_start, unsigned long row_end)
    {
        CONTEXT("When calculating defect of SparseMatrixELL<double> with DenseVector<double> (SSE):");
        PROFILER_START("Defect SMELL double tags::CPU::SSE");

        if (b.size() != a.columns())
        {
            throw VectorSizeDoesNotMatch(b.size(), a.columns());
        }
        if (right_hand_side.size() != result.size())
        {
            throw VectorSizeDoesNotMatch(result.size(), right_hand_side.size());
        }


        if (row_end == 0)
        {
            fill<tags::CPU::SSE>(result, double(0));
            honei::sse::defect_smell_dv(result.elements(), right_hand_side.elements(), a.Aj().elements(), a.Ax().elements(), b.elements(),
                    a.stride(), a.rows(), a.num_cols_per_row(), a.threads());
        }

        else
        {
            honei::sse::defect_smell_dv(result.elements(), right_hand_side.elements(), a.Aj().elements(), a.Ax().elements(), a.Arl().elements(), b.elements(),
                    a.stride(), a.rows(), a.num_cols_per_row(), row_start, row_end, a.threads());
        }

        PROFILER_STOP("Defect SMELL double tags::CPU::SSE");
        return result;
    }
}
