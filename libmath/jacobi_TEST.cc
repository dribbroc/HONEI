/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the Math C++ library. LibMath is free software;
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

#include <libmath/jacobi.hh>
#include <unittest/unittest.hh>
#include <libutil/stringify.hh>
#include <iostream>

using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DT1_>
class JacobiTestDense:
    public BaseTest
{
    public:
        JacobiTestDense(const std::string & tag) :
            BaseTest("Jacobi solver test (Dense system)<" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            DenseMatrix<DT1_> A(3, 3, DT1_(1));
            DenseVector<DT1_> b(3, DT1_(1));
            A[0][0] = DT1_(7);
            A[0][1] = DT1_(-2);
            A[0][2] = DT1_(0);
            A[1][0] = DT1_(-2);
            A[1][1] = DT1_(6);
            A[1][2] = DT1_(2);
            A[2][0] = DT1_(0);
            A[2][1] = DT1_(2);
            A[2][2] = DT1_(5);

            b[0] = DT1_(3);
            b[1] = DT1_(3);
            b[2] = DT1_(0);

            std::cout<<"A:"<<A<<endl;
            std::cout<<"b:"<<b<<endl;
            DenseVector<DT1_> result = Jacobi<Tag_>::value(A,b,long(12));
            DT1_ x_n = Norm< vnt_l_two, false, DT1_>::value(result);
            DenseVector<DT1_> x_analytical(3, DT1_(0));
            cout<<"RESULT(v1):"<<result<<endl;

            x_analytical[0] = DT1_(2./3.);
            x_analytical[1] = DT1_(5./6.);
            x_analytical[2] = DT1_(-1./3.);
            DT1_ x_analytical_n = Norm< vnt_l_two, false, DT1_>::value(x_analytical);
            TEST_CHECK_EQUAL_WITHIN_EPS(x_analytical_n, x_n , double(0.1));

/*            DenseVector<DT1_> result_2 = Jacobi<Tag_>::value(A,b,double(0.1));
            cout<<"RESULT(v2):"<<result_2<<endl;

            DT1_ x_n_2 = Norm< vnt_l_two, false, DT1_>::value(result_2);
            TEST_CHECK_EQUAL_WITHIN_EPS(x_analytical_n, x_n_2 , double(0.1));
*/
        }
};

JacobiTestDense<tags::CPU, float> jacobi_test_float_dense("float");
JacobiTestDense<tags::CPU, double> jacobi_test_double_dense("double");

#ifdef HONEI_SSE
JacobiTestDense<tags::CPU::SSE, float> sse_jacobi_test_float_dense("SSE float");
JacobiTestDense<tags::CPU::SSE, double> sse_jacobi_test_double_dense("SSE double");
#endif
