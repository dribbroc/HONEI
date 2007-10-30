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

#include <libmath/iterative_refinement.hh>
#include <unittest/unittest.hh>
#include <libutil/stringify.hh>
#include <iostream>

using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DT1_>
class IterativeRefinementTestDenseCG:
    public BaseTest
{
    public:
        IterativeRefinementTestDenseCG(const std::string & tag) :
            BaseTest("Iterative Refinement with CG solver test (Dense system)<" + tag + ">")
        {
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
            DenseVector<DT1_> result(IterativeRefinement<CG, tags::CPU>::value(A,b,double(0.01), double(0.01)));
            DT1_ x_n = Norm< vnt_l_two, false, DT1_>::value(result);
            DenseVector<DT1_> x_analytical(3, DT1_(0));
            cout<<"RESULT(v1):"<<result<<endl;

            x_analytical[0] = DT1_(2./3.);
            x_analytical[1] = DT1_(5./6.);
            x_analytical[2] = DT1_(-1./3.);
            DT1_ x_analytical_n = Norm< vnt_l_two, false, DT1_>::value(x_analytical);
            TEST_CHECK_EQUAL_WITHIN_EPS(x_analytical_n, x_n , double(0.1));

        }
};

template <typename Tag_, typename DT1_>
class IterativeRefinementTestDenseJAC:
    public BaseTest
{
    public:
        IterativeRefinementTestDenseJAC(const std::string & tag) :
            BaseTest("Iterative Refinement with JAC solver test (Dense system)<" + tag + ">")
        {
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
            DenseVector<DT1_> result(IterativeRefinement<JAC, tags::CPU>::value(A,b,double(0.0001), double(0.0001)));
            DT1_ x_n = Norm< vnt_l_two, false, DT1_>::value(result);
            DenseVector<DT1_> x_analytical(3, DT1_(0));
            cout<<"RESULT(v1):"<<result<<endl;

            x_analytical[0] = DT1_(2./3.);
            x_analytical[1] = DT1_(5./6.);
            x_analytical[2] = DT1_(-1./3.);
            DT1_ x_analytical_n = Norm< vnt_l_two, false, DT1_>::value(x_analytical);
            TEST_CHECK_EQUAL_WITHIN_EPS(x_analytical_n, x_n , double(0.1));

        }
};

IterativeRefinementTestDenseCG<tags::CPU, double> iterref_test_double_denseCG("double");
IterativeRefinementTestDenseJAC<tags::CPU, double> iterref_test_double_denseJAC("double");

