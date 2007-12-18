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
            DenseVector<DT1_> result(IterativeRefinement<CG, tags::CPU>::value(A,b,double(std::numeric_limits<double>::epsilon()), double(std::numeric_limits<double>::epsilon())));
            DT1_ x_n = Norm< vnt_l_two, false, Tag_>::value(result);
            DenseVector<DT1_> x_analytical(3, DT1_(0));
            cout<<"RESULT(v1):"<<result<<endl;

            x_analytical[0] = DT1_(2./3.);
            x_analytical[1] = DT1_(5./6.);
            x_analytical[2] = DT1_(-1./3.);
            DT1_ x_analytical_n = Norm< vnt_l_two, false, Tag_>::value(x_analytical);
            TEST_CHECK_EQUAL_WITHIN_EPS(x_analytical_n, x_n , double(0.001));

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
            DenseVector<DT1_> result(IterativeRefinement<JAC, tags::CPU>::value(A,b,double(std::numeric_limits<double>::epsilon()), double(std::numeric_limits<double>::epsilon())));
            DT1_ x_n = Norm< vnt_l_two, false, Tag_>::value(result);
            DenseVector<DT1_> x_analytical(3, DT1_(0));
            cout<<"RESULT(v1):"<<result<<endl;

            x_analytical[0] = DT1_(2./3.);
            x_analytical[1] = DT1_(5./6.);
            x_analytical[2] = DT1_(-1./3.);
            DT1_ x_analytical_n = Norm< vnt_l_two, false, Tag_>::value(x_analytical);
            TEST_CHECK_EQUAL_WITHIN_EPS(x_analytical_n, x_n , double(0.001));

        }
};

template <typename Tag_, typename DT1_>
class IterativeRefinementTestBandedCG:
    public BaseTest
{
    public:
        IterativeRefinementTestBandedCG(const std::string & tag) :
            BaseTest("Iterative Refinement with CG solver test (Banded system)<" + tag + ">")
        {
        }

        virtual void run() const
        {

            DenseVector<DT1_> diag(3, DT1_(1));
            DenseVector<DT1_> bp1(3, DT1_(1));
            DenseVector<DT1_> bp2(3, DT1_(1));
            DenseVector<DT1_> bm1(3, DT1_(1));
            DenseVector<DT1_> bm2(3, DT1_(1));

            diag[0] = DT1_(7);
            diag[1] = DT1_(6);
            diag[2] = DT1_(5);
            bp1[0] = DT1_(-2);
            bp1[1] = DT1_(2);
            bp1[2] = DT1_(0);
            bp2[0] = DT1_(0);
            bp2[1] = DT1_(0);
            bp2[2] = DT1_(0);
            bm1[0] = DT1_(0);
            bm1[1] = DT1_(-2);
            bm1[2] = DT1_(2);
            bm2[0] = DT1_(0);
            bm2[1] = DT1_(0);
            bm2[2] = DT1_(0);

            BandedMatrix<DT1_> A(3);
            A.insert_band(0, diag);
            A.insert_band(1, bp1);
            A.insert_band(2, bp2);
            A.insert_band(-1, bm1);
            A.insert_band(-2, bm2);

            DenseVector<DT1_> b(3, DT1_(1));
            b[0] = DT1_(3);
            b[1] = DT1_(3);
            b[2] = DT1_(0);


            std::cout<<"A:"<<A<<endl;
            std::cout<<"b:"<<b<<endl;
            DenseVector<DT1_> result(IterativeRefinement<CG, tags::CPU>::value(A,b,double(std::numeric_limits<double>::epsilon()), double(std::numeric_limits<double>::epsilon())));
            DT1_ x_n = Norm< vnt_l_two, false, Tag_>::value(result);
            DenseVector<DT1_> x_analytical(3, DT1_(0));
            cout<<"RESULT(v1):"<<result<<endl;

            x_analytical[0] = DT1_(2./3.);
            x_analytical[1] = DT1_(5./6.);
            x_analytical[2] = DT1_(-1./3.);
            DT1_ x_analytical_n = Norm< vnt_l_two, false, Tag_>::value(x_analytical);
            TEST_CHECK_EQUAL_WITHIN_EPS(x_analytical_n, x_n , double(0.001));

        }
};

template <typename Tag_, typename DT1_>
class IterativeRefinementTestBandedJacobi:
    public BaseTest
{
    public:
        IterativeRefinementTestBandedJacobi(const std::string & tag) :
            BaseTest("Iterative Refinement with JAC solver test (Banded system)<" + tag + ">")
        {
        }

        virtual void run() const
        {

            DenseVector<DT1_> diag(3, DT1_(1));
            DenseVector<DT1_> bp1(3, DT1_(1));
            DenseVector<DT1_> bp2(3, DT1_(1));
            DenseVector<DT1_> bm1(3, DT1_(1));
            DenseVector<DT1_> bm2(3, DT1_(1));

            diag[0] = DT1_(7);
            diag[1] = DT1_(6);
            diag[2] = DT1_(5);
            bp1[0] = DT1_(-2);
            bp1[1] = DT1_(2);
            bp1[2] = DT1_(0);
            bp2[0] = DT1_(0);
            bp2[1] = DT1_(0);
            bp2[2] = DT1_(0);
            bm1[0] = DT1_(0);
            bm1[1] = DT1_(-2);
            bm1[2] = DT1_(2);
            bm2[0] = DT1_(0);
            bm2[1] = DT1_(0);
            bm2[2] = DT1_(0);

            BandedMatrix<DT1_> A(3);
            A.insert_band(0, diag);
            A.insert_band(1, bp1);
            A.insert_band(2, bp2);
            A.insert_band(-1, bm1);
            A.insert_band(-2, bm2);

            DenseVector<DT1_> b(3, DT1_(1));
            b[0] = DT1_(3);
            b[1] = DT1_(3);
            b[2] = DT1_(0);


            std::cout<<"A:"<<A<<endl;
            std::cout<<"b:"<<b<<endl;
            DenseVector<DT1_> result(IterativeRefinement<JAC, tags::CPU>::value(A,b,double(std::numeric_limits<double>::epsilon()), double(std::numeric_limits<double>::epsilon())));
            DT1_ x_n = Norm< vnt_l_two, false, Tag_>::value(result);
            DenseVector<DT1_> x_analytical(3, DT1_(0));
            cout<<"RESULT(v1):"<<result<<endl;

            x_analytical[0] = DT1_(2./3.);
            x_analytical[1] = DT1_(5./6.);
            x_analytical[2] = DT1_(-1./3.);
            DT1_ x_analytical_n = Norm< vnt_l_two, false, Tag_>::value(x_analytical);
            TEST_CHECK_EQUAL_WITHIN_EPS(x_analytical_n, x_n , double(0.001));

        }
};

template <typename Tag_, typename DT1_>
class IterativeRefinementTestDensePCGJAC:
    public BaseTest
{
    public:
        IterativeRefinementTestDensePCGJAC(const std::string & tag) :
            BaseTest("Iterative Refinement with PCG JAC solver test (Dense system)<" + tag + ">")
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
            DenseVector<DT1_> result(IterativeRefinement<PCG::JAC, tags::CPU>::value(A,b,double(std::numeric_limits<double>::epsilon()), double(std::numeric_limits<double>::epsilon())));
            DT1_ x_n = Norm< vnt_l_two, false, Tag_>::value(result);
            DenseVector<DT1_> x_analytical(3, DT1_(0));
            cout<<"RESULT(v1):"<<result<<endl;

            x_analytical[0] = DT1_(2./3.);
            x_analytical[1] = DT1_(5./6.);
            x_analytical[2] = DT1_(-1./3.);
            DT1_ x_analytical_n = Norm< vnt_l_two, false, Tag_>::value(x_analytical);
            TEST_CHECK_EQUAL_WITHIN_EPS(x_analytical_n, x_n , double(0.001));

        }
};

template <typename Tag_, typename DT1_>
class IterativeRefinementTestBandedPCGJAC:
    public BaseTest
{
    public:
        IterativeRefinementTestBandedPCGJAC(const std::string & tag) :
            BaseTest("Iterative Refinement with PCG (Jacobi) solver test (Banded system)<" + tag + ">")
        {
        }

        virtual void run() const
        {

            DenseVector<DT1_> diag(3, DT1_(1));
            DenseVector<DT1_> bp1(3, DT1_(1));
            DenseVector<DT1_> bp2(3, DT1_(1));
            DenseVector<DT1_> bm1(3, DT1_(1));
            DenseVector<DT1_> bm2(3, DT1_(1));

            diag[0] = DT1_(7);
            diag[1] = DT1_(6);
            diag[2] = DT1_(5);
            bp1[0] = DT1_(-2);
            bp1[1] = DT1_(2);
            bp1[2] = DT1_(0);
            bp2[0] = DT1_(0);
            bp2[1] = DT1_(0);
            bp2[2] = DT1_(0);
            bm1[0] = DT1_(0);
            bm1[1] = DT1_(-2);
            bm1[2] = DT1_(2);
            bm2[0] = DT1_(0);
            bm2[1] = DT1_(0);
            bm2[2] = DT1_(0);

            BandedMatrix<DT1_> A(3);
            A.insert_band(0, diag);
            A.insert_band(1, bp1);
            A.insert_band(2, bp2);
            A.insert_band(-1, bm1);
            A.insert_band(-2, bm2);

            DenseVector<DT1_> b(3, DT1_(1));
            b[0] = DT1_(3);
            b[1] = DT1_(3);
            b[2] = DT1_(0);


            std::cout<<"A:"<<A<<endl;
            std::cout<<"b:"<<b<<endl;
            DenseVector<DT1_> result(IterativeRefinement<PCG::JAC, tags::CPU>::value(A,b,double(std::numeric_limits<double>::epsilon()), double(std::numeric_limits<double>::epsilon())));
            DT1_ x_n = Norm< vnt_l_two, false, Tag_>::value(result);
            DenseVector<DT1_> x_analytical(3, DT1_(0));
            cout<<"RESULT(v1):"<<result<<endl;

            x_analytical[0] = DT1_(2./3.);
            x_analytical[1] = DT1_(5./6.);
            x_analytical[2] = DT1_(-1./3.);
            DT1_ x_analytical_n = Norm< vnt_l_two, false, Tag_>::value(x_analytical);
            TEST_CHECK_EQUAL_WITHIN_EPS(x_analytical_n, x_n , double(0.001));

        }
};
template <typename Tag_, typename DT1_>
class IterativeRefinementTestDensePCGJAC_big:
    public BaseTest
{
    public:
        IterativeRefinementTestDensePCGJAC_big(const std::string & tag) :
            BaseTest("Iterative Refinement with PCG (Jacobi) solver BIG!!!! test (dense system)<" + tag + ">")
        {
        }

        virtual void run() const
        {

            unsigned long size = 1000;
            DenseMatrix<DT1_> A(size,size, DT1_(1));
            DenseVector<DT1_> b(size, DT1_(1));

            //std::cout<<"A:"<<A<<endl;
            //std::cout<<"b:"<<b<<endl;
            DenseVector<DT1_> result(IterativeRefinement<PCG::JAC, tags::CPU>::value(A,b,double(std::numeric_limits<double>::epsilon()), double(std::numeric_limits<double>::epsilon())));
            DT1_ x_n = Norm< vnt_l_two, false, Tag_>::value(result);
            DenseVector<DT1_> x_analytical(size, DT1_(0.001));
            //cout<<"RESULT(v1):"<<result<<endl;

            DT1_ x_analytical_n = Norm< vnt_l_two, false, Tag_>::value(x_analytical);
            TEST_CHECK_EQUAL_WITHIN_EPS(x_analytical_n, x_n , double(0.001));

        }
};

template <typename Tag_, typename DT1_>
class IterativeRefinementTestDenseCG_big:
    public BaseTest
{
    public:
        IterativeRefinementTestDenseCG_big(const std::string & tag) :
            BaseTest("Iterative Refinement with CG  solver BIG!!!! test (dense system)<" + tag + ">")
        {
        }

        virtual void run() const
        {

            unsigned long size = 1000;
            DenseMatrix<DT1_> A(size,size, DT1_(1));
            DenseVector<DT1_> b(size, DT1_(1));

            //std::cout<<"A:"<<A<<endl;
            //std::cout<<"b:"<<b<<endl;
            DenseVector<DT1_> result(IterativeRefinement<CG, tags::CPU>::value(A,b,double(std::numeric_limits<double>::epsilon()), double(std::numeric_limits<double>::epsilon())));
            DT1_ x_n = Norm< vnt_l_two, false, Tag_>::value(result);
            DenseVector<DT1_> x_analytical(size, DT1_(0.001));
            //cout<<"RESULT(v1):"<<result<<endl;

            DT1_ x_analytical_n = Norm< vnt_l_two, false, Tag_>::value(x_analytical);
            TEST_CHECK_EQUAL_WITHIN_EPS(x_analytical_n, x_n , double(0.001));

        }
};
template <typename Tag_, typename DT1_>
class IterativeRefinementTestSparseCG:
    public BaseTest
{
    public:
        IterativeRefinementTestSparseCG(const std::string & tag) :
            BaseTest("Iterative Refinement with CG solver test (sparse system)<" + tag + ">")
        {
        }

        virtual void run() const
        {
            SparseMatrix<DT1_> A(3, 3);
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
            DenseVector<DT1_> result(IterativeRefinement<CG, tags::CPU>::value(A,b,double(std::numeric_limits<double>::epsilon()), double(std::numeric_limits<double>::epsilon())));
            DT1_ x_n = Norm< vnt_l_two, false, Tag_>::value(result);
            DenseVector<DT1_> x_analytical(3, DT1_(0));
            cout<<"RESULT(v1):"<<result<<endl;

            x_analytical[0] = DT1_(2./3.);
            x_analytical[1] = DT1_(5./6.);
            x_analytical[2] = DT1_(-1./3.);
            DT1_ x_analytical_n = Norm< vnt_l_two, false, Tag_>::value(x_analytical);
            TEST_CHECK_EQUAL_WITHIN_EPS(x_analytical_n, x_n , double(0.001));

        }
};

template <typename Tag_, typename DT1_>
class IterativeRefinementTestSparsePCGJAC:
    public BaseTest
{
    public:
        IterativeRefinementTestSparsePCGJAC(const std::string & tag) :
            BaseTest("Iterative Refinement with PCG (JAC) solver test (sparse system)<" + tag + ">")
        {
        }

        virtual void run() const
        {
            SparseMatrix<DT1_> A(3, 3);
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
            DenseVector<DT1_> result(IterativeRefinement<PCG::JAC, tags::CPU>::value(A,b,double(std::numeric_limits<double>::epsilon()), double(std::numeric_limits<double>::epsilon())));
            DT1_ x_n = Norm< vnt_l_two, false, Tag_>::value(result);
            DenseVector<DT1_> x_analytical(3, DT1_(0));
            cout<<"RESULT(v1):"<<result<<endl;

            x_analytical[0] = DT1_(2./3.);
            x_analytical[1] = DT1_(5./6.);
            x_analytical[2] = DT1_(-1./3.);
            DT1_ x_analytical_n = Norm< vnt_l_two, false, Tag_>::value(x_analytical);
            TEST_CHECK_EQUAL_WITHIN_EPS(x_analytical_n, x_n , double(0.001));

        }
};
template <typename Tag_, typename DT1_>
class IterativeRefinementTestSparseJAC:
    public BaseTest
{
    public:
        IterativeRefinementTestSparseJAC(const std::string & tag) :
            BaseTest("Iterative Refinement with Jacobi solver test (sparse system)<" + tag + ">")
        {
        }

        virtual void run() const
        {
            SparseMatrix<DT1_> A(3, 3);
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
            DenseVector<DT1_> result(IterativeRefinement<JAC, tags::CPU>::value(A,b,double(std::numeric_limits<double>::epsilon()), double(std::numeric_limits<double>::epsilon())));
            DT1_ x_n = Norm< vnt_l_two, false, Tag_>::value(result);
            DenseVector<DT1_> x_analytical(3, DT1_(0));
            cout<<"RESULT(v1):"<<result<<endl;

            x_analytical[0] = DT1_(2./3.);
            x_analytical[1] = DT1_(5./6.);
            x_analytical[2] = DT1_(-1./3.);
            DT1_ x_analytical_n = Norm< vnt_l_two, false, Tag_>::value(x_analytical);
            TEST_CHECK_EQUAL_WITHIN_EPS(x_analytical_n, x_n , double(0.001));

        }
};
IterativeRefinementTestSparseCG<tags::CPU, double> iterref_test_double_sparseCG("double");
IterativeRefinementTestSparsePCGJAC<tags::CPU, double> iterref_test_double_sparsePCGJAC("double");
IterativeRefinementTestSparseJAC<tags::CPU, double> iterref_test_double_sparseJAC("double");
IterativeRefinementTestDenseCG<tags::CPU, double> iterref_test_double_denseCG("double");
IterativeRefinementTestDenseJAC<tags::CPU, double> iterref_test_double_denseJAC("double");
IterativeRefinementTestBandedCG<tags::CPU, double> iterref_test_double_bandedCG("double");
IterativeRefinementTestBandedJacobi<tags::CPU, double> iterref_test_double_bandedJAC("double");
IterativeRefinementTestDensePCGJAC<tags::CPU, double> iterref_test_double_densePCGJAC("double");
IterativeRefinementTestBandedPCGJAC<tags::CPU, double> iterref_test_double_bandedPCGJAC("double");
IterativeRefinementTestDensePCGJAC_big<tags::CPU, double> iterref_test_double_densePCGJAC_big("double");
IterativeRefinementTestDenseCG_big<tags::CPU, double> iterref_test_double_denseCG_big("double");

IterativeRefinementTestSparseCG<tags::CPU, float> iterref_test_float_sparseCG("float");
IterativeRefinementTestSparsePCGJAC<tags::CPU, float> iterref_test_float_sparsePCGJAC("float");
IterativeRefinementTestSparseJAC<tags::CPU, float> iterref_test_float_sparseJAC("float");
IterativeRefinementTestDenseCG<tags::CPU, float> iterref_test_float_denseCG("float");
IterativeRefinementTestDenseJAC<tags::CPU, float> iterref_test_float_denseJAC("float");
IterativeRefinementTestBandedCG<tags::CPU, float> iterref_test_float_bandedCG("float");
IterativeRefinementTestBandedJacobi<tags::CPU, float> iterref_test_float_bandedJAC("float");
IterativeRefinementTestDensePCGJAC<tags::CPU, float> iterref_test_float_densePCGJAC("float");
IterativeRefinementTestBandedPCGJAC<tags::CPU, float> iterref_test_float_bandedPCGJAC("float");
IterativeRefinementTestDensePCGJAC_big<tags::CPU, float> iterref_test_float_densePCGJAC_big("float");
IterativeRefinementTestDenseCG_big<tags::CPU, float> iterref_test_float_denseCG_big("float");

