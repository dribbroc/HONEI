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

#include <honei/math/jacobi.hh>
#include <honei/math/matrix_io.hh>
#include <unittest/unittest.hh>
#include <honei/util/stringify.hh>
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

            std::cout<<"A:"<<A<<std::endl;
            std::cout<<"b:"<<b<<std::endl;
            DenseVector<DT1_> result = Jacobi<Tag_>::value(A,b,12ul);
            DT1_ x_n = Norm< vnt_l_two, false, Tag_>::value(result);
            DenseVector<DT1_> x_analytical(3, DT1_(0));
            std::cout << "RESULT(v1):" << result << std::endl;

            x_analytical[0] = DT1_(2./3.);
            x_analytical[1] = DT1_(5./6.);
            x_analytical[2] = DT1_(-1./3.);
            DT1_ x_analytical_n = Norm< vnt_l_two, false, Tag_>::value(x_analytical);
            TEST_CHECK_EQUAL_WITHIN_EPS(x_analytical_n, x_n , double(0.001));

            DenseVector<DT1_> result_2 = Jacobi<Tag_>::value(A,b,double(std::numeric_limits<double>::epsilon()));
            std::cout << "RESULT(v2):" << result_2 << std::endl;

            DT1_ x_n_2 = Norm< vnt_l_two, false, Tag_>::value(result_2);
            TEST_CHECK_EQUAL_WITHIN_EPS(x_analytical_n, x_n_2 , double(0.001));

        }
};

template <typename Tag_, typename DT1_>
class JacobiTestBanded:
    public BaseTest
{
    public:
        JacobiTestBanded(const std::string & tag) :
            BaseTest("Jacobi solver test (Banded system)<" + tag + ">")
        {
            register_tag(Tag_::name);
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

            std::cout << "A:" << A << std::endl;
            std::cout << "b:" << b << std::endl;

            DenseVector<DT1_> result = Jacobi<Tag_>::value(A,b,10ul);
            DT1_ x_n = Norm< vnt_l_two, false, Tag_>::value(result);
            DenseVector<DT1_> x_analytical(3, DT1_(1));
            x_analytical[0] = DT1_(2./3.);
            x_analytical[1] = DT1_(5./6.);
            x_analytical[2] = DT1_(-1./3.);
            DT1_ x_analytical_n = Norm< vnt_l_two, false, Tag_>::value(x_analytical);

            std::cout << "RESULT(v1):" << result << std::endl;
            TEST_CHECK_EQUAL_WITHIN_EPS(x_analytical_n, x_n , double(0.1));

            DenseVector<DT1_> result_2 = Jacobi<Tag_>::value(A,b,double(std::numeric_limits<double>::epsilon()));
            DT1_ x_n_2 = Norm< vnt_l_two, false, Tag_>::value(result_2);
            TEST_CHECK_EQUAL_WITHIN_EPS(x_analytical_n, x_n_2 , double(0.1));

            std::cout << "RESULT(v2):" << result_2 << std::endl;

        }
};
template <typename Tag_, typename DT1_>
class JacobiTestSparse:
    public BaseTest
{
    public:
        JacobiTestSparse(const std::string & tag) :
            BaseTest("Jacobi solver test (sparse system)<" + tag + ">")
        {
            register_tag(Tag_::name);
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

            std::cout << "A:" << A << std::endl;
            std::cout << "b:" << b << std::endl;
            DenseVector<DT1_> result = Jacobi<Tag_>::value(A,b, 12ul);
            DT1_ x_n = Norm< vnt_l_two, false, Tag_>::value(result);
            DenseVector<DT1_> x_analytical(3, DT1_(0));
            std::cout << "RESULT(v1):" << result << std::endl;

            x_analytical[0] = DT1_(2./3.);
            x_analytical[1] = DT1_(5./6.);
            x_analytical[2] = DT1_(-1./3.);
            DT1_ x_analytical_n = Norm< vnt_l_two, false, Tag_>::value(x_analytical);
            TEST_CHECK_EQUAL_WITHIN_EPS(x_analytical_n, x_n , double(0.001));

            DenseVector<DT1_> result_2 = Jacobi<Tag_>::value(A,b,double(std::numeric_limits<double>::epsilon()));
            std::cout << "RESULT(v2):" << result_2 << std::endl;

            DT1_ x_n_2 = Norm< vnt_l_two, false, Tag_>::value(result_2);
            TEST_CHECK_EQUAL_WITHIN_EPS(x_analytical_n, x_n_2 , double(0.001));

        }
};
JacobiTestSparse<tags::CPU, float> jacobi_test_float_sparse("float");
JacobiTestSparse<tags::CPU, double> jacobi_test_double_sparse("double");

JacobiTestDense<tags::CPU, float> jacobi_test_float_dense("float");
JacobiTestDense<tags::CPU, double> jacobi_test_double_dense("double");
JacobiTestBanded<tags::CPU, float> jacobi_test_float_banded("float");
JacobiTestBanded<tags::CPU, double> jacobi_test_double_banded("double");


#ifdef HONEI_SSE
JacobiTestDense<tags::CPU::SSE, float> sse_jacobi_test_float_dense("SSE float");
JacobiTestDense<tags::CPU::SSE, double> sse_jacobi_test_double_dense("SSE double");
JacobiTestBanded<tags::CPU::SSE, float> sse_jacobi_test_float_banded("SSE float");
JacobiTestBanded<tags::CPU::SSE, double> sse_jacobi_test_double_banded("SSE double");

#endif

template <typename Tag_, typename DT1_>
class JacobiTestSparseELL:
    public BaseTest
{
    public:
        JacobiTestSparseELL(const std::string & tag) :
            BaseTest("Jacobi solver test (sparse ELL system)<" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {

            std::string filename(HONEI_SOURCEDIR);
            filename += "/honei/math/testdata/5pt_10x10.mtx";
            unsigned long non_zeros(MatrixIO<io_formats::MTX>::get_non_zeros(filename));
            unsigned long rows, columns, ax, bx;
            DenseVector<unsigned long> r(non_zeros);
            DenseVector<unsigned long> c(non_zeros);
            DenseVector<DT1_> data(non_zeros);

            MatrixIO<io_formats::MTX>::read_matrix(filename, r, c, data);
            MatrixIO<io_formats::MTX>::get_sizes(filename, rows, columns, ax, bx);
            SparseMatrixELL<DT1_> smatrix2(rows, columns, r, c, data);

            DenseVector<DT1_> x(rows, DT1_(1.2345));

            DenseVector<DT1_> rhs(rows, DT1_(0));

            Product<Tag_>::value(rhs, smatrix2, x);

            DenseVector<DT1_> initial_guess(rows, DT1_(1));
            DenseVector<DT1_> diag_inverted(rows, DT1_(0));
            for(unsigned long i(0) ; i < data.size() ; ++i)
            {
                if(r[i] == c[i])
                    diag_inverted[r[i]] = DT1_(1)/data[i];
            }

            //Initial defect (for convergence purposes later):
            DenseVector<DT1_> initial_defect(Defect<Tag_>::value(rhs, smatrix2, x));

            DenseVector<DT1_> result(Jacobi<Tag_>::value(initial_guess, smatrix2, rhs, 1000ul, DT1_(0.7), diag_inverted));

            result.lock(lm_read_only);
            x.lock(lm_read_only);
            for(unsigned long i(0) ; i < rhs.size() ; ++i)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(result[i], x[i], std::numeric_limits<DT1_>::epsilon());
            }
            result.unlock(lm_read_only);
            x.unlock(lm_read_only);
        }
};
JacobiTestSparseELL<tags::CPU, float> jacobi_test_float_sparse_ell("float");
JacobiTestSparseELL<tags::CPU, double> jacobi_test_double_sparse_ell("double");
#ifdef HONEI_SSE
JacobiTestSparseELL<tags::CPU::SSE, float> sse_jacobi_test_float_sparse_ell("float");
JacobiTestSparseELL<tags::CPU::SSE, double> sse_jacobi_test_double_sparse_ell("double");
#endif
#ifdef HONEI_CUDA
JacobiTestSparseELL<tags::GPU::CUDA, float> cuda_jacobi_test_float_sparse_ell("float");
#ifdef HONEI_CUDA_DOUBLE
JacobiTestSparseELL<tags::GPU::CUDA, double> cuda_jacobi_test_double_sparse_ell("double");
#endif
#endif
