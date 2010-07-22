/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Dirk Ribbrock <dirk.ribbrock@math.uni-dortmund.de
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

#include <honei/math/conjugate_gradients.hh>
#include <honei/math/matrix_io.hh>
#include <honei/math/vector_io.hh>
#include <unittest/unittest.hh>
#include <honei/util/stringify.hh>
#include <honei/math/ludecomposition.hh>
#include <honei/la/product.hh>
#include <limits>

using namespace honei;
using namespace tests;

template<typename DT_, typename Tag_>
class LUQuickTest:
    public QuickTest
{
    public:
        LUQuickTest(const std::string & tag) :
            QuickTest("lu_decomp_quick_test_"+tag)
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            DenseMatrix<DT_> a(3, 3, DT_(1));
            DenseVector<DT_> b(3, DT_(1));
            a[0][0] = DT_(7);
            a[0][1] = DT_(-2);
            a[0][2] = DT_(0);
            a[1][0] = DT_(-2);
            a[1][1] = DT_(6);
            a[1][2] = DT_(2);
            a[2][0] = DT_(0);
            a[2][1] = DT_(2);
            a[2][2] = DT_(5);

            b[0] = DT_(3);
            b[1] = DT_(3);
            b[2] = DT_(0);


            DenseVector<DT_> x_analytical(3);
            x_analytical[0] = DT_(2./3.);
            x_analytical[1] = DT_(5./6.);
            x_analytical[2] = DT_(-1./3.);

            DenseVector<DT_> result(3);
            LUDecomposition<Tag_>::value(a, b, result);
            TEST_CHECK_EQUAL(result, x_analytical);
        }
};
LUQuickTest<double, tags::CPU> lu_quick_test_double("double");

template<typename DT_, typename Tag_>
class LUQuickTest2:
    public QuickTest
{
    public:
        LUQuickTest2(const std::string & tag) :
            QuickTest("lu_decomp_quick_test_2_"+tag)
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
                DenseMatrix<DT_> a(4, 4);
                a(0, 0) = 2;
                a(0, 1) = 1;
                a(0, 2) = 1;
                a(0, 3) = 0;
                a(1, 0) = 4;
                a(1, 1) = 3;
                a(1, 2) = 3;
                a(1, 3) = 1;
                a(2, 0) = 8;
                a(2, 1) = 7;
                a(2, 2) = 9;
                a(2, 3) = 5;
                a(3, 0) = 6;
                a(3, 1) = 7;
                a(3, 2) = 9;
                a(3, 3) = 8;

                DenseVector<DT_> x(4);
                for (unsigned long i(0) ; i < x.size() ; ++i)
                {
                    x.elements()[i] = DT_(i * 1.5);
                }

                DenseVector<DT_> b(Product<Tag_>::value(a, x));

                DenseVector<DT_> result(4);
                LUDecomposition<Tag_>::value(a, b, result);
                for (unsigned long i(0) ; i < x.size() ; ++i)
                    TEST_CHECK_EQUAL_WITHIN_EPS(result[i], x[i], std::numeric_limits<double>::epsilon() * 1e3);
        }
};
LUQuickTest2<double, tags::CPU> lu_quick_2_test_double("double");

template<typename DT_, typename Tag_>
class PLUQuickTest:
    public QuickTest
{
    public:
        PLUQuickTest(const std::string & tag) :
            QuickTest("plu_decomp_quick_test_"+tag)
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            DenseMatrix<DT_> a(2, 2);
            DenseVector<DT_> b(2);
            a[0][0] = DT_(0);
            a[0][1] = DT_(1);
            a[1][0] = DT_(1);
            a[1][1] = DT_(1);

            b[0] = DT_(1);
            b[1] = DT_(0);


            DenseVector<DT_> x_analytical(2);
            x_analytical[0] = DT_(-1);
            x_analytical[1] = DT_(1);

            DenseVector<DT_> result(2);
            LUDecomposition<Tag_>::value(a, b, result);
            TEST_CHECK_EQUAL(result, x_analytical);
        }
};
PLUQuickTest<double, tags::CPU> plu_quick_test_double("double");

template <typename Tag_, typename DT1_>
class LUTestDenseELL:
    public BaseTest
{
    private:
        std::string _m_f, _v_f, _i_f;
    public:
        LUTestDenseELL(const std::string & tag,
                std::string m_file,
                std::string v_file,
                std::string init_file) :
            BaseTest("LU Test (Dense ELL system)<" + tag + ">")
        {
            register_tag(Tag_::name);
            _m_f = m_file;
            _v_f = v_file;
            _i_f = init_file;
        }

        virtual void run() const
        {

            std::string filename(HONEI_SOURCEDIR);
            filename += "/honei/math/testdata/";
            filename += _m_f;
            unsigned long non_zeros(MatrixIO<io_formats::M>::get_non_zeros(filename));
            unsigned long rows, columns, ax, bx;
            DenseVector<unsigned long> r(non_zeros);
            DenseVector<unsigned long> c(non_zeros);
            DenseVector<DT1_> data(non_zeros);

            MatrixIO<io_formats::M>::read_matrix(filename, r, c, data);
            MatrixIO<io_formats::M>::get_sizes(filename, rows, columns, ax, bx);
            SparseMatrix<DT1_> tsmatrix(rows, columns, r, c, data);
            SparseMatrixELL<DT1_> smatrix(tsmatrix);

            std::string filename_2(HONEI_SOURCEDIR);
            filename_2 += "/honei/math/testdata/";
            filename_2 += _v_f;
            DenseVector<DT1_> rhs(rows, DT1_(0));
            VectorIO<io_formats::EXP>::read_vector(filename_2, rhs);

            DenseVector<DT1_> diag_inverted(rows, DT1_(0));
            for(unsigned long i(0) ; i < data.size() ; ++i)
            {
                if(r[i] == c[i])
                    diag_inverted[r[i]] = DT1_(1)/data[i];
            }
            for(unsigned long i(0) ; i < data.size() ; ++i)
            {
                if(r[i] == c[i])
                     data[i] = DT1_(0);
            }
            SparseMatrix<DT1_> tdifference(rows, columns, r, c, data);
            SparseMatrixELL<DT1_> difference(tdifference);

            DenseVector<DT1_> init(rhs.size(), DT1_(0));
            std::string filename_4(HONEI_SOURCEDIR);
            filename_4 += "/honei/math/testdata/";
            filename_4 += _i_f;
            VectorIO<io_formats::EXP>::read_vector(filename_4, init);
            DenseVector<DT1_> result(init.copy());

            unsigned long used_iters(0);
            ConjugateGradients<Tag_, JAC>::value(smatrix, rhs, result, diag_inverted, 10000ul, used_iters, DT1_(1e-8));
            //Jacobi<Tag_>::value(smatrix, difference, rhs, result, diag_inverted, 10000ul, used_iters, DT1_(1e-8));


            DenseVector<DT1_> lu_result(rhs.size());
            DenseMatrix<DT1_> dmatrix(tsmatrix);
            LUDecomposition<Tag_>::value(dmatrix, rhs, lu_result);

            DT1_ base_digits(3);
            DT1_ additional_digits(2);

            DT1_ base_eps(1 / pow(10, base_digits));
            DT1_ add_eps(base_eps / pow(10, additional_digits));

            DT1_ m((add_eps - base_eps) / DT1_(4));
            DT1_ b(base_eps - (DT1_(4) * m));

            DT1_ eps(m * sizeof(DT1_) + b);
            eps *= DT1_(3);
            for (unsigned long i(0) ; i < result.size() ; ++i)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(lu_result[i], result[i], eps);
            }
        }
};
LUTestDenseELL<tags::CPU, double> lu_test_dense_ell_double_2("double", "l2/area51_full_0.m", "l2/area51_rhs_0", "l2/area51_init_0");

template <typename Tag_, typename DT1_>
class LUTestSparseELL:
    public BaseTest
{
    private:
        std::string _m_f, _v_f, _r_f;
    public:
        LUTestSparseELL(const std::string & tag,
                std::string m_file,
                std::string v_file,
                std::string ref_file) :
            BaseTest("LU Test (sparse ELL system)<" + tag + ">")
        {
            register_tag(Tag_::name);
            _m_f = m_file;
            _v_f = v_file;
            _r_f = ref_file;
        }

        virtual void run() const
        {

            std::string filename(HONEI_SOURCEDIR);
            filename += "/honei/math/testdata/";
            filename += _m_f;

            SparseMatrixELL<DT1_> smatrix(MatrixIO<io_formats::ELL>::read_matrix(filename, DT1_(0)));

            std::string filename_2(HONEI_SOURCEDIR);
            filename_2 += "/honei/math/testdata/";
            filename_2 += _v_f;
            DenseVector<DT1_> rhs(smatrix.rows(), DT1_(0));
            VectorIO<io_formats::EXP>::read_vector(filename_2, rhs);

            std::string filename_4(HONEI_SOURCEDIR);
            filename_4 += "/honei/math/testdata/";
            filename_4 += _r_f;
            DenseVector<DT1_> result_ref(rhs.size());
            VectorIO<io_formats::EXP>::read_vector(filename_4, result_ref);

            DenseVector<DT1_> lu_result(rhs.size());
            SparseMatrix<DT1_> tsmatrix(smatrix);
            LUDecomposition<Tag_>::value(tsmatrix, rhs, lu_result);

            DT1_ base_digits(3);
            DT1_ additional_digits(2);

            DT1_ base_eps(1 / pow(10, base_digits));
            DT1_ add_eps(base_eps / pow(10, additional_digits));

            DT1_ m((add_eps - base_eps) / DT1_(4));
            DT1_ b(base_eps - (DT1_(4) * m));

            DT1_ eps(m * sizeof(DT1_) + b);
            eps *= DT1_(3);
            for (unsigned long i(0) ; i < result_ref.size() ; ++i)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(lu_result[i], result_ref[i], eps);
            }
        }
};
LUTestSparseELL<tags::CPU, double> lu_test_sparse_ell_double_1("double", "l2/area51_full_0.ell", "l2/area51_rhs_0", "l2/area51_sol_0");
//LUTestSparseELL<tags::CPU, double> lu_test_sparse_ell_double_2("double", "poisson_advanced/sort_0/A_7.ell", "poisson_advanced/sort_0/rhs_7", "poisson_advanced/sort_0/sol_7");
