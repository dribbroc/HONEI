/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2009 Markus Geveler <apryde@gmx.de>
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

#include <honei/math/matrix_io.hh>
#include <honei/math/vector_io.hh>
#include <honei/math/defect.hh>
#include <unittest/unittest.hh>
#include <honei/util/stringify.hh>
#include <honei/la/product.hh>

using namespace honei;
using namespace tests;
using namespace std;

template<typename DT_, typename Tag_>
class DefectTest:
    public BaseTest
{
    public:
        DefectTest(const std::string & tag) :
            BaseTest("Defect Test " + tag)
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            std::string filename(HONEI_SOURCEDIR);
            filename += "/honei/math/testdata/5pt_10x10.mtx";
            unsigned long non_zeros(0);
            DenseMatrix<DT_> matrix = MatrixIO<io_formats::MTX>::read_matrix(filename, DT_(0), non_zeros);

            DenseVector<DT_> x(matrix.rows());
            DenseVector<DT_> b(matrix.rows(), DT_(1.234));
            DenseVector<DT_> y(matrix.rows());
            for (unsigned long i(0) ; i < x.size() ; ++i)
            {
                x[i] = DT_(i) / 1.234;
            }
            SparseMatrix<DT_> ssmatrix(matrix);
            SparseMatrixELL<DT_> smatrix(ssmatrix);

            DenseVector<DT_> y1(Defect<Tag_>::value(b, smatrix, x));
            DenseVector<DT_> y2(b.size());
            Defect<Tag_>::value(y2, b, smatrix, x);
            DenseVector<DT_> yref(b.copy());
            Difference<tags::CPU>::value(yref ,Product<tags::CPU>::value(matrix, x) );

            y1.lock(lm_read_only);
            y1.unlock(lm_read_only);
            y2.lock(lm_read_only);
            y2.unlock(lm_read_only);
            TEST_CHECK_EQUAL(y1, yref);
            TEST_CHECK_EQUAL(y2, yref);

        }
};
DefectTest<float, tags::CPU> defect_test_float_sparse("float");
DefectTest<double, tags::CPU> defect_test_double_sparse("double");
#ifdef HONEI_SSE
DefectTest<float, tags::CPU::SSE> sse_defect_test_float_sparse("float");
DefectTest<double, tags::CPU::SSE> sse_defect_test_double_sparse("double");
#endif
#ifdef HONEI_CUDA
DefectTest<float, tags::GPU::CUDA> cuda_defect_test_float_sparse("float");
DefectTest<float, tags::GPU::MultiCore::CUDA> mc_cuda_defect_test_float_sparse("float");
#ifdef HONEI_CUDA_DOUBLE
DefectTest<double, tags::GPU::CUDA> cuda_defect_test_double_sparse("double");
DefectTest<double, tags::GPU::MultiCore::CUDA> mc_cuda_defect_test_double_sparse("double");
#endif
#endif

template<typename DT_, typename Tag_>
class DefectRegressionTest:
    public BaseTest
{
    private:
        std::string _m_f, _v_f;
    public:
        DefectRegressionTest(const std::string & tag, std::string m_file, std::string v_file) :
            BaseTest("Defect Regression Test " + tag)
        {
            register_tag(Tag_::name);
            _m_f = m_file;
            _v_f = v_file;
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
            DenseVector<DT_> data(non_zeros);

            MatrixIO<io_formats::M>::read_matrix(filename, r, c, data);
            MatrixIO<io_formats::M>::get_sizes(filename, rows, columns, ax, bx);
            SparseMatrix<DT_> tsmatrix2(rows, columns, r, c, data);
            SparseMatrixELL<DT_> smatrix2(tsmatrix2);

            std::string filename_2(HONEI_SOURCEDIR);
            filename_2 += "/honei/math/testdata/";
            filename_2 += _v_f;
            DenseVector<DT_> rhs(rows, DT_(0));
            VectorIO<io_formats::EXP>::read_vector(filename_2, rhs);

            DenseVector<DT_> x(rows, DT_(0));
            for (unsigned long i(0) ; i < x.size() ; ++i)
                if (i%2 == 0) x[i] = 1;

            DenseVector<DT_> ref_result(Defect<tags::CPU>::value(rhs, smatrix2, x));

            DenseVector<DT_> result(Defect<Tag_>::value(rhs, smatrix2, x));
            DenseVector<DT_> result2(result.size());
            DenseVector<DT_> result3(result.size());
            Defect<Tag_>::value(result2, rhs, smatrix2, x);
            Defect<Tag_>::value(result3, rhs, smatrix2, result2);
            DenseVector<DT_> ref_result3(Defect<tags::CPU>::value(rhs, smatrix2, ref_result));

            result.lock(lm_read_only);
            result2.lock(lm_read_only);
            result3.lock(lm_read_only);
            TEST_CHECK_EQUAL(result, ref_result);
            TEST_CHECK_EQUAL(result2, ref_result);
            for (unsigned long i(0) ; i < x.size() ; ++i)
            TEST_CHECK_EQUAL_WITHIN_EPS(result3[i], ref_result3[i], 1e-3);
            result2.unlock(lm_read_only);
            result.unlock(lm_read_only);
            result3.unlock(lm_read_only);
        }
};
DefectRegressionTest<float, tags::CPU> regression_defect_test_float_sparse("Regression float", "l2/area51_full_0.m", "l2/area51_rhs_0");
DefectRegressionTest<double, tags::CPU> regression_defect_test_double_sparse("Regression double", "l2/area51_full_0.m", "l2/area51_rhs_0");
#ifdef HONEI_SSE
DefectRegressionTest<float, tags::CPU::SSE> sse_regression_defect_test_float_sparse("Regression float", "l2/area51_full_0.m", "l2/area51_rhs_0");
DefectRegressionTest<double, tags::CPU::SSE> sse_regression_defect_test_double_sparse("Regression double", "l2/area51_full_0.m", "l2/area51_rhs_0");
DefectRegressionTest<float, tags::CPU::MultiCore::SSE> mc_sse_regression_defect_test_float_sparse("Regression float", "l2/area51_full_0.m", "l2/area51_rhs_0");
DefectRegressionTest<double, tags::CPU::MultiCore::SSE> mc_sse_regression_defect_test_double_sparse("Regression double", "l2/area51_full_0.m", "l2/area51_rhs_0");
#endif
#ifdef HONEI_CUDA
DefectRegressionTest<float, tags::GPU::CUDA> cuda_regression_defect_test_float_sparse("CUDA Regression float", "l2/area51_full_0.m", "l2/area51_rhs_0");
DefectRegressionTest<float, tags::GPU::MultiCore::CUDA> mc_cuda_regression_defect_test_float_sparse("CUDA Regression float", "l2/area51_full_0.m", "l2/area51_rhs_0");
#ifdef HONEI_CUDA_DOUBLE
DefectRegressionTest<double, tags::GPU::CUDA> cuda_regression_defect_test_double_sparse("CUDA Regression double", "area51_full_0.m", "area51_rhs_0");
DefectRegressionTest<double, tags::GPU::MultiCore::CUDA> mc_cuda_regression_defect_test_double_sparse("CUDA Regression double", "area51_full_0.m", "area51_rhs_0");
#endif
#endif
