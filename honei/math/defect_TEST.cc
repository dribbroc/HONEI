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
#include <honei/util/unittest.hh>
#include <honei/util/stringify.hh>
#include <honei/util/configuration.hh>
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
            unsigned long old_threads = Configuration::instance()->get_value("ell::threads", 1);
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
            for (unsigned long threads(1) ; threads <= 16 ; threads *= 2)
            {
                Configuration::instance()->set_value("ell::threads", threads);
                SparseMatrixELL<DT_> smatrix(ssmatrix);

                DenseVector<DT_> y1(Defect<Tag_>::value(b, smatrix, x));
                DenseVector<DT_> y2(b.size());
                Defect<Tag_>::value(y2, b, smatrix, x);
                DenseVector<DT_> yref(b.copy());
                Difference<tags::CPU>::value(yref ,Product<tags::CPU>::value(matrix, x));

                y1.lock(lm_read_only);
                y1.unlock(lm_read_only);
                y2.lock(lm_read_only);
                y2.unlock(lm_read_only);
                for (unsigned long i(0) ; i < x.size() ; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(y1[i], yref[i], 1e-3);
                    TEST_CHECK_EQUAL_WITHIN_EPS(y2[i], yref[i], 1e-3);
                }
            }

            Configuration::instance()->set_value("ell::threads", old_threads);
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
#ifdef HONEI_OPENCL
DefectTest<float, tags::OpenCL::CPU> ocl_cpu_defect_test_float_sparse("float");
DefectTest<double, tags::OpenCL::CPU> ocl_cpu_defect_test_double_sparse("double");
DefectTest<float, tags::OpenCL::GPU> ocl_gpu_defect_test_float_sparse("float");
#ifdef HONEI_CUDA_DOUBLE
DefectTest<double, tags::OpenCL::GPU> ocl_gpu_defect_test_double_sparse("double");
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
            SparseMatrix<DT_> tsmatrix2(MatrixIO<io_formats::M>::read_matrix(filename, DT_(0)));
            SparseMatrixELL<DT_> smatrix2(tsmatrix2);

            std::string filename_2(HONEI_SOURCEDIR);
            filename_2 += "/honei/math/testdata/";
            filename_2 += _v_f;
            DenseVector<DT_> rhs(VectorIO<io_formats::EXP>::read_vector(filename_2, DT_(0)));

            DenseVector<DT_> x(tsmatrix2.rows(), DT_(0));
            for (unsigned long i(0) ; i < x.size() ; ++i)
                if (i%2 == 0) x[i] = 1;

            DenseVector<DT_> ref_result(Defect<tags::CPU>::value(rhs, smatrix2, x));

            DenseVector<DT_> result(Defect<Tag_>::value(rhs, smatrix2, x));
            DenseVector<DT_> result2(result.size());
            DenseVector<DT_> result3(result.size());
            DenseVector<DT_> result4(result.size());
            Defect<Tag_>::value(result2, rhs, smatrix2, x);
            Defect<Tag_>::value(result3, rhs, smatrix2, result2);
            Defect<Tag_>::value(result4, rhs, smatrix2, result3);
            DenseVector<DT_> ref_result3(Defect<tags::CPU>::value(rhs, smatrix2, ref_result));
            DenseVector<DT_> ref_result4(Defect<tags::CPU>::value(rhs, smatrix2, ref_result3));

            result.lock(lm_read_only);
            result2.lock(lm_read_only);
            result3.lock(lm_read_only);
            result4.lock(lm_read_only);
            for (unsigned long i(0) ; i < x.size() ; ++i)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(result[i], ref_result[i], 1e-3);
                TEST_CHECK_EQUAL_WITHIN_EPS(result2[i], ref_result[i], 1e-3);
            }
            for (unsigned long i(0) ; i < x.size() ; ++i)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(result3[i], ref_result3[i], 1e-3);
                TEST_CHECK_EQUAL_WITHIN_EPS(result4[i], ref_result4[i], 1e-3);
            }
            result2.unlock(lm_read_only);
            result.unlock(lm_read_only);
            result3.unlock(lm_read_only);
            result4.unlock(lm_read_only);
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
#ifdef HONEI_OPENCL
DefectRegressionTest<float, tags::OpenCL::CPU> ocl_cpu_regression_defect_test_float_sparse("Regression float", "l2/area51_full_0.m", "l2/area51_rhs_0");
DefectRegressionTest<double, tags::OpenCL::CPU> ocl_cpu_regression_defect_test_double_sparse("Regression double", "l2/area51_full_0.m", "l2/area51_rhs_0");
DefectRegressionTest<float, tags::OpenCL::GPU> ocl_gpu_regression_defect_test_float_sparse("Regression float", "l2/area51_full_0.m", "l2/area51_rhs_0");
#ifdef HONEI_CUDA_DOUBLE
DefectRegressionTest<double, tags::OpenCL::GPU> ocl_gpu_regression_defect_test_double_sparse("Regression double", "l2/area51_full_0.m", "l2/area51_rhs_0");
#endif
#endif

template <typename Tag_, typename DataType_>
class Q1MatrixDenseVectorDefectTest :
    public BaseTest
{
    public:
        Q1MatrixDenseVectorDefectTest(const std::string & type) :
            BaseTest("q1_matrix_dense_vector_defect_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size;
            for (unsigned long level(1) ; level <= 10 ; ++level)
            {
                size  = (unsigned long)pow((pow(double(2), double(level)) + 1), 2);
                unsigned long num_limit(311); //value of used elements will be <= num_limit* size

                DenseVector<DataType_> dv1(size, DataType_(1));
                DenseVector<DataType_> dv2(size, DataType_(1));
                DenseVector<DataType_> dv3(size, DataType_(1));
                DenseVector<DataType_> dv4(size, DataType_(1));
                DenseVector<DataType_> b(size, DataType_(0.4711));

                for (unsigned long i(0); i < size; ++i)
                {
                    (dv1)[i]= DataType_((i + 2) % num_limit);
                }
                for (unsigned long i(0) ; i < size ; ++i)
                {
                    (dv2)[i]= DataType_((i + 1) % num_limit);
                }
                for (unsigned long i(0) ; i < size ; ++i)
                {
                    (dv3)[i]= DataType_((i + 25) % num_limit);
                }
                for (unsigned long i(0) ; i < size ; ++i)
                {
                    (dv4)[i]= DataType_((i + 7) % num_limit);
                }

                BandedMatrixQx<Q1Type, DataType_> bm1(size, dv3, dv2, dv4, dv2, dv3, dv4, dv4, dv3, dv2);

                DenseVector<DataType_> prod2(b.size());
                Defect<Tag_>::value(prod2, b, bm1, dv1);
#ifdef HONEI_SSE
            DenseVector<DataType_> y_ref(b.copy());
            Difference<tags::CPU::SSE>::value(y_ref ,Product<tags::CPU::SSE>::value(bm1, dv1));
#else
            DenseVector<DataType_> y_ref(b.copy());
            Difference<tags::CPU>::value(y_ref ,Product<tags::CPU>::value(bm1, dv1));
#endif

            TEST(prod2.lock(lm_read_only),
                    for (typename DenseVector<DataType_>::ConstElementIterator dit(y_ref.begin_elements()), it(prod2.begin_elements()), i_end(prod2.end_elements()) ;
                        it != i_end ; ++it, ++dit)
                    {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*it, *dit, std::numeric_limits<DataType_>::epsilon());
                    },
                    prod2.unlock(lm_read_only));
            }

            DenseVector<DataType_> dv01(4, DataType_(1));
            DenseVector<DataType_> dv02(1089, DataType_(1));
            BandedMatrixQx<Q1Type, DataType_> bm01(1089, dv02, dv02, dv02, dv02, dv02, dv02, dv02, dv02, dv02);
            TEST_CHECK_THROWS(Product<Tag_>::value(bm01, dv01), VectorSizeDoesNotMatch);
        }
};
Q1MatrixDenseVectorDefectTest<tags::CPU, float> q1_defect_test_float("float");
Q1MatrixDenseVectorDefectTest<tags::CPU, double> q1_defect_test_double("double");
Q1MatrixDenseVectorDefectTest<tags::CPU::MultiCore, float> q1_prod_mc_test_float("MC float");
Q1MatrixDenseVectorDefectTest<tags::CPU::MultiCore, double> q1_prod_mc_test_double("MC double");
#ifdef HONEI_SSE
Q1MatrixDenseVectorDefectTest<tags::CPU::SSE, float> sse_q1_defect_test_float("float");
Q1MatrixDenseVectorDefectTest<tags::CPU::SSE, double> sse_q1_defect_test_double("double");
Q1MatrixDenseVectorDefectTest<tags::CPU::MultiCore::SSE, float> q1_prod_mc_sse_test_float("MC SSE float");
Q1MatrixDenseVectorDefectTest<tags::CPU::MultiCore::SSE, double> q1_prod_mc_sse_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
Q1MatrixDenseVectorDefectTest<tags::GPU::CUDA, float> cuda_q1_defect_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
Q1MatrixDenseVectorDefectTest<tags::GPU::CUDA, double> cuda_q1_defect_test_double("double");
#endif
#endif
#ifdef HONEI_OPENCL
Q1MatrixDenseVectorDefectTest<tags::OpenCL::CPU, float> ocl_cpu_q1_defect_test_float("float");
Q1MatrixDenseVectorDefectTest<tags::OpenCL::CPU, double> ocl_cpu_q1_defect_test_double("double");
Q1MatrixDenseVectorDefectTest<tags::OpenCL::GPU, float> ocl_gpu_q1_defect_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
Q1MatrixDenseVectorDefectTest<tags::OpenCL::GPU, double> ocl_gpu_q1_defect_test_double("double");
#endif
#endif
