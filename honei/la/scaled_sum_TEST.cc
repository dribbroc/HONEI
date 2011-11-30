/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
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

#include <honei/la/dense_vector.hh>
#include <honei/la/norm.hh>
#include <honei/la/scaled_sum.hh>
#include <honei/la/sparse_vector.hh>
#include <honei/util/unittest.hh>

#include <limits>

using namespace honei;
using namespace tests;

template <typename Tag_, typename DataType_>
class DenseVectorScaledSumTest :
    public BaseTest
{
    public:
        DenseVectorScaledSumTest(const std::string & type) :
            BaseTest("dense_vector_scaled_sum_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 10) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size, DataType_(2));
                DenseVector<DataType_> dv2(size, DataType_(3));
                DataType_ scal(DataType_(2));

                ScaledSum<Tag_>::value(dv1, dv2, scal);
                DenseVector<DataType_> dv3(size, DataType_(8));

                TEST(dv1.lock(lm_read_only), TEST_CHECK_EQUAL(dv1, dv3), dv1.unlock(lm_read_only));
            }

            DenseVector<DataType_> dv00(1, DataType_(1));
            DenseVector<DataType_> dv01(2, DataType_(1));
            DataType_ scal00(DataType_(2));

            TEST_CHECK_THROWS(ScaledSum<Tag_>::value(dv00, dv01, scal00), VectorSizeDoesNotMatch);
        }
};
DenseVectorScaledSumTest<tags::CPU, float> dense_vector_scaled_sum_test_float("float");
DenseVectorScaledSumTest<tags::CPU, double> dense_vector_scaled_sum_test_double("double");
DenseVectorScaledSumTest<tags::CPU::MultiCore, float> mc_dense_vector_scaled_sum_test_float("MC float");
DenseVectorScaledSumTest<tags::CPU::MultiCore, double> mc_dense_vector_scaled_sum_test_double("MC double");
DenseVectorScaledSumTest<tags::CPU::Generic, float> generic_dense_vector_scaled_sum_test_float("float");
DenseVectorScaledSumTest<tags::CPU::Generic, double> generic_dense_vector_scaled_sum_test_double("double");
DenseVectorScaledSumTest<tags::CPU::MultiCore::Generic, float> generic_mc_dense_vector_scaled_sum_test_float("MC float");
DenseVectorScaledSumTest<tags::CPU::MultiCore::Generic, double> generic_mc_dense_vector_scaled_sum_test_double("MC double");
#ifdef HONEI_SSE
DenseVectorScaledSumTest<tags::CPU::SSE, float> sse_dense_vector_scaled_sum_test_float("SSE float");
DenseVectorScaledSumTest<tags::CPU::SSE, double> sse_dense_vector_scaled_sum_test_double("SSE double");
DenseVectorScaledSumTest<tags::CPU::MultiCore::SSE, float>
    mc_sse_dense_vector_scaled_sum_test_float("MC SSE float");
DenseVectorScaledSumTest<tags::CPU::MultiCore::SSE, double>
    mc_sse_dense_vector_scaled_sum_test_double("MC SSE double");
#endif
#ifdef HONEI_ITANIUM
DenseVectorScaledSumTest<tags::CPU::Itanium, float> itanium_dense_vector_scaled_sum_test_float("float");
DenseVectorScaledSumTest<tags::CPU::Itanium, double> itanium_dense_vector_scaled_sum_test_double("double");
#endif
#ifdef HONEI_CUDA
DenseVectorScaledSumTest<tags::GPU::CUDA, float> cuda_dense_vector_scaled_sum_test_float("float");
DenseVectorScaledSumTest<tags::GPU::MultiCore::CUDA, float> mc_cuda_dense_vector_scaled_sum_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
DenseVectorScaledSumTest<tags::GPU::CUDA, double> cuda_dense_vector_scaled_sum_test_double("double");
DenseVectorScaledSumTest<tags::GPU::MultiCore::CUDA, double> mc_cuda_dense_vector_scaled_sum_test_double("double");
#endif
#endif
#ifdef HONEI_OPENCL
DenseVectorScaledSumTest<tags::OpenCL::CPU, float> opencl_cpu_dense_vector_scaled_sum_test_float("float");
DenseVectorScaledSumTest<tags::OpenCL::CPU, double> opencl_cpu_dense_vector_scaled_sum_test_double("double");
DenseVectorScaledSumTest<tags::OpenCL::GPU, float> opencl_gpu_dense_vector_scaled_sum_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
DenseVectorScaledSumTest<tags::OpenCL::GPU, double> opencl_gpu_dense_vector_scaled_sum_test_double("double");
#endif
#endif
#ifdef HONEI_CELL
DenseVectorScaledSumTest<tags::Cell, float> cell_dense_vector_scaled_sum_test_float("Cell float");
DenseVectorScaledSumTest<tags::Cell, double> cell_dense_vector_scaled_sum_test_double("Cell double");
#endif

template <typename Tag_, typename DataType_>
class DenseVectorScaledSumQuickTest :
    public QuickTest
{
    public:
        DenseVectorScaledSumQuickTest(const std::string & type) :
            QuickTest("dense_vector_scaled_sum_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(65);

            DenseVector<DataType_> dv1(size, DataType_(2));
            DenseVector<DataType_> dv2(size, DataType_(3));
            DataType_ scal(DataType_(2));

            ScaledSum<Tag_>::value(dv1, dv2, scal);
            DenseVector<DataType_> dv3(size, DataType_(8));

            TEST(dv1.lock(lm_read_only), TEST_CHECK_EQUAL(dv1, dv3), dv1.unlock(lm_read_only));

            DenseVector<DataType_> dv00(1, DataType_(1));
            DenseVector<DataType_> dv01(2, DataType_(1));
            DataType_ scal00(DataType_(2));

            TEST_CHECK_THROWS(ScaledSum<Tag_>::value(dv00, dv01, scal00), VectorSizeDoesNotMatch);
        }
};
DenseVectorScaledSumQuickTest<tags::CPU, float> dense_vector_scaled_sum_quick_test_float("float");
DenseVectorScaledSumQuickTest<tags::CPU, double> dense_vector_scaled_sum_quick_test_double("double");
DenseVectorScaledSumQuickTest<tags::CPU::MultiCore, float> mc_dense_vector_scaled_sum_quick_test_float("MC float");
DenseVectorScaledSumQuickTest<tags::CPU::MultiCore, double> mc_dense_vector_scaled_sum_quick_test_double("MC double");
DenseVectorScaledSumQuickTest<tags::CPU::Generic, float> generic_dense_vector_scaled_sum_quick_test_float("float");
DenseVectorScaledSumQuickTest<tags::CPU::Generic, double> generic_dense_vector_scaled_sum_quick_test_double("double");
DenseVectorScaledSumQuickTest<tags::CPU::MultiCore::Generic, float> generic_mc_dense_vector_scaled_sum_quick_test_float("MC float");
DenseVectorScaledSumQuickTest<tags::CPU::MultiCore::Generic, double> generic_mc_dense_vector_scaled_sum_quick_test_double("MC double");
#ifdef HONEI_SSE
DenseVectorScaledSumQuickTest<tags::CPU::SSE, float> sse_dense_vector_scaled_sum_quick_test_float("SSE float");
DenseVectorScaledSumQuickTest<tags::CPU::SSE, double> sse_dense_vector_scaled_sum_quick_test_double("SSE double");
DenseVectorScaledSumQuickTest<tags::CPU::MultiCore::SSE, float> mc_sse_dense_vector_scaled_sum_quick_test_float("MC SSE float");
DenseVectorScaledSumQuickTest<tags::CPU::MultiCore::SSE, double> mc_sse_dense_vector_scaled_sum_quick_test_double("MC SSE double");
#endif
#ifdef HONEI_ITANIUM
DenseVectorScaledSumQuickTest<tags::CPU::Itanium, float> itanium_dense_vector_scaled_sum_quick_test_float("float");
DenseVectorScaledSumQuickTest<tags::CPU::Itanium, double> itanium_dense_vector_scaled_sum_quick_test_double("double");
#endif
#ifdef HONEI_CUDA
DenseVectorScaledSumQuickTest<tags::GPU::CUDA, float> cuda_dense_vector_scaled_sum_quick_test_float("float");
DenseVectorScaledSumQuickTest<tags::GPU::MultiCore::CUDA, float> mc_cuda_dense_vector_scaled_sum_quick_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
DenseVectorScaledSumQuickTest<tags::GPU::CUDA, double> cuda_dense_vector_scaled_sum_quick_test_double("double");
DenseVectorScaledSumQuickTest<tags::GPU::MultiCore::CUDA, double> mc_cuda_dense_vector_scaled_sum_quick_test_double("double");
#endif
#endif
#ifdef HONEI_OPENCL
DenseVectorScaledSumQuickTest<tags::OpenCL::CPU, float> opencl_cpu_dense_vector_scaled_sum_quick_test_float("float");
DenseVectorScaledSumQuickTest<tags::OpenCL::CPU, double> opencl_cpu_dense_vector_scaled_sum_quick_test_double("double");
DenseVectorScaledSumQuickTest<tags::OpenCL::GPU, float> opencl_gpu_dense_vector_scaled_sum_quick_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
DenseVectorScaledSumQuickTest<tags::OpenCL::GPU, double> opencl_gpu_dense_vector_scaled_sum_quick_test_double("double");
#endif
#endif
#ifdef HONEI_CELL
DenseVectorScaledSumQuickTest<tags::Cell, float> cell_dense_vector_scaled_sum_quick_test_float("Cell float");
DenseVectorScaledSumQuickTest<tags::Cell, double> cell_dense_vector_scaled_sum_quick_test_double("Cell double");
#endif

template <typename Tag_, typename DataType_>
class DenseVectorResScaledSumTest :
    public BaseTest
{
    public:
        DenseVectorResScaledSumTest(const std::string & type) :
            BaseTest("dense_vector_res_scaled_sum_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 10) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size, DataType_(2));
                DenseVector<DataType_> dv2(size, DataType_(3));
                DataType_ scal(DataType_(2));
                DenseVector<DataType_> dv4(size);
                ScaledSum<Tag_>::value(dv4, dv1, dv2, scal);

                ScaledSum<Tag_>::value(dv1, dv2, scal);
                DenseVector<DataType_> dv3(size, DataType_(8));

                TEST(dv1.lock(lm_read_only), TEST_CHECK_EQUAL(dv1, dv3), dv1.unlock(lm_read_only));
                TEST(dv4.lock(lm_read_only), TEST_CHECK_EQUAL(dv4, dv3), dv4.unlock(lm_read_only));
            }

            DenseVector<DataType_> dv00(1, DataType_(1));
            DenseVector<DataType_> dv01(2, DataType_(1));
            DataType_ scal00(DataType_(2));

            TEST_CHECK_THROWS(ScaledSum<Tag_>::value(dv00, dv01, scal00), VectorSizeDoesNotMatch);
        }
};
DenseVectorResScaledSumTest<tags::CPU, float> dense_vector_scaled_res_sum_test_float("float");
DenseVectorResScaledSumTest<tags::CPU, double> dense_vector_scaled_res_sum_test_double("double");
DenseVectorResScaledSumTest<tags::CPU::Generic, float> generic_dense_vector_scaled_res_sum_test_float("float");
DenseVectorResScaledSumTest<tags::CPU::Generic, double> generic_dense_vector_scaled_res_sum_test_double("double");
DenseVectorResScaledSumTest<tags::CPU::MultiCore::Generic, float> generic_mc_dense_vector_scaled_res_sum_test_float("float");
DenseVectorResScaledSumTest<tags::CPU::MultiCore::Generic, double> generic_mc_dense_vector_scaled_res_sum_test_double("double");
#ifdef HONEI_SSE
DenseVectorResScaledSumTest<tags::CPU::SSE, float> sse_dense_vector_scaled_res_sum_test_float("float");
DenseVectorResScaledSumTest<tags::CPU::SSE, double> sse_dense_vector_scaled_res_sum_test_double("double");
DenseVectorResScaledSumTest<tags::CPU::MultiCore::SSE, float> mc_sse_dense_vector_scaled_res_sum_test_float("float");
DenseVectorResScaledSumTest<tags::CPU::MultiCore::SSE, double> mc_sse_dense_vector_scaled_res_sum_test_double("double");
#endif
#ifdef HONEI_CUDA
DenseVectorResScaledSumTest<tags::GPU::CUDA, float> cuda_dense_vector_res_scaled_sum_test_float("float");
DenseVectorResScaledSumTest<tags::GPU::MultiCore::CUDA, float> mc_cuda_dense_vector_res_scaled_sum_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
DenseVectorResScaledSumTest<tags::GPU::CUDA, double> cuda_dense_vector_res_scaled_sum_test_double("double");
DenseVectorResScaledSumTest<tags::GPU::MultiCore::CUDA, double> mc_cuda_dense_vector_res_scaled_sum_test_double("double");
#endif
#endif
#ifdef HONEI_OPENCL
DenseVectorResScaledSumTest<tags::OpenCL::CPU, float> ocl_cpu_dense_vector_scaled_res_sum_test_float("float");
DenseVectorResScaledSumTest<tags::OpenCL::CPU, double> ocl_cpu_dense_vector_scaled_res_sum_test_double("double");
DenseVectorResScaledSumTest<tags::OpenCL::GPU, float> ocl_gpu_dense_vector_scaled_res_sum_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
DenseVectorResScaledSumTest<tags::OpenCL::GPU, double> ocl_gpu_dense_vector_scaled_res_sum_test_double("double");
#endif
#endif

template <typename Tag_, typename DataType_>
class DenseVectorResScaledSumQuickTest :
    public QuickTest
{
    public:
        DenseVectorResScaledSumQuickTest(const std::string & type) :
            QuickTest("dense_vector_Res_scaled_sum_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(65);

            DenseVector<DataType_> dv1(size, DataType_(2));
            DenseVector<DataType_> dv2(size, DataType_(3));
            DataType_ scal(DataType_(2));
            DenseVector<DataType_> dv4(size);
            ScaledSum<Tag_>::value(dv4, dv1, dv2, scal);

            ScaledSum<Tag_>::value(dv1, dv2, scal);
            DenseVector<DataType_> dv3(size, DataType_(8));

            TEST(dv1.lock(lm_read_only), TEST_CHECK_EQUAL(dv1, dv3), dv1.unlock(lm_read_only));
            TEST(dv4.lock(lm_read_only), TEST_CHECK_EQUAL(dv4, dv3), dv4.unlock(lm_read_only));

            DenseVector<DataType_> dv00(1, DataType_(1));
            DenseVector<DataType_> dv01(2, DataType_(1));
            DataType_ scal00(DataType_(2));

            TEST_CHECK_THROWS(ScaledSum<Tag_>::value(dv00, dv01, scal00), VectorSizeDoesNotMatch);
        }
};
DenseVectorResScaledSumQuickTest<tags::CPU, float> dense_vector_res_scaled_sum_quick_test_float("float");
DenseVectorResScaledSumQuickTest<tags::CPU, double> dense_vector_res_scaled_sum_quick_test_double("double");
DenseVectorResScaledSumQuickTest<tags::CPU::Generic, float> generic_dense_vector_res_scaled_sum_quick_test_float("float");
DenseVectorResScaledSumQuickTest<tags::CPU::Generic, double> generic_dense_vector_res_scaled_sum_quick_test_double("double");
DenseVectorResScaledSumQuickTest<tags::CPU::MultiCore::Generic, float> mc_generic_dense_vector_res_scaled_sum_quick_test_float("float");
DenseVectorResScaledSumQuickTest<tags::CPU::MultiCore::Generic, double> mc_generic_dense_vector_res_scaled_sum_quick_test_double("double");
#ifdef HONEI_SSE
DenseVectorResScaledSumQuickTest<tags::CPU::SSE, float> sse_dense_vector_res_scaled_sum_quick_test_float("float");
DenseVectorResScaledSumQuickTest<tags::CPU::SSE, double> sse_dense_vector_res_scaled_sum_quick_test_double("double");
DenseVectorResScaledSumQuickTest<tags::CPU::MultiCore::SSE, float> mc_sse_dense_vector_res_scaled_sum_quick_test_float("float");
DenseVectorResScaledSumQuickTest<tags::CPU::MultiCore::SSE, double> mc_sse_dense_vector_res_scaled_sum_quick_test_double("double");
#endif
#ifdef HONEI_CUDA
DenseVectorResScaledSumQuickTest<tags::GPU::CUDA, float> cuda_dense_vector_res_scaled_sum_quick_test_float("float");
DenseVectorResScaledSumQuickTest<tags::GPU::MultiCore::CUDA, float> mc_cuda_dense_vector_res_scaled_sum_quick_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
DenseVectorResScaledSumQuickTest<tags::GPU::CUDA, double> cuda_dense_vector_res_scaled_sum_quick_test_double("double");
DenseVectorResScaledSumQuickTest<tags::GPU::MultiCore::CUDA, double> mc_cuda_dense_vector_res_scaled_sum_quick_test_double("double");
#endif
#endif
#ifdef HONEI_OPENCL
DenseVectorResScaledSumQuickTest<tags::OpenCL::CPU, float> ocl_cpu_dense_vector_scaled_res_sum_quick_test_float("float");
DenseVectorResScaledSumQuickTest<tags::OpenCL::CPU, double> ocl_cpu_dense_vector_scaled_res_sum_quick_test_double("double");
DenseVectorResScaledSumQuickTest<tags::OpenCL::GPU, float> ocl_gpu_dense_vector_scaled_res_sum_quick_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
DenseVectorResScaledSumQuickTest<tags::OpenCL::GPU, double> ocl_gpu_dense_vector_scaled_res_sum_quick_test_double("double");
#endif
#endif

template <typename Tag_, typename DataType_>
class DenseVectorRangeScaledSumTest :
    public BaseTest
{
    public:
        DenseVectorRangeScaledSumTest(const std::string & type) :
            BaseTest("dense_vector_range_scaled_sum_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 10) ; size <<= 1)
            {
                for (int i(0) ; i < 4 ; ++i)
                {
                    for (int j(0) ; j < 4 ; ++j)
                    {
                        DenseVector<DataType_> dv1_source(size * 2, DataType_(2));
                        DenseVectorRange<DataType_> dv1(dv1_source, size, i);

                        DenseVector<DataType_> dv2_source(size * 2, DataType_(3));
                        DenseVectorRange<DataType_> dv2(dv2_source, size, j);
                        DataType_ scal(DataType_(2));

                        ScaledSum<Tag_>::value(dv1, dv2, scal);
                        DenseVector<DataType_> dv3(size, DataType_(8));

                        TEST(dv1.lock(lm_read_only), TEST_CHECK_EQUAL(dv1, dv3), dv1.unlock(lm_read_only));
                    }
                }
            }

            DenseVector<DataType_> dv00(1, DataType_(1));
            DenseVector<DataType_> dv01(2, DataType_(1));
            DataType_ scal00(DataType_(2));

            TEST_CHECK_THROWS(ScaledSum<Tag_>::value(dv00, dv01, scal00), VectorSizeDoesNotMatch);
        }
};
DenseVectorRangeScaledSumTest<tags::CPU, float> dense_vector_range_scaled_sum_test_float("float");
DenseVectorRangeScaledSumTest<tags::CPU, double> dense_vector_range_scaled_sum_test_double("double");
DenseVectorRangeScaledSumTest<tags::CPU::MultiCore, float> mc_dense_vector_range_scaled_sum_test_float("MC float");
DenseVectorRangeScaledSumTest<tags::CPU::MultiCore, double> mc_dense_vector_range_scaled_sum_test_double("MC double");
DenseVectorRangeScaledSumTest<tags::CPU::Generic, float> generic_dense_vector_range_scaled_sum_test_float("float");
DenseVectorRangeScaledSumTest<tags::CPU::Generic, double> generic_dense_vector_range_scaled_sum_test_double("double");
DenseVectorRangeScaledSumTest<tags::CPU::MultiCore::Generic, float> generic_mc_dense_vector_range_scaled_sum_test_float("MC float");
DenseVectorRangeScaledSumTest<tags::CPU::MultiCore::Generic, double> generic_mc_dense_vector_range_scaled_sum_test_double("MC double");
#ifdef HONEI_SSE
DenseVectorRangeScaledSumTest<tags::CPU::SSE, float> sse_dense_vector_range_scaled_sum_test_float("SSE float");
DenseVectorRangeScaledSumTest<tags::CPU::SSE, double> sse_dense_vector_range_scaled_sum_test_double("SSE double");
DenseVectorRangeScaledSumTest<tags::CPU::MultiCore::SSE, float> mc_sse_dense_vector_range_scaled_sum_test_float("MC SSE float");
DenseVectorRangeScaledSumTest<tags::CPU::MultiCore::SSE, double> mc_sse_dense_vector_range_scaled_sum_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
DenseVectorRangeScaledSumTest<tags::GPU::CUDA, float> cuda_dense_vector_range_scaled_sum_test_float("float");
DenseVectorRangeScaledSumTest<tags::GPU::MultiCore::CUDA, float> mc_cuda_dense_vector_range_scaled_sum_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
DenseVectorRangeScaledSumTest<tags::GPU::CUDA, double> cuda_dense_vector_range_scaled_sum_test_double("double");
DenseVectorRangeScaledSumTest<tags::GPU::MultiCore::CUDA, double> mc_cuda_dense_vector_range_scaled_sum_test_double("double");
#endif
#endif
#ifdef HONEI_OPENCL
DenseVectorRangeScaledSumTest<tags::OpenCL::CPU, float> opencl_cpu_dense_vector_range_scaled_sum_test_float("float");
DenseVectorRangeScaledSumTest<tags::OpenCL::CPU, double> opencl_cpu_dense_vector_range_scaled_sum_test_double("double");
DenseVectorRangeScaledSumTest<tags::OpenCL::GPU, float> opencl_gpu_dense_vector_range_scaled_sum_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
DenseVectorRangeScaledSumTest<tags::OpenCL::GPU, double> opencl_gpu_dense_vector_range_scaled_sum_test_double("double");
#endif
#endif
#ifdef HONEI_CELL
DenseVectorRangeScaledSumTest<tags::Cell, float> cell_dense_vector_range_scaled_sum_test_float("Cell float");
DenseVectorRangeScaledSumTest<tags::Cell, double> cell_dense_vector_range_scaled_sum_test_double("Cell double");
#endif

template <typename Tag_, typename DataType_>
class DenseVectorRangeScaledSumQuickTest :
    public QuickTest
{
    public:
        DenseVectorRangeScaledSumQuickTest(const std::string & type) :
            QuickTest("dense_vector_range_scaled_sum_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(65);

            DenseVector<DataType_> dv1_source(size * 2, DataType_(2));
            DenseVectorRange<DataType_> dv1(dv1_source, size, 1);

            DenseVector<DataType_> dv2_source(size * 2, DataType_(3));
            DenseVectorRange<DataType_> dv2(dv2_source, size, 3);

            DataType_ scal(DataType_(2));
            DenseVector<DataType_> dv3(size, DataType_(8));
            ScaledSum<Tag_>::value(dv1, dv2, scal);
            TEST(dv1.lock(lm_read_only), TEST_CHECK_EQUAL(dv1, dv3), dv1.unlock(lm_read_only));

            DenseVector<DataType_> dv00(1, DataType_(1));
            DenseVector<DataType_> dv01(2, DataType_(1));
            DataType_ scal00(DataType_(2));

            TEST_CHECK_THROWS(ScaledSum<Tag_>::value(dv00, dv01, scal00), VectorSizeDoesNotMatch);
        }
};
DenseVectorRangeScaledSumQuickTest<tags::CPU, float> dense_vector_range_scaled_sum_quick_test_float("float");
DenseVectorRangeScaledSumQuickTest<tags::CPU, double> dense_vector_range_scaled_sum_quick_test_double("double");
DenseVectorRangeScaledSumQuickTest<tags::CPU::MultiCore, float> mc_dense_vector_range_scaled_sum_quick_test_float("MC float");
DenseVectorRangeScaledSumQuickTest<tags::CPU::MultiCore, double> mc_dense_vector_range_scaled_sum_quick_test_double("MC double");
DenseVectorRangeScaledSumQuickTest<tags::CPU::Generic, float> generic_dense_vector_range_scaled_sum_quick_test_float("float");
DenseVectorRangeScaledSumQuickTest<tags::CPU::Generic, double> generic_dense_vector_range_scaled_sum_quick_test_double("double");
DenseVectorRangeScaledSumQuickTest<tags::CPU::MultiCore::Generic, float> generic_mc_dense_vector_range_scaled_sum_quick_test_float("MC float");
DenseVectorRangeScaledSumQuickTest<tags::CPU::MultiCore::Generic, double> generic_mc_dense_vector_range_scaled_sum_quick_test_double("MC double");
#ifdef HONEI_SSE
DenseVectorRangeScaledSumQuickTest<tags::CPU::SSE, float> sse_dense_vector_range_scaled_sum_quick_test_float("SSE float");
DenseVectorRangeScaledSumQuickTest<tags::CPU::SSE, double> sse_dense_vector_range_scaled_sum_quick_test_double("SSE double");
DenseVectorRangeScaledSumQuickTest<tags::CPU::MultiCore::SSE, float> mc_sse_dense_vector_range_scaled_sum_quick_test_float("MC SSE float");
DenseVectorRangeScaledSumQuickTest<tags::CPU::MultiCore::SSE, double> mc_sse_dense_vector_range_scaled_sum_quick_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
DenseVectorRangeScaledSumQuickTest<tags::GPU::CUDA, float> cuda_dense_vector_range_scaled_sum_quick_test_float("float");
DenseVectorRangeScaledSumQuickTest<tags::GPU::MultiCore::CUDA, float> mc_cuda_dense_vector_range_scaled_sum_quick_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
DenseVectorRangeScaledSumQuickTest<tags::GPU::CUDA, double> cuda_dense_vector_range_scaled_sum_quick_test_double("double");
DenseVectorRangeScaledSumQuickTest<tags::GPU::MultiCore::CUDA, double> mc_cuda_dense_vector_range_scaled_sum_quick_test_double("double");
#endif
#endif
#ifdef HONEI_OPENCL
DenseVectorRangeScaledSumQuickTest<tags::OpenCL::CPU, float> opencl_cpu_dense_vector_range_scaled_sum_quick_test_float("float");
DenseVectorRangeScaledSumQuickTest<tags::OpenCL::CPU, double> opencl_cpu_dense_vector_range_scaled_sum_quick_test_double("double");
DenseVectorRangeScaledSumQuickTest<tags::OpenCL::GPU, float> opencl_gpu_dense_vector_range_scaled_sum_quick_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
DenseVectorRangeScaledSumQuickTest<tags::OpenCL::GPU, double> opencl_gpu_dense_vector_range_scaled_sum_quick_test_double("double");
#endif
#endif
#ifdef HONEI_CELL
DenseVectorRangeScaledSumQuickTest<tags::Cell, float> cell_dense_vector_range_scaled_sum_quick_test_float("Cell float");
DenseVectorRangeScaledSumQuickTest<tags::Cell, double> cell_dense_vector_range_scaled_sum_quick_test_double("Cell double");
#endif

template <typename Tag_, typename DataType_>
class DenseVectorSliceScaledSumTest :
    public BaseTest
{
    public:
        DenseVectorSliceScaledSumTest(const std::string & type) :
            BaseTest("dense_vector_Slice_scaled_sum_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 10) ; size <<= 1)
            {
                for (int i(1) ; i < 4 ; ++i)
                {
                    for (int j(1) ; j < 4 ; ++j)
                    {
                        DenseVector<DataType_> dv1_source(size * 6, DataType_(2));
                        DenseVectorSlice<DataType_> dv1(dv1_source, size, i, j);

                        DenseVector<DataType_> dv2_source(size * 6, DataType_(3));
                        DenseVectorSlice<DataType_> dv2(dv2_source, size, i, j);
                        DataType_ scal(DataType_(2));

                        ScaledSum<Tag_>::value(dv1, dv2, scal);
                        DenseVector<DataType_> dv3(size, DataType_(8));

                        TEST(dv1.lock(lm_read_only), TEST_CHECK_EQUAL(dv1, dv3), dv1.unlock(lm_read_only));
                    }
                }
            }

            DenseVector<DataType_> dv00(1, DataType_(1));
            DenseVector<DataType_> dv01(2, DataType_(1));
            DataType_ scal00(DataType_(2));

            TEST_CHECK_THROWS(ScaledSum<Tag_>::value(dv00, dv01, scal00), VectorSizeDoesNotMatch);
        }
};
DenseVectorSliceScaledSumTest<tags::CPU, float> dense_vector_Slice_scaled_sum_test_float("float");
DenseVectorSliceScaledSumTest<tags::CPU, double> dense_vector_Slice_scaled_sum_test_double("double");
DenseVectorSliceScaledSumTest<tags::CPU::MultiCore, float> mc_dense_vector_Slice_scaled_sum_test_float("MC float");
DenseVectorSliceScaledSumTest<tags::CPU::MultiCore, double> mc_dense_vector_Slice_scaled_sum_test_double("MC double");
/*
#ifdef HONEI_SSE
DenseVectorSliceScaledSumTest<tags::CPU::SSE, float> sse_dense_vector_Slice_scaled_sum_test_float("SSE float");
DenseVectorSliceScaledSumTest<tags::CPU::SSE, double> sse_dense_vector_Slice_scaled_sum_test_double("SSE double");
DenseVectorSliceScaledSumTest<tags::CPU::MultiCore::SSE, float> mc_sse_dense_vector_Slice_scaled_sum_test_float("MC SSE float");
DenseVectorSliceScaledSumTest<tags::CPU::MultiCore::SSE, double> mc_sse_dense_vector_Slice_scaled_sum_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
DenseVectorSliceScaledSumTest<tags::GPU::CUDA, float> cuda_dense_vector_Slice_scaled_sum_test_float("float");
#endif
#ifdef HONEI_CELL
DenseVectorSliceScaledSumTest<tags::Cell, float> cell_dense_vector_Slice_scaled_sum_test_float("Cell float");
DenseVectorSliceScaledSumTest<tags::Cell, double> cell_dense_vector_Slice_scaled_sum_test_double("Cell double");
#endif
*/

template <typename Tag_, typename DataType_>
class DenseVectorSliceScaledSumQuickTest :
    public QuickTest
{
    public:
        DenseVectorSliceScaledSumQuickTest(const std::string & type) :
            QuickTest("dense_vector_Slice_scaled_sum_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(165);

            DenseVector<DataType_> dv1_source(size * 3, DataType_(2));
            DenseVectorSlice<DataType_> dv1(dv1_source, size , 5, 2);

            DenseVector<DataType_> dv2_source(size * 2, DataType_(3));
            DenseVectorSlice<DataType_> dv2(dv2_source, size, 2, 1);

            DataType_ scal(DataType_(2));
            DenseVector<DataType_> dv3(size, DataType_(8));
            ScaledSum<Tag_>::value(dv1, dv2, scal);
            TEST(dv1.lock(lm_read_only), TEST_CHECK_EQUAL(dv1, dv3), dv1.unlock(lm_read_only));

            DenseVector<DataType_> dv00(1, DataType_(1));
            DenseVector<DataType_> dv01(2, DataType_(1));
            DataType_ scal00(DataType_(2));

            TEST_CHECK_THROWS(ScaledSum<Tag_>::value(dv00, dv01, scal00), VectorSizeDoesNotMatch);
        }
};
DenseVectorSliceScaledSumQuickTest<tags::CPU, float> dense_vector_Slice_scaled_sum_quick_test_float("float");
DenseVectorSliceScaledSumQuickTest<tags::CPU, double> dense_vector_Slice_scaled_sum_quick_test_double("double");
DenseVectorSliceScaledSumQuickTest<tags::CPU::MultiCore, float> mc_dense_vector_Slice_scaled_sum_quick_test_float("MC float");
DenseVectorSliceScaledSumQuickTest<tags::CPU::MultiCore, double> mc_dense_vector_Slice_scaled_sum_quick_test_double("MC double");
/*
#ifdef HONEI_SSE
DenseVectorSliceScaledSumQuickTest<tags::CPU::SSE, float> sse_dense_vector_Slice_scaled_sum_quick_test_float("SSE float");
DenseVectorSliceScaledSumQuickTest<tags::CPU::SSE, double> sse_dense_vector_Slice_scaled_sum_quick_test_double("SSE double");
DenseVectorSliceScaledSumQuickTest<tags::CPU::MultiCore::SSE, float> mc_sse_dense_vector_Slice_scaled_sum_quick_test_float("MC SSE float");
DenseVectorSliceScaledSumQuickTest<tags::CPU::MultiCore::SSE, double> mc_sse_dense_vector_Slice_scaled_sum_quick_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
DenseVectorSliceScaledSumQuickTest<tags::GPU::CUDA, float> cuda_dense_vector_Slice_scaled_sum_quick_test_float("float");
#endif
#ifdef HONEI_CELL
DenseVectorSliceScaledSumQuickTest<tags::Cell, float> cell_dense_vector_Slice_scaled_sum_quick_test_float("Cell float");
DenseVectorSliceScaledSumQuickTest<tags::Cell, double> cell_dense_vector_Slice_scaled_sum_quick_test_double("Cell double");
#endif
*/

template <typename Tag_, typename DataType_>
class DenseVector3ScaledSumTest :
    public BaseTest
{
    public:
        DenseVector3ScaledSumTest(const std::string & type) :
            BaseTest("dense_vector_3_scaled_sum_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(11) ; size < (1 << 10) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size, DataType_(2));
                DenseVector<DataType_> dv2(size, DataType_(3));
                DenseVector<DataType_> sum(size, DataType_(2));
                DenseVectorRange<DataType_> dvr1(dv1, size / 2, 3);
                DenseVectorRange<DataType_> dvr2(dv2, size / 2, 2);
                DenseVectorRange<DataType_> sumr(sum, size / 2, 5);

                ScaledSum<Tag_>::value(sumr, dvr1, dvr2);
                sumr.lock(lm_read_only);
                DataType_ v1(Norm<vnt_l_one>::value(sumr));
                sumr.unlock(lm_read_only);

                TEST_CHECK_EQUAL(v1, 8 * (size / 2));
            }

            DenseVector<DataType_> dv00(1, DataType_(1));
            DenseVector<DataType_> dv01(2, DataType_(1));
            DenseVector<DataType_> sum00(1, DataType_(2));

            TEST_CHECK_THROWS(ScaledSum<Tag_>::value(sum00, dv00, dv01), VectorSizeDoesNotMatch);
        }
};
DenseVector3ScaledSumTest<tags::CPU, float> dense_vector_3_scaled_sum_test_float("float");
DenseVector3ScaledSumTest<tags::CPU, double> dense_vector_scaled_3_sum_test_double("double");
DenseVector3ScaledSumTest<tags::CPU::MultiCore, float> mc_dense_vector_3_scaled_sum_test_float("MC float");
DenseVector3ScaledSumTest<tags::CPU::MultiCore, double> mc_dense_vector_scaled_3_sum_test_double("MC double");
#ifdef HONEI_SSE
DenseVector3ScaledSumTest<tags::CPU::SSE, float> sse_dense_vector_3_scaled_sum_test_float("SSE float");
DenseVector3ScaledSumTest<tags::CPU::SSE, double> sse_dense_vector_3_scaled_sum_test_double("SSE double");
DenseVector3ScaledSumTest<tags::CPU::MultiCore::SSE, float> mc_sse_dense_vector_3_scaled_sum_test_float("MC SSE float");
DenseVector3ScaledSumTest<tags::CPU::MultiCore::SSE, double> mc_sse_dense_vector_scaled_3_sum_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
DenseVector3ScaledSumTest<tags::GPU::CUDA, float> cuda_dense_vector_3_scaled_sum_test_float("float");
DenseVector3ScaledSumTest<tags::GPU::MultiCore::CUDA, float> mc_cuda_dense_vector_3_scaled_sum_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
DenseVector3ScaledSumTest<tags::GPU::CUDA, double> cuda_dense_vector_3_scaled_sum_test_double("double");
DenseVector3ScaledSumTest<tags::GPU::MultiCore::CUDA, double> mc_cuda_dense_vector_3_scaled_sum_test_double("double");
#endif
#endif
#ifdef HONEI_OPENCL
DenseVector3ScaledSumTest<tags::OpenCL::CPU, float> ocl_cpu_dense_vector_3_scaled_sum_test_float("float");
DenseVector3ScaledSumTest<tags::OpenCL::CPU, double> ocl_cpu_dense_vector_scaled_3_sum_test_double("double");
DenseVector3ScaledSumTest<tags::OpenCL::GPU, float> ocl_gpu_dense_vector_3_scaled_sum_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
DenseVector3ScaledSumTest<tags::OpenCL::GPU, double> ocl_gpu_dense_vector_scaled_3_sum_test_double("double");
#endif
#endif

template <typename Tag_, typename DataType_>
class DenseVector3ScaledSumQuickTest :
    public QuickTest
{
    public:
        DenseVector3ScaledSumQuickTest(const std::string & type) :
            QuickTest("dense_vector_3_scaled_sum_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(65);
            DenseVector<DataType_> dv1(size, DataType_(2));
            DenseVector<DataType_> dv2(size, DataType_(3));
            DenseVector<DataType_> sum(size, DataType_(2));
            DenseVectorRange<DataType_> dvr1(dv1, size / 2, 3);
            DenseVectorRange<DataType_> dvr2(dv2, size / 2, 2);
            DenseVectorRange<DataType_> sumr(sum, size / 2, 5);

            ScaledSum<Tag_>::value(sumr, dvr1, dvr2);
            sumr.lock(lm_read_only);
            DataType_ v1(Norm<vnt_l_one>::value(sumr));
            sumr.unlock(lm_read_only);
            TEST_CHECK_EQUAL(v1, 8 * (size / 2));

            DenseVector<DataType_> dv00(1, DataType_(1));
            DenseVector<DataType_> dv01(2, DataType_(1));
            DenseVector<DataType_> sum00(1, DataType_(2));

            TEST_CHECK_THROWS(ScaledSum<Tag_>::value(sum00, dv00, dv01), VectorSizeDoesNotMatch);
        }
};
DenseVector3ScaledSumQuickTest<tags::CPU, float> dense_vector_3_scaled_sum_quick_test_float("float");
DenseVector3ScaledSumQuickTest<tags::CPU, double> dense_vector_scaled_3_sum_quick_test_double("double");
DenseVector3ScaledSumQuickTest<tags::CPU::MultiCore, float> mc_dense_vector_3_scaled_sum_quick_test_float("MC float");
DenseVector3ScaledSumQuickTest<tags::CPU::MultiCore, double> mc_dense_vector_scaled_3_sum_quick_test_double("MC double");
#ifdef HONEI_SSE
DenseVector3ScaledSumQuickTest<tags::CPU::SSE, float> sse_dense_vector_3_scaled_sum_quick_test_float("SSE float");
DenseVector3ScaledSumQuickTest<tags::CPU::SSE, double> sse_dense_vector_3_scaled_sum_quick_test_double("SSE double");
DenseVector3ScaledSumQuickTest<tags::CPU::MultiCore::SSE, float> mc_sse_dense_vector_3_scaled_sum_quick_test_float("MC SSE float");
DenseVector3ScaledSumQuickTest<tags::CPU::MultiCore::SSE, double> mc_sse_dense_vector_scaled_3_sum_quick_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
DenseVector3ScaledSumQuickTest<tags::GPU::CUDA, float> cuda_dense_vector_3_scaled_sum_quick_test_float("float");
DenseVector3ScaledSumQuickTest<tags::GPU::MultiCore::CUDA, float> mc_cuda_dense_vector_3_scaled_sum_quick_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
DenseVector3ScaledSumQuickTest<tags::GPU::CUDA, double> cuda_dense_vector_3_scaled_sum_quick_test_double("double");
DenseVector3ScaledSumQuickTest<tags::GPU::MultiCore::CUDA, double> mc_cuda_dense_vector_3_scaled_sum_quick_test_double("double");
#endif
#endif
#ifdef HONEI_OPENCL
DenseVector3ScaledSumQuickTest<tags::OpenCL::CPU, float> ocl_cpu_dense_vector_3_scaled_sum_quick_test_float("float");
DenseVector3ScaledSumQuickTest<tags::OpenCL::CPU, double> ocl_cpu_dense_vector_scaled_3_sum_quick_test_double("double");
DenseVector3ScaledSumQuickTest<tags::OpenCL::GPU, float> ocl_gpu_dense_vector_3_scaled_sum_quick_test_float("float");
#ifdef HONEI_CUDA_DOUBLE
DenseVector3ScaledSumQuickTest<tags::OpenCL::GPU, double> ocl_gpu_dense_vector_scaled_3_sum_quick_test_double("double");
#endif
#endif

template <typename Tag_, typename DataType_>
class DenseVectorSlice3ScaledSumTest :
    public BaseTest
{
    public:
        DenseVectorSlice3ScaledSumTest(const std::string & type) :
            BaseTest("dense_vector_slice_3_scaled_sum_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 10) ; size <<= 1)
            {
                for (int i(1) ; i < 4 ; ++i)
                {
                    for (int j(1) ; j < 4 ; ++j)
                    {
                        DenseVector<DataType_> dv1_source(size * 6, DataType_(2));
                        DenseVectorSlice<DataType_> dv1(dv1_source, size, i, j);

                        DenseVector<DataType_> dv2_source(size * 6, DataType_(3));
                        DenseVectorSlice<DataType_> dv2(dv2_source, size, i, j);

                        DenseVector<DataType_> dv3_source(size * 6, DataType_(4));
                        DenseVectorSlice<DataType_> dv3(dv3_source, size, i, j);

                        ScaledSum<Tag_>::value(dv1, dv2, dv3);
                        DenseVector<DataType_> dv4(size, DataType_(14));

                        TEST(dv1.lock(lm_read_only), TEST_CHECK_EQUAL(dv1, dv4), dv1.unlock(lm_read_only));
                    }
                }
            }

            DenseVector<DataType_> dv00(1, DataType_(1));
            DenseVector<DataType_> dv01(2, DataType_(1));
            DenseVector<DataType_> dv02(2, DataType_(1));

            TEST_CHECK_THROWS(ScaledSum<Tag_>::value(dv00, dv01, dv02), VectorSizeDoesNotMatch);
        }
};
DenseVectorSlice3ScaledSumTest<tags::CPU, float> dense_vector_Slice_3scaled_sum_test_float("float");
DenseVectorSlice3ScaledSumTest<tags::CPU, double> dense_vector_Slice_3scaled_sum_test_double("double");
DenseVectorSlice3ScaledSumTest<tags::CPU::MultiCore, float> mc_dense_vector_Slice_3scaled_sum_test_float("MC float");
DenseVectorSlice3ScaledSumTest<tags::CPU::MultiCore, double> mc_dense_vector_Slice_3scaled_sum_test_double("MC double");
/*
#ifdef HONEI_SSE
DenseVectorSlice3ScaledSumTest<tags::CPU::SSE, float> sse_dense_vector_Slice_3scaled_sum_test_float("SSE float");
DenseVectorSlice3ScaledSumTest<tags::CPU::SSE, double> sse_dense_vector_Slice_3scaled_sum_test_double("SSE double");
DenseVectorSlice3ScaledSumTest<tags::CPU::MultiCore::SSE, float> mc_sse_dense_vector_Slice_3scaled_sum_test_float("MC SSE float");
DenseVectorSlice3ScaledSumTest<tags::CPU::MultiCore::SSE, double> mc_sse_dense_vector_Slice_3scaled_sum_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
DenseVectorSlice3ScaledSumTest<tags::GPU::CUDA, float> cuda_dense_vector_Slice_3scaled_sum_test_float("float");
#endif
#ifdef HONEI_CELL
DenseVectorSlice3ScaledSumTest<tags::Cell, float> cell_dense_vector_Slice_3scaled_sum_test_float("Cell float");
DenseVectorSlice3ScaledSumTest<tags::Cell, double> cell_dense_vector_Slice_3scaled_sum_test_double("Cell double");
#endif
*/

template <typename Tag_, typename DataType_>
class DenseVectorSlice3ScaledSumQuickTest :
    public QuickTest
{
    public:
        DenseVectorSlice3ScaledSumQuickTest(const std::string & type) :
            QuickTest("dense_vector_slice_3_scaled_sum_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(165);

            DenseVector<DataType_> dv1_source(size * 3, DataType_(2));
            DenseVectorSlice<DataType_> dv1(dv1_source, size , 5, 2);

            DenseVector<DataType_> dv2_source(size * 2, DataType_(3));
            DenseVectorSlice<DataType_> dv2(dv2_source, size, 2, 1);

            DenseVector<DataType_> dv3_source(size * 2, DataType_(4));
            DenseVectorSlice<DataType_> dv3(dv3_source, size, 2, 1);

            DenseVector<DataType_> dv4(size, DataType_(14));
            ScaledSum<Tag_>::value(dv1, dv2, dv3);
            TEST(dv1.lock(lm_read_only), TEST_CHECK_EQUAL(dv1, dv4), dv1.unlock(lm_read_only));

            DenseVector<DataType_> dv00(1, DataType_(1));
            DenseVector<DataType_> dv01(2, DataType_(1));
            DenseVector<DataType_> dv02(2, DataType_(1));

            TEST_CHECK_THROWS(ScaledSum<Tag_>::value(dv00, dv01, dv02), VectorSizeDoesNotMatch);
        }
};
DenseVectorSlice3ScaledSumQuickTest<tags::CPU, float> dense_vector_Slice_3scaled_sum_quick_test_float("float");
DenseVectorSlice3ScaledSumQuickTest<tags::CPU, double> dense_vector_Slice_3scaled_sum_quick_test_double("double");
DenseVectorSlice3ScaledSumQuickTest<tags::CPU::MultiCore, float> mc_dense_vector_Slice_3scaled_sum_quick_test_float("MC float");
DenseVectorSlice3ScaledSumQuickTest<tags::CPU::MultiCore, double> mc_dense_vector_Slice_3scaled_sum_quick_test_double("MC double");
/*
#ifdef HONEI_SSE
DenseVectorSlice3ScaledSumQuickTest<tags::CPU::SSE, float> sse_dense_vector_Slice_3scaled_sum_quick_test_float("SSE float");
DenseVectorSlice3ScaledSumQuickTest<tags::CPU::SSE, double> sse_dense_vector_Slice_3scaled_sum_quick_test_double("SSE double");
DenseVectorSlice3ScaledSumQuickTest<tags::CPU::MultiCore::SSE, float> mc_sse_dense_vector_Slice_3scaled_sum_quick_test_float("MC SSE float");
DenseVectorSlice3ScaledSumQuickTest<tags::CPU::MultiCore::SSE, double> mc_sse_dense_vector_Slice_3scaled_sum_quick_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
DenseVectorSlice3ScaledSumQuickTest<tags::GPU::CUDA, float> cuda_dense_vector_Slice_3scaled_sum_quick_test_float("float");
#endif
#ifdef HONEI_CELL
DenseVectorSlice3ScaledSumQuickTest<tags::Cell, float> cell_dense_vector_Slice_3scaled_sum_quick_test_float("Cell float");
DenseVectorSlice3ScaledSumQuickTest<tags::Cell, double> cell_dense_vector_Slice_3scaled_sum_quick_test_double("Cell double");
#endif
*/

template <typename Tag_, typename DataType_>
class DenseVectorSparseVectorScaledSumTest :
    public BaseTest
{
    public:
        DenseVectorSparseVectorScaledSumTest(const std::string & type) :
            BaseTest("dense_vector_sparse_vector_scaled_sum_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size, DataType_(0));
                for (typename DenseVector<DataType_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0) *i = DataType_(2);
                }
                SparseVector<DataType_> sv2(size, size / 8 + 1);
                for (typename SparseVector<DataType_>::ElementIterator i(sv2.begin_elements()), i_end(sv2.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 7 == 0) *i = DataType_(3);
                }

                DataType_ scal(2);

                DenseVector<DataType_> sum2(size, DataType_(0));
                for (typename DenseVector<DataType_>::ElementIterator i(sum2.begin_elements()), i_end(sum2.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 7 == 0) *i = DataType_(6);
                    if (i.index() % 10 == 0) *i = DataType_(2);
                    if (i.index() % 10 == 0 && i.index() % 7 == 0) *i = DataType_(8);
                }

                ScaledSum<Tag_>::value(dv1, sv2, scal);

                TEST_CHECK_EQUAL(dv1, sum2);
            }

            DenseVector<DataType_> sv00(1);
            SparseVector<DataType_> sv01(2, 1);
            DataType_ scal00(DataType_(2));

            TEST_CHECK_THROWS(ScaledSum<Tag_>::value(sv00, sv01, scal00), VectorSizeDoesNotMatch);
        }
};
DenseVectorSparseVectorScaledSumTest<tags::CPU, float> dense_vector_sparse_vector_scaled_sum_test_float("float");
DenseVectorSparseVectorScaledSumTest<tags::CPU, double> dense_vector_sparse_vector_scaled__sum_test_double("double");

template <typename Tag_, typename DataType_>
class DenseVectorSparseVectorScaledSumQuickTest :
    public QuickTest
{
    public:
        DenseVectorSparseVectorScaledSumQuickTest(const std::string & type) :
            QuickTest("dense_vector_sparse_vector_scaled_sum_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size (22);
            DenseVector<DataType_> dv1(size, DataType_(0));
            for (typename DenseVector<DataType_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 10 == 0) *i = DataType_(2);
            }
            SparseVector<DataType_> sv2(size, size / 8 + 1);
            for (typename SparseVector<DataType_>::ElementIterator i(sv2.begin_elements()), i_end(sv2.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 7 == 0) *i = DataType_(3);
            }

            DataType_ scal(2);

            DenseVector<DataType_> sum2(size, DataType_(0));
            for (typename DenseVector<DataType_>::ElementIterator i(sum2.begin_elements()), i_end(sum2.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 7 == 0) *i = DataType_(6);
                if (i.index() % 10 == 0) *i = DataType_(2);
                if (i.index() % 10 == 0 && i.index() % 7 == 0) *i = DataType_(8);
            }

            ScaledSum<Tag_>::value(dv1, sv2, scal);

            TEST_CHECK_EQUAL(dv1, sum2);

            DenseVector<DataType_> sv00(1);
            SparseVector<DataType_> sv01(2, 1);
            DataType_ scal00(DataType_(2));

            TEST_CHECK_THROWS(ScaledSum<Tag_>::value(sv00, sv01, scal00), VectorSizeDoesNotMatch);
        }
};
DenseVectorSparseVectorScaledSumQuickTest<tags::CPU, float> dense_vector_sparse_vector_scaled_sum_quick_test_float("float");
DenseVectorSparseVectorScaledSumQuickTest<tags::CPU, double> dense_vector_sparse_vector_scaled__sum_quick_test_double("double");

template <typename Tag_, typename DataType_>
class SparseVectorScaledSumTest :
    public BaseTest
{
    public:
        SparseVectorScaledSumTest(const std::string & type) :
            BaseTest("sparse_vector_scaled_sum_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                SparseVector<DataType_> sv1(size, size / 7 + 1);
                for (typename SparseVector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0) *i = DataType_(2);
                }
                SparseVector<DataType_> sv2(size, size / 8 + 1);
                for (typename SparseVector<DataType_>::ElementIterator i(sv2.begin_elements()), i_end(sv2.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 7 == 0) *i = DataType_(3);
                }

                DataType_ scal(2);

                SparseVector<DataType_> sum2(size, size / 5 + 1);
                for (typename SparseVector<DataType_>::ElementIterator i(sum2.begin_elements()), i_end(sum2.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 7 == 0) *i = DataType_(6);
                    if (i.index() % 10 == 0) *i = DataType_(2);
                    if (i.index() % 10 == 0 && i.index() % 7 == 0) *i = DataType_(8);
                }

                ScaledSum<Tag_>::value(sv1, sv2, scal);

                TEST_CHECK_EQUAL(sv1, sum2);
            }

            SparseVector<DataType_> sv00(1, 1);
            SparseVector<DataType_> sv01(2, 1);
            DataType_ scal00(DataType_(2));

            TEST_CHECK_THROWS(ScaledSum<Tag_>::value(sv00, sv01, scal00), VectorSizeDoesNotMatch);
        }
};
SparseVectorScaledSumTest<tags::CPU, float> sparse_vector_scaled_sum_test_float("float");
SparseVectorScaledSumTest<tags::CPU, double> sparse_vector_scaled__sum_test_double("double");

template <typename Tag_, typename DataType_>
class SparseVectorScaledSumQuickTest :
    public QuickTest
{
    public:
        SparseVectorScaledSumQuickTest(const std::string & type) :
            QuickTest("sparse_vector_scaled_sum_quick_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            unsigned long size(15);
            SparseVector<DataType_> sv1(size, size / 7 + 1);
            for (typename SparseVector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 10 == 0) *i = DataType_(2);
            }
            SparseVector<DataType_> sv2(size, size / 8 + 1);
            for (typename SparseVector<DataType_>::ElementIterator i(sv2.begin_elements()), i_end(sv2.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 7 == 0) *i = DataType_(3);
            }
            DataType_ scal(2);

            SparseVector<DataType_> sum2(size, size / 5 + 1);
            for (typename SparseVector<DataType_>::ElementIterator i(sum2.begin_elements()), i_end(sum2.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 7 == 0) *i = DataType_(6);
                if (i.index() % 10 == 0) *i = DataType_(2);
                if (i.index() % 10 == 0 && i.index() % 7 == 0) *i = DataType_(8);
            }

            ScaledSum<Tag_>::value(sv1, sv2, scal);

            TEST_CHECK_EQUAL(sv1, sum2);

            SparseVector<DataType_> sv00(1, 1);
            SparseVector<DataType_> sv01(2, 1);
            DataType_ scal00(DataType_(2));

            TEST_CHECK_THROWS(ScaledSum<Tag_>::value(sv00, sv01, scal00), VectorSizeDoesNotMatch);
        }
};
SparseVectorScaledSumQuickTest<tags::CPU, float> sparse_vector_scaled_sum_quick_test_float("float");
SparseVectorScaledSumQuickTest<tags::CPU, double> sparse_vector_scaled_sum_quick_test_double("double");
