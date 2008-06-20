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
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>

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
                MemoryArbiter::instance()->write<tags::CPU>(dv1.memid(), dv1.address(), dv1.size() * sizeof(DataType_));
                MemoryArbiter::instance()->release_write<tags::CPU>(dv1.memid());
                MemoryArbiter::instance()->write<tags::CPU>(dv2.memid(), dv2.address(), dv2.size() * sizeof(DataType_));
                MemoryArbiter::instance()->release_write<tags::CPU>(dv2.memid());
                DataType_ scal(DataType_(2));

                ScaledSum<Tag_>::value(dv1, dv2, scal);
                MemoryArbiter::instance()->read<tags::CPU>(dv1.memid(), dv1.address(), dv1.size() * sizeof(DataType_));
                DenseVector<DataType_> dv3(size, DataType_(8));
                MemoryArbiter::instance()->release_read<tags::CPU>(dv1.memid());

                TEST_CHECK_EQUAL(dv1, dv3);
            }

            DenseVector<DataType_> dv00(1, DataType_(1));
            DenseVector<DataType_> dv01(2, DataType_(1));
            DataType_ scal00(DataType_(2));

            TEST_CHECK_THROWS(ScaledSum<Tag_>::value(dv00, dv01, scal00), VectorSizeDoesNotMatch);
        }
};
DenseVectorScaledSumTest<tags::CPU, float> dense_vector_scaled_sum_test_float("float");
DenseVectorScaledSumTest<tags::CPU, double> dense_vector_scaled_sum_test_double("double");
DenseVectorScaledSumTest<tags::CPU::MultiCore, float> mc_dense_vector_scaled_sum_test_float("float");
DenseVectorScaledSumTest<tags::CPU::MultiCore, double> mc_dense_vector_scaled_sum_test_double("double");
#ifdef HONEI_SSE
DenseVectorScaledSumTest<tags::CPU::SSE, float> sse_dense_vector_scaled_sum_test_float("SSE float");
DenseVectorScaledSumTest<tags::CPU::SSE, double> sse_dense_vector_scaled_sum_test_double("SSE double");
DenseVectorScaledSumTest<tags::CPU::MultiCore::SSE, float>
    mc_sse_dense_vector_scaled_sum_test_float("MC SSE float");
DenseVectorScaledSumTest<tags::CPU::MultiCore::SSE, double>
    mc_sse_dense_vector_scaled_sum_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
DenseVectorScaledSumTest<tags::GPU::CUDA, float> cuda_dense_vector_scaled_sum_test_float("float");
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
            MemoryArbiter::instance()->write<tags::CPU>(dv1.memid(), dv1.address(), dv1.size() * sizeof(DataType_));
            MemoryArbiter::instance()->release_write<tags::CPU>(dv1.memid());
            MemoryArbiter::instance()->write<tags::CPU>(dv2.memid(), dv2.address(), dv2.size() * sizeof(DataType_));
            MemoryArbiter::instance()->release_write<tags::CPU>(dv2.memid());
            DataType_ scal(DataType_(2));

            ScaledSum<Tag_>::value(dv1, dv2, scal);
            DenseVector<DataType_> dv3(size, DataType_(8));

            MemoryArbiter::instance()->read<tags::CPU>(dv1.memid(), dv1.address(), dv1.size() * sizeof(DataType_));
            TEST_CHECK_EQUAL(dv1, dv3);
            MemoryArbiter::instance()->release_read<tags::CPU>(dv1.memid());

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
#ifdef HONEI_SSE
DenseVectorScaledSumQuickTest<tags::CPU::SSE, float> sse_dense_vector_scaled_sum_quick_test_float("SSE float");
DenseVectorScaledSumQuickTest<tags::CPU::SSE, double> sse_dense_vector_scaled_sum_quick_test_double("SSE double");
DenseVectorScaledSumQuickTest<tags::CPU::MultiCore::SSE, float> mc_sse_dense_vector_scaled_sum_quick_test_float("MC SSE float");
DenseVectorScaledSumQuickTest<tags::CPU::MultiCore::SSE, double> mc_sse_dense_vector_scaled_sum_quick_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
DenseVectorScaledSumQuickTest<tags::GPU::CUDA, float> cuda_dense_vector_scaled_sum_quick_test_float("float");
#endif
#ifdef HONEI_CELL
DenseVectorScaledSumQuickTest<tags::Cell, float> cell_dense_vector_scaled_sum_quick_test_float("Cell float");
DenseVectorScaledSumQuickTest<tags::Cell, double> cell_dense_vector_scaled_sum_quick_test_double("Cell double");
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
                for (int i(0) ; i < 4 ; i++)
                {
                    for (int j(0) ; j < 4 ; j++)
                    {
                        DenseVector<DataType_> dv1_source(size * 2, DataType_(2));
                        DenseVectorRange<DataType_> dv1(dv1_source, size, i);

                        DenseVector<DataType_> dv2_source(size * 2, DataType_(3));
                        DenseVectorRange<DataType_> dv2(dv2_source, size, j);
                        MemoryArbiter::instance()->write<tags::CPU>(dv1.memid(), dv1.address(), dv1.size() * sizeof(DataType_));
                        MemoryArbiter::instance()->release_write<tags::CPU>(dv1.memid());
                        MemoryArbiter::instance()->write<tags::CPU>(dv2.memid(), dv2.address(), dv2.size() * sizeof(DataType_));
                        MemoryArbiter::instance()->release_write<tags::CPU>(dv2.memid());
                        DataType_ scal(DataType_(2));

                        ScaledSum<Tag_>::value(dv1, dv2, scal);
                        MemoryArbiter::instance()->read<tags::CPU>(dv1.memid(), dv1.address(), dv1.size() * sizeof(DataType_));
                        DenseVector<DataType_> dv3(size, DataType_(8));

                        TEST_CHECK_EQUAL(dv1, dv3);
                        MemoryArbiter::instance()->release_read<tags::CPU>(dv1.memid());
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
#ifdef HONEI_SSE
DenseVectorRangeScaledSumTest<tags::CPU::SSE, float> sse_dense_vector_range_scaled_sum_test_float("SSE float");
DenseVectorRangeScaledSumTest<tags::CPU::SSE, double> sse_dense_vector_range_scaled_sum_test_double("SSE double");
DenseVectorRangeScaledSumTest<tags::CPU::MultiCore::SSE, float> mc_sse_dense_vector_range_scaled_sum_test_float("MC SSE float");
DenseVectorRangeScaledSumTest<tags::CPU::MultiCore::SSE, double> mc_sse_dense_vector_range_scaled_sum_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
DenseVectorRangeScaledSumTest<tags::GPU::CUDA, float> cuda_dense_vector_range_scaled_sum_test_float("float");
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
            MemoryArbiter::instance()->write<tags::CPU>(dv1.memid(), dv1.address(), dv1.size() * sizeof(DataType_));
            MemoryArbiter::instance()->release_write<tags::CPU>(dv1.memid());
            MemoryArbiter::instance()->write<tags::CPU>(dv2.memid(), dv2.address(), dv2.size() * sizeof(DataType_));
            MemoryArbiter::instance()->release_write<tags::CPU>(dv2.memid());

            DataType_ scal(DataType_(2));
            DenseVector<DataType_> dv3(size, DataType_(8));
            ScaledSum<Tag_>::value(dv1, dv2, scal);
            MemoryArbiter::instance()->read<tags::CPU>(dv1.memid(), dv1.address(), dv1.size() * sizeof(DataType_));
            TEST_CHECK_EQUAL(dv1, dv3);
            MemoryArbiter::instance()->release_read<tags::CPU>(dv1.memid());

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
#ifdef HONEI_SSE
DenseVectorRangeScaledSumQuickTest<tags::CPU::SSE, float> sse_dense_vector_range_scaled_sum_quick_test_float("SSE float");
DenseVectorRangeScaledSumQuickTest<tags::CPU::SSE, double> sse_dense_vector_range_scaled_sum_quick_test_double("SSE double");
DenseVectorRangeScaledSumQuickTest<tags::CPU::MultiCore::SSE, float> mc_sse_dense_vector_range_scaled_sum_quick_test_float("MC SSE float");
DenseVectorRangeScaledSumQuickTest<tags::CPU::MultiCore::SSE, double> mc_sse_dense_vector_range_scaled_sum_quick_test_double("MC SSE double");
#endif
#ifdef HONEI_CUDA
DenseVectorRangeScaledSumQuickTest<tags::GPU::CUDA, float> cuda_dense_vector_range_scaled_sum_quick_test_float("float");
#endif
#ifdef HONEI_CELL
DenseVectorRangeScaledSumQuickTest<tags::Cell, float> cell_dense_vector_range_scaled_sum_quick_test_float("Cell float");
DenseVectorRangeScaledSumQuickTest<tags::Cell, double> cell_dense_vector_range_scaled_sum_quick_test_double("Cell double");
#endif

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
                DataType_ v1(Norm<vnt_l_one>::value(sumr));

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
            DataType_ v1(Norm<vnt_l_one>::value(sumr));
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

template <typename DataType_>
class DenseVectorSparseVectorScaledSumTest :
    public BaseTest
{
    public:
        DenseVectorSparseVectorScaledSumTest(const std::string & type) :
            BaseTest("dense_vector_sparse_vector_scaled_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size, DataType_(0));
                for (typename Vector<DataType_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0) *i = DataType_(2);
                }
                SparseVector<DataType_> sv2(size, size / 8 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv2.begin_elements()), i_end(sv2.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 7 == 0) *i = DataType_(3);
                }

                DataType_ scal(2);

                DenseVector<DataType_> sum2(size, DataType_(0));
                for (typename Vector<DataType_>::ElementIterator i(sum2.begin_elements()), i_end(sum2.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 7 == 0) *i = DataType_(6);
                    if (i.index() % 10 == 0) *i = DataType_(2);
                    if (i.index() % 10 == 0 && i.index() % 7 == 0) *i = DataType_(8);
                }

                ScaledSum<>::value(dv1, sv2, scal);

                TEST_CHECK_EQUAL(dv1, sum2);
            }

            DenseVector<DataType_> sv00(1);
            SparseVector<DataType_> sv01(2, 1);
            DataType_ scal00(DataType_(2));

            TEST_CHECK_THROWS(ScaledSum<>::value(sv00, sv01, scal00), VectorSizeDoesNotMatch);
        }
};
DenseVectorSparseVectorScaledSumTest<float> dense_vector_sparse_vector_scaled_sum_test_float("float");
DenseVectorSparseVectorScaledSumTest<double> dense_vector_sparse_vector_scaled__sum_test_double("double");

template <typename DataType_>
class DenseVectorSparseVectorScaledSumQuickTest :
    public QuickTest
{
    public:
        DenseVectorSparseVectorScaledSumQuickTest(const std::string & type) :
            QuickTest("dense_vector_sparse_vector_scaled_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size (22);
            DenseVector<DataType_> dv1(size, DataType_(0));
            for (typename Vector<DataType_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 10 == 0) *i = DataType_(2);
            }
            SparseVector<DataType_> sv2(size, size / 8 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv2.begin_elements()), i_end(sv2.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 7 == 0) *i = DataType_(3);
            }

            DataType_ scal(2);

            DenseVector<DataType_> sum2(size, DataType_(0));
            for (typename Vector<DataType_>::ElementIterator i(sum2.begin_elements()), i_end(sum2.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 7 == 0) *i = DataType_(6);
                if (i.index() % 10 == 0) *i = DataType_(2);
                if (i.index() % 10 == 0 && i.index() % 7 == 0) *i = DataType_(8);
            }

            ScaledSum<>::value(dv1, sv2, scal);

            TEST_CHECK_EQUAL(dv1, sum2);

            DenseVector<DataType_> sv00(1);
            SparseVector<DataType_> sv01(2, 1);
            DataType_ scal00(DataType_(2));

            TEST_CHECK_THROWS(ScaledSum<>::value(sv00, sv01, scal00), VectorSizeDoesNotMatch);
        }
};
DenseVectorSparseVectorScaledSumQuickTest<float> dense_vector_sparse_vector_scaled_sum_quick_test_float("float");
DenseVectorSparseVectorScaledSumQuickTest<double> dense_vector_sparse_vector_scaled__sum_quick_test_double("double");

template <typename DataType_>
class SparseVectorScaledSumTest :
    public BaseTest
{
    public:
        SparseVectorScaledSumTest(const std::string & type) :
            BaseTest("sparse_vector_scaled_sum_test<" + type + ">")
    {
    }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                SparseVector<DataType_> sv1(size, size / 7 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0) *i = DataType_(2);
                }
                SparseVector<DataType_> sv2(size, size / 8 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv2.begin_elements()), i_end(sv2.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 7 == 0) *i = DataType_(3);
                }

                DataType_ scal(2);

                SparseVector<DataType_> sum2(size, size / 5 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sum2.begin_elements()), i_end(sum2.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 7 == 0) *i = DataType_(6);
                    if (i.index() % 10 == 0) *i = DataType_(2);
                    if (i.index() % 10 == 0 && i.index() % 7 == 0) *i = DataType_(8);
                }

                ScaledSum<>::value(sv1, sv2, scal);

                TEST_CHECK_EQUAL(sv1, sum2);
            }

            SparseVector<DataType_> sv00(1, 1);
            SparseVector<DataType_> sv01(2, 1);
            DataType_ scal00(DataType_(2));

            TEST_CHECK_THROWS(ScaledSum<>::value(sv00, sv01, scal00), VectorSizeDoesNotMatch);
        }
};
SparseVectorScaledSumTest<float> sparse_vector_scaled_sum_test_float("float");
SparseVectorScaledSumTest<double> sparse_vector_scaled__sum_test_double("double");

template <typename DataType_>
class SparseVectorScaledSumQuickTest :
    public QuickTest
{
    public:
        SparseVectorScaledSumQuickTest(const std::string & type) :
            QuickTest("sparse_vector_scaled_sum_quick_test<" + type + ">")
    {
    }

        virtual void run() const
        {
            unsigned long size(15);
            SparseVector<DataType_> sv1(size, size / 7 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 10 == 0) *i = DataType_(2);
            }
            SparseVector<DataType_> sv2(size, size / 8 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv2.begin_elements()), i_end(sv2.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 7 == 0) *i = DataType_(3);
            }
            DataType_ scal(2);

            SparseVector<DataType_> sum2(size, size / 5 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sum2.begin_elements()), i_end(sum2.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 7 == 0) *i = DataType_(6);
                if (i.index() % 10 == 0) *i = DataType_(2);
                if (i.index() % 10 == 0 && i.index() % 7 == 0) *i = DataType_(8);
            }

            ScaledSum<>::value(sv1, sv2, scal);

            TEST_CHECK_EQUAL(sv1, sum2);

            SparseVector<DataType_> sv00(1, 1);
            SparseVector<DataType_> sv01(2, 1);
            DataType_ scal00(DataType_(2));

            TEST_CHECK_THROWS(ScaledSum<>::value(sv00, sv01, scal00), VectorSizeDoesNotMatch);
        }
};
SparseVectorScaledSumQuickTest<float> sparse_vector_scaled_sum_quick_test_float("float");
SparseVectorScaledSumQuickTest<double> sparse_vector_scaled__sum_quick_test_double("double");
