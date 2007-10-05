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

#include <libla/dense_vector.hh>
#include <libla/norm.hh>
#include <libla/scaled_sum.hh>
#include <libla/sparse_vector.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>

using namespace honei;
using  namespace tests;

template <typename Tag_, typename DataType_>
class DenseVectorScaledSumTest :
    public BaseTest
{
    public:
        DenseVectorScaledSumTest(const std::string & type) :
            BaseTest("dense_vector_scaled_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 10) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size, DataType_(2));
                DenseVector<DataType_> dv2(size, DataType_(3));
                DataType_ scal(DataType_(2));

                DenseVector<DataType_> sum1(ScaledSum<Tag_>::value(dv1, dv2, scal));
                DataType_ v1(Norm<vnt_l_one>::value(sum1));
                TEST_CHECK_EQUAL(v1, 8 * size);
            }

            DenseVector<DataType_> dv00(1, DataType_(1));
            DenseVector<DataType_> dv01(2, DataType_(1));
            DataType_ scal00(DataType_(2));

            TEST_CHECK_THROWS(ScaledSum<Tag_>::value(dv00, dv01, scal00), VectorSizeDoesNotMatch);
        }
};
DenseVectorScaledSumTest<tags::CPU, float> dense_vector_scaled_sum_test_float("float");
DenseVectorScaledSumTest<tags::CPU, double> dense_vector_scaled__sum_test_double("double");
#ifdef HONEI_SSE
DenseVectorScaledSumTest<tags::CPU::SSE, float> sse_dense_vector_scaled_sum_test_float("SSE float");
DenseVectorScaledSumTest<tags::CPU::SSE, double> sse_dense_vector_scaled_sum_test_double("SSE double");
#endif

template <typename Tag_, typename DataType_>
class DenseVectorScaledSumQuickTest :
    public QuickTest
{
    public:
        DenseVectorScaledSumQuickTest(const std::string & type) :
            QuickTest("dense_vector_scaled_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(65);
            DenseVector<DataType_> dv1(size, DataType_(2));
            DenseVector<DataType_> dv2(size, DataType_(3));
            DataType_ scal(DataType_(2));
            DenseVector<DataType_> dv3(size, DataType_(8));

            DenseVector<DataType_> sum1(ScaledSum<Tag_>::value(dv1, dv2, scal));
            TEST_CHECK_EQUAL(sum1, dv3);

            DenseVector<DataType_> dv00(1, DataType_(1));
            DenseVector<DataType_> dv01(2, DataType_(1));
            DataType_ scal00(DataType_(2));

            TEST_CHECK_THROWS(ScaledSum<Tag_>::value(dv00, dv01, scal00), VectorSizeDoesNotMatch);
        }
};
DenseVectorScaledSumQuickTest<tags::CPU, float> dense_vector_scaled_sum_quick_test_float("float");
DenseVectorScaledSumQuickTest<tags::CPU, double> dense_vector_scaled__sum_quick_test_double("double");
#ifdef HONEI_SSE
DenseVectorScaledSumQuickTest<tags::CPU::SSE, float> sse_dense_vector_scaled_sum_quick_test_float("SSE float");
DenseVectorScaledSumQuickTest<tags::CPU::SSE, double> sse_dense_vector_scaled_sum_quick_test_double("SSE double");
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

                DenseVector<DataType_> sum1(ScaledSum<>::value(dv1, sv2, scal));

                TEST_CHECK_EQUAL(sum1, sum2);
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

            DenseVector<DataType_> sum1(ScaledSum<>::value(dv1, sv2, scal));

            TEST_CHECK_EQUAL(sum1, sum2);

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

                SparseVector<DataType_> sum1(ScaledSum<>::value(sv1, sv2, scal));

                TEST_CHECK_EQUAL(sum1, sum2);
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

            SparseVector<DataType_> sum1(ScaledSum<>::value(sv1, sv2, scal));

            TEST_CHECK_EQUAL(sum1, sum2);

            SparseVector<DataType_> sv00(1, 1);
            SparseVector<DataType_> sv01(2, 1);
            DataType_ scal00(DataType_(2));

            TEST_CHECK_THROWS(ScaledSum<>::value(sv00, sv01, scal00), VectorSizeDoesNotMatch);
        }
};
SparseVectorScaledSumQuickTest<float> sparse_vector_scaled_sum_quick_test_float("float");
SparseVectorScaledSumQuickTest<double> sparse_vector_scaled__sum_quick_test_double("double");
