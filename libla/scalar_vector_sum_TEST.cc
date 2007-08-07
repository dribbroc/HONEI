/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
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
#include <libla/sparse_vector.hh>
#include <libla/scalar_vector_sum.hh>
#include <libla/vector_norm.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
class ScalarDenseVectorSumTest :
    public BaseTest
{
    public:
        ScalarDenseVectorSumTest(const std::string & type) :
            BaseTest("scalar_dense_vector_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size, DataType_(3));

                DenseVector<DataType_> sum1(ScalarVectorSum<DataType_>::value(DataType_(2), dv1));
                DataType_ v1(VectorNorm<DataType_, vnt_l_one>::value(sum1));
                TEST_CHECK_EQUAL(v1, 5 * size);
            }
        }
};

ScalarDenseVectorSumTest<float> scalar_dense_vector_sum_test_float("float");
ScalarDenseVectorSumTest<double> scalar_dense_vector_sum_test_double("double");

template <typename DataType_>
class ScalarDenseVectorSumQuickTest :
    public QuickTest
{
    public:
        ScalarDenseVectorSumQuickTest(const std::string & type) :
            QuickTest("scalar_dense_vector_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseVector<DataType_> dv1(size, DataType_(3));

            DenseVector<DataType_> sum1(ScalarVectorSum<DataType_>::value(DataType_(2), dv1));
            DataType_ v1(VectorNorm<DataType_, vnt_l_one>::value(sum1));
            TEST_CHECK_EQUAL(v1, 5 * size);
        }
};
ScalarDenseVectorSumQuickTest<float>  scalar_dense_vector_sum_quick_test_float("float");
ScalarDenseVectorSumQuickTest<double> scalar_dense_vector_sum_quick_test_double("double");

template <typename DataType_>
class ScalarSparseVectorSumTest :
    public BaseTest
{
    public:
        ScalarSparseVectorSumTest(const std::string & type) :
            BaseTest("scalar_sparse_vector_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                SparseVector<DataType_> sv1(size, size / 8 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0) *i = 3;
                }

                SparseVector<DataType_> sum1(ScalarVectorSum<DataType_>::value(DataType_(2), sv1));
                DataType_ v1(VectorNorm<DataType_, vnt_l_one>::value(sum1));
                TEST_CHECK_EQUAL(v1, 2 * size + 3 * (size / 10 + 1));
            }
        }
};
ScalarSparseVectorSumTest<float> scalar_sparse_vector_sum_test_float("float");
ScalarSparseVectorSumTest<double> scalar_sparse_vector_sum_test_double("double");

template <typename DataType_>
class ScalarSparseVectorSumQuickTest :
    public QuickTest
{
    public:
        ScalarSparseVectorSumQuickTest(const std::string & type) :
            QuickTest("scalar_sparse_vector_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            SparseVector<DataType_> sv1(size, size / 8 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 10 == 0) *i = 3;
            }

            SparseVector<DataType_> sum1(ScalarVectorSum<DataType_>::value(DataType_(2), sv1));
            DataType_ v1(VectorNorm<DataType_, vnt_l_one>::value(sum1));
            TEST_CHECK_EQUAL(v1, 2 * size + 3 * (size / 10 + 1));
        }
};
ScalarSparseVectorSumQuickTest<float> scalar_sparse_vector_sum_quick_test_float("float");
ScalarSparseVectorSumQuickTest<double> scalar_sparse_vector_sum_quick_test_double("double");
