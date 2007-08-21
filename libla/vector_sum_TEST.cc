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
#include <libla/vector_sum.hh>
#include <libla/vector_norm.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
class DenseVectorSumTest :
    public BaseTest
{
    public:
        DenseVectorSumTest(const std::string & type) :
            BaseTest("dense_vector_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size, DataType_(1)), dv2(size, DataType_(2));

                VectorSum<>::value(dv1, dv2);

                DataType_ v1(VectorNorm<DataType_, vnt_l_one>::value(dv1));
                TEST_CHECK_EQUAL_WITHIN_EPS(v1, 3 * size, std::numeric_limits<DataType_>::epsilon());

                DenseVector<DataType_> dv3(size + 1), dv4(size);
                TEST_CHECK_THROWS(VectorSum<>::value(dv3, dv4), VectorSizeDoesNotMatch);
            }
        }
};

DenseVectorSumTest<float> dense_vector_sum_test_float("float");
DenseVectorSumTest<double> dense_vector_sum_test_double("double");

template <typename DataType_>
class DenseVectorSumQuickTest :
    public QuickTest
{
    public:
        DenseVectorSumQuickTest(const std::string & type) :
            QuickTest("dense_vector_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseVector<DataType_> dv1(size, DataType_(1)), dv2(size, DataType_(2));

            VectorSum<>::value(dv1, dv2);
            DataType_ v1(VectorNorm<DataType_, vnt_l_one>::value(dv1));
            TEST_CHECK_EQUAL_WITHIN_EPS(v1, 3 * size, std::numeric_limits<DataType_>::epsilon());

            DenseVector<DataType_> dv3(6, DataType_(1)), dv4(4, DataType_(1));
            TEST_CHECK_THROWS(VectorSum<>::value(dv3, dv4), VectorSizeDoesNotMatch);
        }
};
DenseVectorSumQuickTest<float>  dense_vector_sum_quick_test_float("float");
DenseVectorSumQuickTest<double> dense_vector_sum_quick_test_double("double");

template <typename DataType_>
class SparseVectorSumTest :
    public BaseTest
{
    public:
        SparseVectorSumTest(const std::string & type) :
            BaseTest("sparse_vector_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                SparseVector<DataType_> sv1(size, size / 8 + 1), sv2(size, size / 7 + 1), sv3(size, size / 8 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()),
                    j(sv2.begin_elements()), k(sv3.begin_elements()) ; i != i_end ; ++i, ++j, ++k)
                {
                    if (i.index() % 10 == 0) 
                    {
                        *i = DataType_(1);
                        *k = DataType_(1);
                    }
                    if (i.index() % 7 == 0) 
                    {
                        *j = DataType_(2);
                        *k = DataType_(2);
                    }
                    if (i.index() % 7 == 0 && i.index() % 10 == 0) 
                    {
                        *k = DataType_(3);
                    }                                        
                }

                VectorSum<>::value(sv1, sv2);

                TEST_CHECK_EQUAL(sv1, sv3);
            }
            
            SparseVector<DataType_> sv01(2, 1), sv02(1, 1);
            TEST_CHECK_THROWS(VectorSum<>::value(sv01, sv02), VectorSizeDoesNotMatch);            
        }
};
SparseVectorSumTest<float> sparse_vector_sum_test_float("float");
SparseVectorSumTest<double> sparse_vector_sum_test_double("double");

template <typename DataType_>
class SparseVectorSumQuickTest :
    public QuickTest
{
    public:
        SparseVectorSumQuickTest(const std::string & type) :
            QuickTest("sparse_vector_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(20);
            SparseVector<DataType_> sv1(size, size / 8 + 1), sv2(size, size / 7 + 1), sv3(size, size / 8 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()),
                j(sv2.begin_elements()), k(sv3.begin_elements()) ; i != i_end ; ++i, ++j, ++k)
            {
                if (i.index() % 10 == 0) 
                {
                    *i = DataType_(1);
                    *k = DataType_(1);
                }
                if (i.index() % 7 == 0) 
                {
                    *j = DataType_(2);
                    *k = DataType_(2);
                }
                if (i.index() % 7 == 0 && i.index() % 10 == 0) 
                {
                    *k = DataType_(3);
                }                                        
            }

            VectorSum<>::value(sv1, sv2);

            TEST_CHECK_EQUAL(sv1, sv3);

            SparseVector<DataType_> sv01(21, 1), sv02(20, 1);
            TEST_CHECK_THROWS(VectorSum<>::value(sv01, sv02), VectorSizeDoesNotMatch);
        }
};
SparseVectorSumQuickTest<float> sparse_vector_sum_quick_test_float("float");
SparseVectorSumQuickTest<double> sparse_vector_sum_quick_test_double("double");
