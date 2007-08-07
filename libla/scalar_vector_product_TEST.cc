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
#include <libla/scalar_vector_product.hh>
#include <libla/vector_norm.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
class ScalarDenseVectorProductTest :
    public BaseTest
{
    public:
        ScalarDenseVectorProductTest(const std::string & type) :
            BaseTest("scalar_dense_vector_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size, DataType_(3));

                DenseVector<DataType_> prod1(ScalarVectorProduct<DataType_>::value(DataType_(2), dv1));
                DataType_ v1(VectorNorm<DataType_, vnt_l_one>::value(prod1));
                TEST_CHECK_EQUAL(v1, 6 * size);
            }
        }
};

ScalarDenseVectorProductTest<float> scalar_dense_vector_product_test_float("float");
ScalarDenseVectorProductTest<double> scalar_dense_vector_product_test_double("double");

template <typename DataType_>
class ScalarDenseVectorProductQuickTest :
    public QuickTest
{
    public:
        ScalarDenseVectorProductQuickTest(const std::string & type) :
            QuickTest("scalar_dense_vector_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseVector<DataType_> dv1(size, DataType_(3));

            DenseVector<DataType_> prod1(ScalarVectorProduct<DataType_>::value(DataType_(2), dv1));
            DataType_ v1(VectorNorm<DataType_, vnt_l_one>::value(prod1));
            TEST_CHECK_EQUAL(v1, 6 * size);
        }
};
ScalarDenseVectorProductQuickTest<float>  scalar_dense_vector_product_quick_test_float("float");
ScalarDenseVectorProductQuickTest<double> scalar_dense_vector_product_quick_test_double("double");

template <typename DataType_>
class ScalarSparseVectorProductTest :
    public BaseTest
{
    public:
        ScalarSparseVectorProductTest(const std::string & type) :
            BaseTest("scalar_sparse_vector_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(11) ; size < (1 << 14) ; size <<= 1)
            {
                SparseVector<DataType_> sv1(size, size / 8 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0) *i = 3;
                }
                SparseVector<DataType_> prod1(ScalarVectorProduct<DataType_>::value(static_cast<DataType_>(2), sv1));
                DataType_ v1(VectorNorm<DataType_, vnt_l_one>::value(prod1));
                TEST_CHECK_EQUAL(v1, 6 * (size / 10 + 1));
            }
        }
};

ScalarSparseVectorProductTest<float> scalar_sparse_vector_product_test_float("float");
ScalarSparseVectorProductTest<double> scalar_sparse_vector_product_test_double("double");

template <typename DataType_>
class ScalarSparseVectorProductQuickTest :
    public QuickTest
{
    public:
        ScalarSparseVectorProductQuickTest(const std::string & type) :
            QuickTest("scalar_sparse_vector_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(111);
            SparseVector<DataType_> sv1(size, size / 8 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 10 == 0) *i = 3;
            }
            SparseVector<DataType_> prod1(ScalarVectorProduct<DataType_>::value(static_cast<DataType_>(2), sv1));
            DataType_ v1(VectorNorm<DataType_, vnt_l_one>::value(prod1));
            TEST_CHECK_EQUAL(v1, 6 * (size / 10 + 1));
         
        }
};

ScalarSparseVectorProductQuickTest<float> scalar_sparse_vector_product_quick_test_float("float");
ScalarSparseVectorProductQuickTest<double> scalar_sparse_vector_product_quick_test_double("double");
