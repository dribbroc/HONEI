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
#include <libla/vector_elementwise_product.hh>
#include <libla/vector_error.cc>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using namespace tests;

template <typename DataType_>
class DenseVectorElementwiseProductTest :
    public BaseTest
{
    public:
        DenseVectorElementwiseProductTest(const std::string & type) :
            BaseTest("dense_vector_elementwise_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size, DataType_(2)), dv2(size, DataType_(3)),
                    dv3(size, DataType_(6 * size));
                DenseVector<DataType_> prod(VectorElementwiseProduct<DataType_>::value(dv1, dv2));

                TEST_CHECK_EQUAL(prod, dv3);
            }

            DenseVector<DataType_> dv01(3, static_cast<DataType_>(1)), dv02(4, static_cast<DataType_>(1));

            TEST_CHECK_THROWS(VectorElementwiseProduct<DataType_>::value(dv02, dv01), VectorSizeDoesNotMatch);
        }
};

DenseVectorElementwiseProductTest<float> dense_vector_elementwise_product_test_float("float");
DenseVectorElementwiseProductTest<double> dense_vector_elementwise_product_test_double("double");

template <typename DataType_>
class DenseVectorElementwiseProductQuickTest :
    public QuickTest
{
    public:
        DenseVectorElementwiseProductQuickTest(const std::string & type) :
            QuickTest("dense_vector_elementwise_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseVector<DataType_> dv1(size, DataType_(2)), dv2(size, DataType_(3)),
                dv3(size, DataType_(6));
            DenseVector<DataType_> prod(VectorElementwiseProduct<DataType_>::value(dv1, dv2));

            TEST_CHECK_EQUAL(prod, dv3);


            DenseVector<DataType_> dv01(3, static_cast<DataType_>(1)), dv02(4, static_cast<DataType_>(1));

            TEST_CHECK_THROWS(VectorElementwiseProduct<DataType_>::value(dv02, dv01), VectorSizeDoesNotMatch);

        }
};
DenseVectorElementwiseProductQuickTest<float> dense_vector_elementwise_product_quick_test_float("float");
DenseVectorElementwiseProductQuickTest<double> dense_vector_elementwise_product_quick_test_double("double");
