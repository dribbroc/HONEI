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

#include <libla/dense_matrix.hh>
#include <libla/dense_vector.hh>
#include <libla/matrix_vector_product.hh>
#include <libla/matrix_error.cc>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using namespace tests;

template <typename DataType_>
class DenseMatrixVectorProductTest :
    public BaseTest
{
    public:
        DenseMatrixVectorProductTest(const std::string & type) :
            BaseTest("dense_matrix_vector_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size, size + 1, DataType_(2));
                DenseVector<DataType_> dv1(size, DataType_(3)),  dv2(size, DataType_(6 * size));
                DenseVector<DataType_> prod(MatrixVectorProduct<DataType_>::value(dv1, dm1));

                TEST_CHECK_EQUAL(prod, dv2);
            }

            DenseMatrix<DataType_> dm01(4, 3, static_cast<DataType_>(1));
            DenseVector<DataType_> dv01(4, static_cast<DataType_>(1));

            TEST_CHECK_THROWS(MatrixVectorProduct<DataType_>::value(dv01, dm01), MatrixRowsDoNotMatch);
        }
};

DenseMatrixVectorProductTest<float> dense_matrix_vector_product_test_float("float");
DenseMatrixVectorProductTest<double> dense_matrix_vector_product_test_double("double");

template <typename DataType_>
class DenseMatrixVectorProductQuickTest :
    public QuickTest
{
    public:
        DenseMatrixVectorProductQuickTest(const std::string & type) :
            QuickTest("dense_matrix_vector_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseMatrix<DataType_> dm1(size + 1, size, DataType_(2));
            DenseVector<DataType_> dv1(size, DataType_(3)),  dv2(size, DataType_(6 * size));
            DenseVector<DataType_> prod(MatrixVectorProduct<DataType_>::value(dv1, dm1));

            TEST_CHECK_EQUAL(prod, dv2);

            DenseMatrix<DataType_> dm01(4, 3, static_cast<DataType_>(1));
            DenseVector<DataType_> dv01(4, static_cast<DataType_>(1));

            TEST_CHECK_THROWS(MatrixVectorProduct<DataType_>::value(dv01, dm01), MatrixRowsDoNotMatch);

        }
};
DenseMatrixVectorProductQuickTest<float> dense_matrix_vector_product_quick_test_float("float");
DenseMatrixVectorProductQuickTest<double> dense_matrix_vector_product_quick_test_double("double");
