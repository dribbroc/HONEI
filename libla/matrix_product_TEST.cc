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

#include <libla/dense_matrix.hh>
#include <libla/matrix_product.hh>
#include <libla/matrix_error.cc>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using namespace tests;

template <typename DataType_>
class DenseMatrixProductTest :
    public BaseTest
{
    public:
        DenseMatrixProductTest(const std::string & type) :
            BaseTest("dense_matrix_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size, size + 1, DataType_(2)), dm2(size + 1, size, DataType_(3)),
                    dm3(size, size, DataType_(6 * size));
                DenseMatrix<DataType_> & prod(MatrixProduct<DataType_>::value(dm1, dm2));

                for (typename Matrix<DataType_>::ConstElementIterator i(prod.begin_elements()),
                        i_end(prod.end_elements()), j(dm3.begin_elements()) ; i != i_end ; ++i, ++j)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(*i, *j, std::numeric_limits<DataType_>::epsilon());
                }
            }

            DenseMatrix<DataType_> dm01(3, 3, static_cast<DataType_>(1)), dm02(4, 3, static_cast<DataType_>(1));

            TEST_CHECK_THROWS(MatrixProduct<DataType_>::value(dm02, dm01), MatrixRowsDoNotMatch);
        }
};

DenseMatrixProductTest<float> dense_matrix_product_test_float("float");
DenseMatrixProductTest<double> dense_matrix_product_test_double("double");

template <typename DataType_>
class DenseMatrixProductQuickTest :
    public QuickTest
{
    public:
        DenseMatrixProductQuickTest(const std::string & type) :
            QuickTest("dense_matrix_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseMatrix<DataType_> dm1(size, size + 1, DataType_(2)), dm2(size + 1, size, DataType_(3)),
                dm3(size, size, DataType_(6 * size));
            DenseMatrix<DataType_> & prod(MatrixProduct<DataType_>::value(dm1, dm2));

            for (typename Matrix<DataType_>::ConstElementIterator i(prod.begin_elements()),
                    i_end(prod.end_elements()), j(dm3.begin_elements()) ; i != i_end ; ++i, ++j)
            {
                TEST_CHECK_EQUAL_WITHIN_EPS(*i, *j, std::numeric_limits<DataType_>::epsilon());
            }

            DenseMatrix<DataType_> dm01(3, 3, static_cast<DataType_>(1)), dm02(4, 3, static_cast<DataType_>(1));

            TEST_CHECK_THROWS(MatrixProduct<DataType_>::value(dm02, dm01), MatrixRowsDoNotMatch);
        }
};
DenseMatrixProductQuickTest<float> dense_matrix_product_quick_test_float("float");
DenseMatrixProductQuickTest<double> dense_matrix_product_quick_test_double("double");
