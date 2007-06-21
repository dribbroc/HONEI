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
#include <libla/matrix_element_sum.hh>
#include <matrix_error.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using namespace tests;

template <typename DataType_>
class DenseMatrixElementSumTest :
    public BaseTest
{
    public:
        DenseMatrixElementSumTest(const std::string & type) :
            BaseTest("dense_matrix_element_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size, size + 1, DataType_(1));
                DataType_ sum(MatrixElementSum<DataType_>::value(dm1));

				TEST_CHECK_EQUAL(sum, size * (size + 1));
            }
        }
};

DenseMatrixElementSumTest<float> dense_matrix_element_sum_test_float("float");
DenseMatrixElementSumTest<double> dense_matrix_element_sum_test_double("double");

template <typename DataType_>
class DenseMatrixElementSumQuickTest :
    public QuickTest
{
    public:
        DenseMatrixElementSumQuickTest(const std::string & type) :
            QuickTest("dense_matrix_element_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseMatrix<DataType_> dm1(size, size + 1, DataType_(1));
            DataType_ sum(MatrixElementSum<DataType_>::value(dm1));

			TEST_CHECK_EQUAL(sum, size * (size + 1));
        }
};
DenseMatrixElementSumQuickTest<float> dense_matrix_element_sum_quick_test_float("float");
DenseMatrixElementSumQuickTest<double> dense_matrix_element_sum_quick_test_double("double");
