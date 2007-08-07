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
class BandedMatrixElementSumTest :
    public BaseTest
{
    public:
        BandedMatrixElementSumTest(const std::string & type) :
            BaseTest("banded_matrix_element_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseVector<DataType_> * dv1 (new DenseVector<DataType_>(size, DataType_(2)));                    
                BandedMatrix<DataType_> bm1(size, dv1);
                DataType_ sum(MatrixElementSum<>::value(bm1));            

				TEST_CHECK_EQUAL(sum, size * 2);
            }
        }
};
BandedMatrixElementSumTest<float> banded_matrix_element_sum_test_float("float");
BandedMatrixElementSumTest<double> banded_matrix_element_sum_test_double("double");

template <typename DataType_>
class BandedMatrixElementSumQuickTest :
    public QuickTest
{
    public:
        BandedMatrixElementSumQuickTest(const std::string & type) :
            QuickTest("banded_matrix_element_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(20);
            DenseVector<DataType_> * dv1 (new DenseVector<DataType_>(size, DataType_(2)));                    
            BandedMatrix<DataType_> bm1(size, dv1);
            DataType_ sum(MatrixElementSum<>::value(bm1));            

			TEST_CHECK_EQUAL(sum, size * 2);
        }
};
BandedMatrixElementSumQuickTest<float> banded_matrix_element_sum_quick_test_float("float");
BandedMatrixElementSumQuickTest<double> banded_matrix_element_sum_quick_test_double("double");

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
                DataType_ sum(MatrixElementSum<>::value(dm1));

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
            DataType_ sum(MatrixElementSum<>::value(dm1));

			TEST_CHECK_EQUAL(sum, size * (size + 1));
        }
};
DenseMatrixElementSumQuickTest<float> dense_matrix_element_sum_quick_test_float("float");
DenseMatrixElementSumQuickTest<double> dense_matrix_element_sum_quick_test_double("double");

template <typename DataType_>
class SparseMatrixElementSumTest :
    public BaseTest
{
    public:
        SparseMatrixElementSumTest(const std::string & type) :
            BaseTest("dense_matrix_element_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(11) ; size < (1 << 12) ; size <<= 1)
            {
                SparseMatrix<DataType_> sm1(size, size + 1, size + 8 / 1);
                for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()), 
                    i_end(sm1.end_elements()) ; i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0) 
                    {
                        *i = DataType_(2);
                    }
                }                 
                DataType_ sum(MatrixElementSum<>::value(sm1));

				TEST_CHECK_EQUAL(sum, (size / 10 + 1) * 2);
            }
        }
};
SparseMatrixElementSumTest<float> sparse_matrix_element_sum_test_float("float");
SparseMatrixElementSumTest<double> sparse_matrix_element_sum_test_double("double");

template <typename DataType_>
class SparseMatrixElementSumQuickTest :
    public QuickTest
{
    public:
        SparseMatrixElementSumQuickTest(const std::string & type) :
            QuickTest("dense_matrix_element_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(22);
            SparseMatrix<DataType_> sm1(size, size + 1, size / 8 + 1);
            for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()), 
                i_end(sm1.end_elements()) ; i != i_end ; ++i)
            {
                if (i.index() % 10 == 0) 
                {
                    *i = DataType_(2);
                }
            }                 
            DataType_ sum(MatrixElementSum<>::value(sm1));

			TEST_CHECK_EQUAL(sum, ((size * size + 1) / 10 + 1) * 2);
        }
};
SparseMatrixElementSumQuickTest<float> sparse_matrix_element_sum_quick_test_float("float");
SparseMatrixElementSumQuickTest<double> sparse_matrix_element_sum_quick_test_double("double");
