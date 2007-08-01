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
#include <libla/scalar_matrix_sum.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using namespace tests;

template <typename DataType_>
class ScalarBandedMatrixSumTest :
    public BaseTest
{
    public:
        ScalarBandedMatrixSumTest(const std::string & type) :
            BaseTest("scalar_banded_matrix_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseVector<DataType_> * dv1 (new DenseVector<DataType_>(size, static_cast<DataType_>(2)));                 
                DenseVector<DataType_> * dv2 (new DenseVector<DataType_>(size, static_cast<DataType_>(5)));                   
                BandedMatrix<DataType_> bm1(size, dv1), bm2(size, dv2);     
                BandedMatrix<DataType_> & prod(ScalarMatrixSum<DataType_>::value(DataType_(3), bm1));

                TEST_CHECK_EQUAL(prod, bm2);
            }
        }
};
ScalarBandedMatrixSumTest<float> scalar_banded_matrix_sum_test_float("float");
ScalarBandedMatrixSumTest<double> scalar_banded_matrix_sum_test_double("double");

template <typename DataType_>
class ScalarBandedMatrixSumQuickTest :
    public QuickTest
{
    public:
        ScalarBandedMatrixSumQuickTest(const std::string & type) :
            QuickTest("scalar_banded_matrix_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseVector<DataType_> * dv1 (new DenseVector<DataType_>(size, static_cast<DataType_>(2)));                 
                DenseVector<DataType_> * dv2 (new DenseVector<DataType_>(size, static_cast<DataType_>(5)));                   
                BandedMatrix<DataType_> bm1(size, dv1), bm2(size, dv2);     
                BandedMatrix<DataType_> & prod(ScalarMatrixSum<DataType_>::value(DataType_(3), bm1));

                TEST_CHECK_EQUAL(prod, bm2);
            }
        }
};
ScalarBandedMatrixSumQuickTest<float> scalar_banded_matrix_sum_quick_test_float("float");
ScalarBandedMatrixSumQuickTest<double> scalar_banded_matrix_sum_quick_test_double("double");

template <typename DataType_>
class ScalarDenseMatrixSumTest :
    public BaseTest
{
    public:
        ScalarDenseMatrixSumTest(const std::string & type) :
            BaseTest("scalar_dense_matrix_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size, size + 1, DataType_(2)), dm2(size, size + 1, DataType_(5));
                DenseMatrix<DataType_> & prod(ScalarMatrixSum<DataType_>::value(DataType_(3), dm1));

                TEST_CHECK_EQUAL(prod, dm2);
            }
        }
};
ScalarDenseMatrixSumTest<float> scalar_dense_matrix_sum_test_float("float");
ScalarDenseMatrixSumTest<double> scalar_dense_matrix_sum_test_double("double");

template <typename DataType_>
class ScalarDenseMatrixSumQuickTest :
    public QuickTest
{
    public:
        ScalarDenseMatrixSumQuickTest(const std::string & type) :
            QuickTest("scalar_dense_matrix_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseMatrix<DataType_> dm(size, size + 1, DataType_(2));
            DenseMatrix<DataType_> prod1(ScalarMatrixSum<DataType_>::value(DataType_(3), dm));

            DataType_ vsum(0), ssum(2 * DataType_(size) * DataType_(size + 1));
            for (typename MutableMatrix<DataType_>::ElementIterator i(dm.begin_elements()),
                    i_end(dm.end_elements()) ; i != i_end ; ++i)
            {
                vsum += *i;
            }
            TEST_CHECK_EQUAL_WITHIN_EPS(vsum, static_cast<DataType_>(5 * size * (size +1)), 
                std::numeric_limits<DataType_>::epsilon());
        }
};
ScalarDenseMatrixSumQuickTest<float> scalar_dense_matrix_sum_quick_test_float("float");
ScalarDenseMatrixSumQuickTest<double> scalar_dense_matrix_sum_quick_test_double("double");
/*
template <typename DataType_>
class ScalarSparseMatrixSumTest :
    public BaseTest
{
    public:
        ScalarSparseMatrixSumTest(const std::string & type) :
            BaseTest("scalar_sparse_matrix_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                SparseMatrix<DataType_> sm1(size, size + 1, size / 8 + 1);
                DenseMatrix<DataType_> dm1(size, size + 1, DataType_(3));                    
                for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()), 
                    i_end(sm1.end_elements()), j(dm1.begin_elements()) ; i != i_end ; ++i, ++j)
                {
                    if (i.index() % 10 == 0) 
                    {
                        *i = DataType_(2);
                        *j += DataType_(2);
                    }
                }              
                DenseMatrix<DataType_> & prod(ScalarMatrixSum<DataType_>::value(DataType_(3), sm1));

                TEST_CHECK_EQUAL(prod, dm1);
            }
        }
};
ScalarSparseMatrixSumTest<float> scalar_sparse_matrix_sum_test_float("float");
ScalarSparseMatrixSumTest<double> scalar_sparse_matrix_sum_test_double("double");

template <typename DataType_>
class ScalarSparseMatrixQuickSumTest :
    public QuickTest
{
    public:
        ScalarSparseMatrixQuickSumTest(const std::string & type) :
            QuickTest("scalar_sparse_matrix_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                SparseMatrix<DataType_> sm1(size, size + 1, size / 8 + 1);
                DenseMatrix<DataType_> dm1(size, size + 1, DataType_(3));                    
                for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()), 
                    i_end(sm1.end_elements()), j(dm1.begin_elements()) ; i != i_end ; ++i, ++j)
                {
                    if (i.index() % 10 == 0) 
                    {
                        *i = DataType_(2);
                        *j += DataType_(2);
                    }
                }              
                DenseMatrix<DataType_> & prod(ScalarMatrixSum<DataType_>::value(DataType_(3), sm1));

                TEST_CHECK_EQUAL(prod, dm1);
            }
        }
};
ScalarSparseMatrixQuickSumTest<float> scalar_sparse_matrix_sum_quick_test_float("float");
ScalarSparseMatrixQuickSumTest<double> scalar_sparse_matrix_sum_quick_test_double("double");*/
