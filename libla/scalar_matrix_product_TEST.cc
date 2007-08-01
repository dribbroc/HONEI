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
#include <libla/scalar_matrix_product.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using namespace tests;

template <typename DataType_>
class ScalarBandedMatrixProductTest :
    public BaseTest
{
    public:
        ScalarBandedMatrixProductTest(const std::string & type) :
            BaseTest("scalar_banded_matrix_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseVector<DataType_> * dv1 (new DenseVector<DataType_>(size, static_cast<DataType_>(2)));                 
                DenseVector<DataType_> * dv2 (new DenseVector<DataType_>(size, static_cast<DataType_>(6)));                   
                BandedMatrix<DataType_> bm1(size, dv1), bm2(size, dv2);            
                BandedMatrix<DataType_> & prod(ScalarMatrixProduct<DataType_>::value(DataType_(3), bm1));

                TEST_CHECK_EQUAL(prod, bm2);
            }
        }
};
ScalarBandedMatrixProductTest<float> scalar_banded_matrix_product_test_float("float");
ScalarBandedMatrixProductTest<double> scalar_banded_matrix_product_test_double("double");

template <typename DataType_>
class ScalarBandedMatrixProductQuickTest :
    public QuickTest
{
    public:
        ScalarBandedMatrixProductQuickTest(const std::string & type) :
            QuickTest("scalar_banded_matrix_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(20);
            DenseVector<DataType_> * dv1 (new DenseVector<DataType_>(size, static_cast<DataType_>(2)));                 
            DenseVector<DataType_> * dv2 (new DenseVector<DataType_>(size, static_cast<DataType_>(6)));                   
            BandedMatrix<DataType_> bm1(size, dv1), bm2(size, dv2);            
            BandedMatrix<DataType_> & prod(ScalarMatrixProduct<DataType_>::value(DataType_(3), bm1));

            TEST_CHECK_EQUAL(prod, bm2);
        }
};
ScalarBandedMatrixProductQuickTest<float> scalar_banded_matrix_product_quick_test_float("float");
ScalarBandedMatrixProductQuickTest<double> scalar_banded_matrix_product_quick_test_double("double");

template <typename DataType_>
class ScalarDenseMatrixProductTest :
    public BaseTest
{
    public:
        ScalarDenseMatrixProductTest(const std::string & type) :
            BaseTest("scalar_dense_matrix_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size, size + 1, DataType_(2)), dm2(size, size + 1, DataType_(6));
                DenseMatrix<DataType_> & prod(ScalarMatrixProduct<DataType_>::value(DataType_(3), dm1));

                TEST_CHECK_EQUAL(prod, dm2);
            }
        }
};
ScalarDenseMatrixProductTest<float> scalar_dense_matrix_product_test_float("float");
ScalarDenseMatrixProductTest<double> scalar_dense_matrix_product_test_double("double");

template <typename DataType_>
class ScalarDenseMatrixProductQuickTest :
    public QuickTest
{
    public:
        ScalarDenseMatrixProductQuickTest(const std::string & type) :
            QuickTest("scalar_dense_matrix_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseMatrix<DataType_> dm(size, size + 1, DataType_(2));
            DenseMatrix<DataType_> prod1(ScalarMatrixProduct<DataType_>::value(DataType_(3), dm));

            DataType_ vsum(0), ssum(2 * DataType_(size) * DataType_(size + 1));
            for (typename MutableMatrix<DataType_>::ElementIterator i(dm.begin_elements()),
                    i_end(dm.end_elements()) ; i != i_end ; ++i)
            {
                vsum += *i;
            }
            TEST_CHECK_EQUAL_WITHIN_EPS(vsum, static_cast<DataType_>(6 * size * (size +1)), 
                std::numeric_limits<DataType_>::epsilon());
        }
};
ScalarDenseMatrixProductQuickTest<float> scalar_dense_matrix_product_quick_test_float("float");
ScalarDenseMatrixProductQuickTest<double> scalar_dense_matrix_product_quick_test_double("double");

template <typename DataType_>
class ScalarSparseMatrixProductTest :
    public BaseTest
{
    public:
        ScalarSparseMatrixProductTest(const std::string & type) :
            BaseTest("scalar_sparse_matrix_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                SparseMatrix<DataType_> sm1(size, size + 1, size / 8 + 1), 
                    sm2(size, size + 1, size / 7 + 1);
                for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()), 
                    i_end(sm1.end_elements()), j(sm2.begin_elements()) ; i != i_end ; ++i, ++j)
                {
                    if (i.index() % 10 == 0) 
                    {
                        *i = DataType_(2);
                        *j = DataType_(6);
                    }
                }             
                SparseMatrix<DataType_> & prod(ScalarMatrixProduct<DataType_>::value(DataType_(3), sm1));

                TEST_CHECK_EQUAL(prod, sm2);
            }
        }
};
ScalarSparseMatrixProductTest<float> scalar_sparse_matrix_product_test_float("float");
ScalarSparseMatrixProductTest<double> scalar_sparse_matrix_product_test_double("double");

template <typename DataType_>
class ScalarSparseMatrixProductQuickTest :
    public QuickTest
{
    public:
        ScalarSparseMatrixProductQuickTest(const std::string & type) :
            QuickTest("scalar_sparse_matrix_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size (22);
            SparseMatrix<DataType_> sm1(size, size + 1, size / 8 + 1), 
                sm2(size, size + 1, size / 7 + 1);
            for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()), 
                i_end(sm1.end_elements()), j(sm2.begin_elements()) ; i != i_end ; ++i, ++j)
            {
                if (i.index() % 10 == 0) 
                {
                    *i = DataType_(2);
                    *j = DataType_(6);
                }
            }             
            SparseMatrix<DataType_> & prod(ScalarMatrixProduct<DataType_>::value(DataType_(3), sm1));

            TEST_CHECK_EQUAL(prod, sm2);
        }
};
ScalarSparseMatrixProductQuickTest<float> scalar_sparse_matrix_product_quick_test_float("float");
ScalarSparseMatrixProductQuickTest<double> scalar_sparse_matrix_product_quick_test_double("double");
