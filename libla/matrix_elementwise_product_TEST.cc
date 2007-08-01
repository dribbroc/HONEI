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

#include <libla/matrix_elementwise_product.hh>
#include <libla/matrix_error.cc>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using namespace tests;

template <typename DataType_>
class BandedMatrixElementwiseProductTest :
    public BaseTest
{
    public:
        BandedMatrixElementwiseProductTest(const std::string & type) :
            BaseTest("banded_matrix_elementwise_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseVector<DataType_> * dv1 (new DenseVector<DataType_>(size, DataType_(3)));                    
                DenseVector<DataType_> * dv2 (new DenseVector<DataType_>(size, DataType_(2)));
                DenseVector<DataType_> * dv3 (new DenseVector<DataType_>(size, DataType_(6)));                  
                BandedMatrix<DataType_> bm1(size, dv1), bm2(size, dv2),  bm3(size, dv3); 
                BandedMatrix<DataType_> & prod(MatrixElementwiseProduct<DataType_>::value(bm1, bm2));

                TEST_CHECK_EQUAL(prod, bm3);
            }

            BandedMatrix<DataType_> bm01(5), bm02(6);

            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(bm02, bm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(bm02, bm01), MatrixColumnsDoNotMatch); 
        }
};
BandedMatrixElementwiseProductTest<float> banded_matrix_elementwise_product_test_float("float");
BandedMatrixElementwiseProductTest<double> banded_matrix_elementwise_product_test_double("double");

template <typename DataType_>
class BandedMatrixElementwiseProductQuickTest :
    public QuickTest
{
    public:
        BandedMatrixElementwiseProductQuickTest(const std::string & type) :
            QuickTest("banded_matrix_elementwise_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(20);
            DenseVector<DataType_> * dv1 (new DenseVector<DataType_>(size, DataType_(3)));                    
            DenseVector<DataType_> * dv2 (new DenseVector<DataType_>(size, DataType_(2)));
            DenseVector<DataType_> * dv3 (new DenseVector<DataType_>(size, DataType_(6)));                  
            BandedMatrix<DataType_> bm1(size, dv1), bm2(size, dv2),  bm3(size, dv3); 
            BandedMatrix<DataType_> & prod(MatrixElementwiseProduct<DataType_>::value(bm1, bm2));

            TEST_CHECK_EQUAL(prod, bm3);

            BandedMatrix<DataType_> bm01(5), bm02(6);

            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(bm02, bm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(bm02, bm01), MatrixColumnsDoNotMatch); 
        }
};
BandedMatrixElementwiseProductQuickTest<float> banded_matrix_elementwise_product_quick_test_float("float");
BandedMatrixElementwiseProductQuickTest<double> banded_matrix_elementwise_product_quick_test_double("double");

template <typename DataType_>
class DenseMatrixElementwiseProductTest :
    public BaseTest
{
    public:
        DenseMatrixElementwiseProductTest(const std::string & type) :
            BaseTest("dense_matrix_elementwise_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size, size + 1, DataType_(2)), dm2(size, size + 1, DataType_(3)),
                    dm3(size, size + 1, DataType_(6));
                DenseMatrix<DataType_> & prod(MatrixElementwiseProduct<DataType_>::value(dm1, dm2));

                TEST_CHECK_EQUAL(prod, dm3);
            }

            DenseMatrix<DataType_> dm01(2, 3, static_cast<DataType_>(1)), dm02(3, 4, static_cast<DataType_>(1)),
                dm03(2, 4, static_cast<DataType_>(1));

            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(dm03, dm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(dm03, dm02), MatrixColumnsDoNotMatch);
        }
};
DenseMatrixElementwiseProductTest<float> dense_matrix_elementwise_product_test_float("float");
DenseMatrixElementwiseProductTest<double> dense_matrix_elementwise_product_test_double("double");

template <typename DataType_>
class DenseMatrixElementwiseProductQuickTest :
    public QuickTest
{
    public:
        DenseMatrixElementwiseProductQuickTest(const std::string & type) :
            QuickTest("dense_matrix_elementwise_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseMatrix<DataType_> dm1(size, size + 1, DataType_(2)), dm2(size, size + 1, DataType_(3)),
                dm3(size, size + 1, DataType_(6));
            DenseMatrix<DataType_> & prod(MatrixElementwiseProduct<DataType_>::value(dm1, dm2));

            TEST_CHECK_EQUAL(prod, dm3);

            DenseMatrix<DataType_> dm01(2, 3, static_cast<DataType_>(1)), dm02(3, 4, static_cast<DataType_>(1)),
                dm03(2, 4, static_cast<DataType_>(1));

            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(dm03, dm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(dm03, dm02), MatrixColumnsDoNotMatch);
        }
};
DenseMatrixElementwiseProductQuickTest<float> dense_matrix_elementwise_product_quick_test_float("float");
DenseMatrixElementwiseProductQuickTest<double> dense_matrix_elementwise_product_quick_test_double("double");

template <typename DataType_>
class SparseMatrixElementwiseProductTest :
    public BaseTest
{
    public:
        SparseMatrixElementwiseProductTest(const std::string & type) :
            BaseTest("sparse_matrix_elementwise_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                SparseMatrix<DataType_> sm1(size, size + 1, size / 8 + 1), 
                    sm2(size, size + 1, size / 7 + 1), sm3(size, size + 1, size / 8 + 1 );
                for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()), 
                    i_end(sm1.end_elements()), j(sm2.begin_elements()), k(sm3.begin_elements()) ; 
                    i != i_end ; ++i, ++j, ++k)
                {
                    if (i.index() % 10 == 0) 
                    {
                        *i = DataType_(2);
                        *j = DataType_(3);
                        *k = DataType_(6);                        
                    }
                }  
                            
                SparseMatrix<DataType_> & prod(MatrixElementwiseProduct<DataType_>::value(sm1, sm2));

                TEST_CHECK_EQUAL(prod, sm3);
            }

            SparseMatrix<DataType_> sm01(5, 5, 1), sm02(6, 6, 1), sm03(5, 6, 1);

            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(sm03, sm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(sm03, sm02), MatrixColumnsDoNotMatch); 
        }
};
SparseMatrixElementwiseProductTest<float> sparse_matrix_elementwise_product_test_float("float");
SparseMatrixElementwiseProductTest<double> sparse_matrix_elementwise_product_test_double("double");

template <typename DataType_>
class SparseMatrixElementwiseProductQuickTest :
    public QuickTest
{
    public:
        SparseMatrixElementwiseProductQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_elementwise_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(22);
            SparseMatrix<DataType_> sm1(size, size + 1, size / 8 + 1), 
                sm2(size, size + 1, size / 7 + 1), sm3(size, size + 1, size / 8 + 1 );
            for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()), 
                i_end(sm1.end_elements()), j(sm2.begin_elements()), k(sm3.begin_elements()) ; 
                i != i_end ; ++i, ++j, ++k)
            {
                if (i.index() % 10 == 0) 
                {
                    *i = DataType_(2);
                    *j = DataType_(3);
                    *k = DataType_(6);                        
                }
            }  
                        
            SparseMatrix<DataType_> & prod(MatrixElementwiseProduct<DataType_>::value(sm1, sm2));

            TEST_CHECK_EQUAL(prod, sm3);

            SparseMatrix<DataType_> sm01(5, 5, 1), sm02(6, 6, 1), sm03(5, 6, 1);

            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(sm03, sm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(sm03, sm02), MatrixColumnsDoNotMatch); 
        }
};
SparseMatrixElementwiseProductQuickTest<float> sparse_matrix_elementwise_product_quick_test_float("float");
SparseMatrixElementwiseProductQuickTest<double> sparse_matrix_elementwise_product_quick_test_double("double");
