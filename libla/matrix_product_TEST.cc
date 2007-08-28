/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2007 Danny van Dyk <danny.dyk@uni-dortmund.de>
 * Copyright (c) 2007 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 * Copyright (c) 2007 Sven Mallach <sven.mallach@uni-dortmund.de>
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
class BandedMatrixProductTest :
    public BaseTest
{
    public:
        BandedMatrixProductTest(const std::string & type) :
            BaseTest("banded_matrix_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseVector<DataType_> * dv1 (new DenseVector<DataType_>(size, DataType_(2)));                    
                DenseVector<DataType_> * dv2 (new DenseVector<DataType_>(size, DataType_(3))); 
                DenseVector<DataType_> * dv3 (new DenseVector<DataType_>(size, DataType_(6)));                 
                BandedMatrix<DataType_> bm1(size, dv1), bm2(size, dv2), bm3(size, dv3);            
                BandedMatrix<DataType_> prod(MatrixProduct<DataType_>::value(bm1, bm2));

                TEST_CHECK_EQUAL(prod, bm3);
            }
            BandedMatrix<DataType_> bm01(5), bm02(6);

            TEST_CHECK_THROWS(MatrixProduct<DataType_>::value(bm02, bm01), MatrixRowsDoNotMatch);
        }
};
BandedMatrixProductTest<float> banded_matrix_product_test_float("float");
BandedMatrixProductTest<double> banded_matrix_product_test_double("double");

template <typename DataType_>
class BandedMatrixProductQuickTest :
    public QuickTest
{
    public:
        BandedMatrixProductQuickTest(const std::string & type) :
            QuickTest("banded_matrix_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseVector<DataType_> * dv1 (new DenseVector<DataType_>(size, DataType_(2)));                    
            DenseVector<DataType_> * dv2 (new DenseVector<DataType_>(size, DataType_(3))); 
            DenseVector<DataType_> * dv3 (new DenseVector<DataType_>(size, DataType_(6)));                 
            BandedMatrix<DataType_> bm1(size, dv1), bm2(size, dv2), bm3(size, dv3);            
            BandedMatrix<DataType_> prod(MatrixProduct<DataType_>::value(bm1, bm2));

            TEST_CHECK_EQUAL(prod, bm3);

            BandedMatrix<DataType_> bm01(5), bm02(6);

            TEST_CHECK_THROWS(MatrixProduct<DataType_>::value(bm02, bm01), MatrixRowsDoNotMatch);
        }
};
BandedMatrixProductQuickTest<float> banded_matrix_product_quick_test_float("float");
BandedMatrixProductQuickTest<double> banded_matrix_product_quick_test_double("double");

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
                    dm3(size+1, size+1, DataType_(6 * size));
                DenseMatrix<DataType_> prod(MatrixProduct<DataType_>::value(dm1, dm2));

                TEST_CHECK_EQUAL(prod, dm3);
            }

            DenseMatrix<DataType_> dm01(3, 3, DataType_(1)), dm02(4, 3, DataType_(1));

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
                dm3(size+1, size+1, DataType_(6 * size));
            DenseMatrix<DataType_> prod(MatrixProduct<DataType_>::value(dm1, dm2));

            TEST_CHECK_EQUAL(prod, dm3);

            DenseMatrix<DataType_> dm01(3, 3, DataType_(1)), dm02(4, 3, DataType_(1));

            TEST_CHECK_THROWS(MatrixProduct<DataType_>::value(dm02, dm01), MatrixRowsDoNotMatch);

        }
};
DenseMatrixProductQuickTest<float> dense_matrix_product_quick_test_float("float");
DenseMatrixProductQuickTest<double> dense_matrix_product_quick_test_double("double");

template <typename DataType_>
class DenseMatrixSparseMatrixProductTest :
    public BaseTest
{
    public:
        DenseMatrixSparseMatrixProductTest(const std::string & type) :
            BaseTest("dense_matrix_sparse_matrix_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size, size + 1, DataType_(2)), dm3(size+1, size+1, DataType_(6 * size));
                SparseMatrix<DataType_> sm2(size, size + 1, size / 8 + 1);
                for (typename MutableMatrix<DataType_>::ElementIterator i(sm2.begin_elements()), 
                    i_end(sm2.end_elements()) ; i != i_end ; ++i)
                {
                    *i = 3;
                }                    
                DenseMatrix<DataType_> prod(MatrixProduct<DataType_>::value(dm1, sm2));

                TEST_CHECK_EQUAL(prod, dm3);
            }

            SparseMatrix<DataType_> sm01(3, 3, 1);
            DenseMatrix<DataType_> dm02(4, 3);

            TEST_CHECK_THROWS(MatrixProduct<DataType_>::value(dm02, sm01), MatrixRowsDoNotMatch);
        }
};
DenseMatrixSparseMatrixProductTest<float> dense_matrix_sparse_matrix_product_test_float("float");
DenseMatrixSparseMatrixProductTest<double> dense_matrix_sparse_matrix_product_test_double("double");

template <typename DataType_>
class DenseMatrixSparseMatrixProductQuickTest :
    public QuickTest
{
    public:
        DenseMatrixSparseMatrixProductQuickTest(const std::string & type) :
            QuickTest("dense_matrix_sparse_matrix_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size (11);
            DenseMatrix<DataType_> dm1(size, size + 1, DataType_(2)), dm3(size+1, size+1, DataType_(6 * size));
            SparseMatrix<DataType_> sm2(size + 1, size, size / 8 + 1);
            for (typename MutableMatrix<DataType_>::ElementIterator i(sm2.begin_elements()), 
                i_end(sm2.end_elements()) ; i != i_end ; ++i)
            {
                *i = 3;
            }                    
            DenseMatrix<DataType_> prod(MatrixProduct<DataType_>::value(dm1, sm2));

            TEST_CHECK_EQUAL(prod, dm3);

            SparseMatrix<DataType_> sm01(3, 3, 1);
            DenseMatrix<DataType_> dm02(4, 3);

            TEST_CHECK_THROWS(MatrixProduct<DataType_>::value(dm02, sm01), MatrixRowsDoNotMatch);
        }
};
DenseMatrixSparseMatrixProductQuickTest<float> dense_matrix_sparse_matrix_product_quick_test_float("float");
DenseMatrixSparseMatrixProductQuickTest<double> dense_matrix_sparse_matrix_product_quick_test_double("double");

template <typename DataType_>
class SparseMatrixDenseMatrixProductTest :
    public BaseTest
{
    public:
        SparseMatrixDenseMatrixProductTest(const std::string & type) :
            BaseTest("sparse_matrix_dense_matrix_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm2(size + 1, size, DataType_(2)), dm3(size+1, size+1, DataType_(6 * size));
                SparseMatrix<DataType_> sm1(size, size + 1, size / 8 + 1);
                for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()), 
                    i_end(sm1.end_elements()) ; i != i_end ; ++i)
                {
                    *i = 3;
                }                    
                DenseMatrix<DataType_> prod(MatrixProduct<DataType_>::value(sm1, dm2));

                TEST_CHECK_EQUAL(prod, dm3);
            }

            SparseMatrix<DataType_> sm01(4, 3, 1);
            DenseMatrix<DataType_> dm02(3, 3);

            TEST_CHECK_THROWS(MatrixProduct<DataType_>::value(sm01, dm02), MatrixRowsDoNotMatch);
        }
};
SparseMatrixDenseMatrixProductTest<float> sparse_matrix_dense_matrix_product_test_float("float");
SparseMatrixDenseMatrixProductTest<double> sparse_matrix_dense_matrix_product_test_double("double");

template <typename DataType_>
class SparseMatrixDenseMatrixProductQuickTest :
    public QuickTest
{
    public:
        SparseMatrixDenseMatrixProductQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_dense_matrix_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(11);
            DenseMatrix<DataType_> dm2(size + 1, size, DataType_(2)), dm3(size+1, size+1, DataType_(6 * size));
            SparseMatrix<DataType_> sm1(size, size + 1, size / 8 + 1);
            for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()), 
                i_end(sm1.end_elements()) ; i != i_end ; ++i)
            {
                *i = 3;
            }                    
            DenseMatrix<DataType_> prod(MatrixProduct<DataType_>::value(sm1, dm2));

            TEST_CHECK_EQUAL(prod, dm3);

            SparseMatrix<DataType_> sm01(4, 3, 1);
            DenseMatrix<DataType_> dm02(3, 3);

            TEST_CHECK_THROWS(MatrixProduct<DataType_>::value(sm01, dm02), MatrixRowsDoNotMatch);
        }
};
SparseMatrixDenseMatrixProductQuickTest<float> sparse_matrix_dense_matrix_product_quick_test_float("float");
SparseMatrixDenseMatrixProductQuickTest<double> sparse_matrix_dense_matrix_product_quick_test_double("double");

template <typename DataType_>
class SparseMatrixProductTest :
    public BaseTest
{
    public:
        SparseMatrixProductTest(const std::string & type) :
            BaseTest("sparse_matrix_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                SparseMatrix<DataType_> sm1(size, size + 1, size / 8 + 1), 
                    sm2(size + 1, size, size / 7 + 1), sm3(size + 1, size + 1, size / 7 + 1);
                for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()), 
                    i_end(sm1.end_elements()) ; i != i_end ; ++i)
                {
                    *i = 2;
                }
                for (typename MutableMatrix<DataType_>::ElementIterator i(sm2.begin_elements()), 
                    i_end(sm2.end_elements()) ; i != i_end ; ++i)
                {
                    *i = 3;
                } 
                for (typename MutableMatrix<DataType_>::ElementIterator i(sm3.begin_elements()), 
                    i_end(sm3.end_elements()) ; i != i_end ; ++i)
                {
                    *i = 6 * size;
                }                                          

                SparseMatrix<DataType_> prod(MatrixProduct<DataType_>::value(sm1, sm2));

                TEST_CHECK_EQUAL(prod, sm3);
            }
            
            SparseMatrix<DataType_> sm01(3, 3, 1), sm02(4, 3, 1);

            TEST_CHECK_THROWS(MatrixProduct<DataType_>::value(sm02, sm01), MatrixRowsDoNotMatch);
        }
};
SparseMatrixProductTest<float> sparse_matrix_product_test_float("float");
SparseMatrixProductTest<double> sparse_matrix_product_test_double("double");

template <typename DataType_>
class SparseMatrixProductQuickTest :
    public QuickTest
{
    public:
        SparseMatrixProductQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(2);
            SparseMatrix<DataType_> sm1(size, size + 1, size / 8 + 1), 
                sm2(size + 1, size, size / 7 + 1), sm3(size + 1, size + 1, size / 7 + 1);
            for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()), 
                i_end(sm1.end_elements()) ; i != i_end ; ++i)
            {
                *i = 2;
            }
            for (typename MutableMatrix<DataType_>::ElementIterator i(sm2.begin_elements()), 
                i_end(sm2.end_elements()) ; i != i_end ; ++i)
            {
                *i = 3;
            } 
            for (typename MutableMatrix<DataType_>::ElementIterator i(sm3.begin_elements()), 
                i_end(sm3.end_elements()) ; i != i_end ; ++i)
            {
                *i = 6 * size;
            }                                          

            SparseMatrix<DataType_> prod(MatrixProduct<DataType_>::value(sm1, sm2));

            TEST_CHECK_EQUAL(prod, sm3);

            SparseMatrix<DataType_> sm01(3, 3, 1), sm02(4, 3, 1);

            TEST_CHECK_THROWS(MatrixProduct<DataType_>::value(sm02, sm01), MatrixRowsDoNotMatch);
        }
};
SparseMatrixProductQuickTest<float> sparse_matrix_product_quick_test_float("float");
SparseMatrixProductQuickTest<double> sparse_matrix_product_quick_test_double("double");
