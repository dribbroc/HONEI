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
#include <libla/matrix_sum.hh>
#include <matrix_error.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using namespace tests;

template <typename DataType_>
class BandedMatrixSumTest :
    public BaseTest
{
    public:
        BandedMatrixSumTest(const std::string & type) :
            BaseTest("banded_matrix_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseVector<DataType_> * dv1 (new DenseVector<DataType_>(size, static_cast<DataType_>(2)));                    
                DenseVector<DataType_> * dv2 (new DenseVector<DataType_>(size, static_cast<DataType_>(3))); 
                BandedMatrix<DataType_> bm1(size, dv1), bm2(size, dv2);
                BandedMatrix<DataType_> & sum(MatrixSum<DataType_>::value(bm1, bm2));
                for (typename BandedMatrix<DataType_>::ConstVectorIterator ce(sum.begin_bands()), ce_end(sum.end_bands()) ;
                    ce != ce_end ; ++ce)
                {
                    DenseVector<DataType_>  dv3 = *ce;
                    if (ce.index() != size - 1 )
                    {
                        for (typename Vector<DataType_>::ConstElementIterator i(dv3.begin_elements()), 
                            i_end(dv3.end_elements()) ; i != i_end ; ++i)                    
                        {
                            TEST_CHECK_EQUAL_WITHIN_EPS(*i, 0, std::numeric_limits<DataType_>::epsilon());
                        }
                    }
                    else
                    {
                        for (typename Vector<DataType_>::ConstElementIterator i(dv3.begin_elements()), 
                            i_end(dv3.end_elements()) ; i != i_end ; ++i)                    
                        {
                            TEST_CHECK_EQUAL_WITHIN_EPS(*i, 5, std::numeric_limits<DataType_>::epsilon());
                        }
                    } 
                }
            }

            BandedMatrix<DataType_> bm01(5), bm02(6);

            TEST_CHECK_THROWS(MatrixSum<DataType_>::value(bm02, bm01), MatrixSizeDoesNotMatch);
        }
};
BandedMatrixSumTest<float> banded_matrix_sum_test_float("float");
BandedMatrixSumTest<double> banded_matrix_sum_test_double("double");

template <typename DataType_>
class BandedMatrixSumQuickTest :
    public QuickTest
{
    public:
        BandedMatrixSumQuickTest(const std::string & type) :
            QuickTest("banded_matrix_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(22);
            DenseVector<DataType_> * dv1 (new DenseVector<DataType_>(size, DataType_(2)));                    
            DenseVector<DataType_> * dv2 (new DenseVector<DataType_>(size, DataType_(3))); 
            BandedMatrix<DataType_> bm1(size, dv1), bm2(size, dv2);
            BandedMatrix<DataType_> & sum(MatrixSum<DataType_>::value(bm1, bm2));
            for (typename BandedMatrix<DataType_>::ConstVectorIterator ce(sum.begin_bands()), ce_end(sum.end_bands()) ;
                ce != ce_end ; ++ce)
            {
                DenseVector<DataType_>  dv3 = *ce;
                if (ce.index() != size - 1 )
                {
                    for (typename Vector<DataType_>::ConstElementIterator i(dv3.begin_elements()), 
                        i_end(dv3.end_elements()) ; i != i_end ; ++i)                    
                    {
                        TEST_CHECK_EQUAL_WITHIN_EPS(*i, 0, std::numeric_limits<DataType_>::epsilon());
                    }
                }
                else
                {
                    for (typename Vector<DataType_>::ConstElementIterator i(dv3.begin_elements()), 
                        i_end(dv3.end_elements()) ; i != i_end ; ++i)                    
                    {
                        TEST_CHECK_EQUAL_WITHIN_EPS(*i, 5, std::numeric_limits<DataType_>::epsilon());
                    }
                } 
            }

            BandedMatrix<DataType_> bm01(5), bm02(6);

            TEST_CHECK_THROWS(MatrixSum<DataType_>::value(bm02, bm01), MatrixSizeDoesNotMatch);              
        }
};
BandedMatrixSumQuickTest<float> banded_matrix_sum_quick_test_float("float");
BandedMatrixSumQuickTest<double> banded_matrix_sum_quick_test_double("double");

template <typename DataType_>
class DenseMatrixSumTest :
    public BaseTest
{
    public:
        DenseMatrixSumTest(const std::string & type) :
            BaseTest("dense_matrix_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size, size + 1, DataType_(2)), dm2(size, size + 1, DataType_(-3)),
                    dm3(size, size + 1, DataType_(-1));
                DenseMatrix<DataType_> & sum(MatrixSum<DataType_>::value(dm1, dm2));

                TEST_CHECK_EQUAL(sum, dm3);
            }

            DenseMatrix<DataType_> dm01(5, 5), dm02(6, 6), dm03(5, 6);

            TEST_CHECK_THROWS(MatrixSum<DataType_>::value(dm03, dm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(MatrixSum<DataType_>::value(dm03, dm02), MatrixColumnsDoNotMatch);            
        }
};
DenseMatrixSumTest<float> dense_matrix_sum_test_float("float");
DenseMatrixSumTest<double> dense_matrix_sum_test_double("double");

template <typename DataType_>
class DenseMatrixSumQuickTest :
    public QuickTest
{
    public:
        DenseMatrixSumQuickTest(const std::string & type) :
            QuickTest("dense_matrix_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseMatrix<DataType_> dm1(size, size + 1, DataType_(2)), dm2(size, size + 1, DataType_(-3)),
                dm3(size, size + 1, DataType_(-1));
            DenseMatrix<DataType_> & sum(MatrixSum<DataType_>::value(dm1, dm2));

            TEST_CHECK_EQUAL(sum, dm3);

            DenseMatrix<DataType_> dm01(5, 5), dm02(6, 6), dm03(5, 6);

            TEST_CHECK_THROWS(MatrixSum<DataType_>::value(dm03, dm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(MatrixSum<DataType_>::value(dm03, dm02), MatrixColumnsDoNotMatch);     
        }
};
DenseMatrixSumQuickTest<float> dense_matrix_sum_quick_test_float("float");
DenseMatrixSumQuickTest<double> dense_matrix_sum_quick_test_double("double");

template <typename DataType_>
class SparseMatrixSumTest :
    public BaseTest
{
    public:
        SparseMatrixSumTest(const std::string & type) :
            BaseTest("sparse_matrix_sum_test<" + type + ">")
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
                        *j = DataType_(-3);
                        *k = DataType_(-1);                        
                    }
                }              
                SparseMatrix<DataType_> & sum(MatrixSum<DataType_>::value(sm1, sm2));

                TEST_CHECK_EQUAL(sum, sm3);
            }

            SparseMatrix<DataType_> sm01(5, 5, 1), sm02(6, 6, 1), sm03(5, 6, 1);

            TEST_CHECK_THROWS(MatrixSum<DataType_>::value(sm03, sm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(MatrixSum<DataType_>::value(sm03, sm02), MatrixColumnsDoNotMatch);           
        }
};
SparseMatrixSumTest<float> sparse_matrix_sum_test_float("float");
SparseMatrixSumTest<double> sparse_matrix_sum_test_double("double");

template <typename DataType_>
class SparseMatrixSumQuickTest :
    public QuickTest
{
    public:
        SparseMatrixSumQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size (22);
            SparseMatrix<DataType_> sm1(size, size + 1, size / 8 + 1), 
                sm2(size, size + 1, size / 7 + 1), sm3(size, size + 1, size / 8 + 1 );
            for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()), 
                i_end(sm1.end_elements()), j(sm2.begin_elements()), k(sm3.begin_elements()) ; 
                i != i_end ; ++i, ++j, ++k)
            {
                if (i.index() % 10 == 0) 
                {
                    *i = DataType_(2);
                    *j = DataType_(-3);
                    *k = DataType_(-1);                        
                }
            }              
            SparseMatrix<DataType_> & sum(MatrixSum<DataType_>::value(sm1, sm2));

            TEST_CHECK_EQUAL(sum, sm3);

            SparseMatrix<DataType_> sm01(5, 5, 1), sm02(6, 6, 1), sm03(5, 6, 1);

            TEST_CHECK_THROWS(MatrixSum<DataType_>::value(sm03, sm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(MatrixSum<DataType_>::value(sm03, sm02), MatrixColumnsDoNotMatch);           
        }
};
SparseMatrixSumQuickTest<float> sparse_matrix_sum_quick_test_float("float");
SparseMatrixSumQuickTest<double> sparse_matrix_sum_quick_test_double("double");
