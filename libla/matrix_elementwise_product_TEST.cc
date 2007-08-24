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
class BandedMatrixDenseMatrixElementwiseProductTest :
    public BaseTest
{
    public:
        BandedMatrixDenseMatrixElementwiseProductTest(const std::string & type) :
            BaseTest("banded_matrix_dense_matrix_elementwise_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size, size, DataType_(3));                  
                DenseVector<DataType_> * dv2 (new DenseVector<DataType_>(size, DataType_(2)));
                DenseVector<DataType_> * dv3 (new DenseVector<DataType_>(size, DataType_(6)));                  
                BandedMatrix<DataType_> bm2(size, dv2),  bm3(size, dv3); 
                
                BandedMatrix<DataType_> & prod(MatrixElementwiseProduct<DataType_>::value(bm2, dm1));

                TEST_CHECK_EQUAL(prod, bm3);
            }

            DenseMatrix<DataType_> dm01(5, 6), dm03(6, 5);
            BandedMatrix<DataType_> bm02(6);

            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(bm02, dm03), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(bm02, dm01), MatrixColumnsDoNotMatch);            
        }
};
BandedMatrixDenseMatrixElementwiseProductTest<float> banded_matrix_dense_matrix_elementwise_product_test_float("float");
BandedMatrixDenseMatrixElementwiseProductTest<double> banded_matrix_dense_matrix_elementwise_product_test_double("double");

template <typename DataType_>
class BandedMatrixDenseMatrixElementwiseProductQuickTest :
    public QuickTest
{
    public:
        BandedMatrixDenseMatrixElementwiseProductQuickTest(const std::string & type) :
            QuickTest("banded_matrix_dense_matrix_elementwise_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size (11);
            DenseMatrix<DataType_> dm1(size, size, DataType_(3));                  
            DenseVector<DataType_> * dv2 (new DenseVector<DataType_>(size, DataType_(2)));
            DenseVector<DataType_> * dv3 (new DenseVector<DataType_>(size, DataType_(6)));                  
            BandedMatrix<DataType_> bm2(size, dv2),  bm3(size, dv3); 
            
            BandedMatrix<DataType_> & prod(MatrixElementwiseProduct<DataType_>::value(bm2, dm1));

            TEST_CHECK_EQUAL(prod, bm3);

            DenseMatrix<DataType_> dm01(5, 6), dm03(6, 5);
            BandedMatrix<DataType_> bm02(6);

            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(bm02, dm03), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(bm02, dm01), MatrixColumnsDoNotMatch);            
        }
};
BandedMatrixDenseMatrixElementwiseProductQuickTest<float> banded_matrix_dense_matrix_elementwise_product_quick_test_float("float");
BandedMatrixDenseMatrixElementwiseProductQuickTest<double> banded_matrix_dense_matrix_elementwise_product_quick_test_double("double");

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

            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(bm02, bm01), MatrixSizeDoesNotMatch);
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

            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(bm02, bm01), MatrixSizeDoesNotMatch); 
        }
};
BandedMatrixElementwiseProductQuickTest<float> banded_matrix_elementwise_product_quick_test_float("float");
BandedMatrixElementwiseProductQuickTest<double> banded_matrix_elementwise_product_quick_test_double("double");

template <typename DataType_>
class BandedMatrixSparseMatrixElementwiseProductTest :
    public BaseTest
{
    public:
        BandedMatrixSparseMatrixElementwiseProductTest(const std::string & type) :
            BaseTest("banded_matrix_sparse_matrix_elementwise_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseVector<DataType_> * dv1 (new DenseVector<DataType_>(size, DataType_(3)));                    
                BandedMatrix<DataType_> bm1(size, dv1),  bm3(size); 
                SparseMatrix<DataType_> sm2(size, size, size / 7 + 1);
                typename BandedMatrix<DataType_>::ConstElementIterator i(bm1.begin_elements());            
                typename BandedMatrix<DataType_>::ElementIterator k(bm3.begin_elements());
                for (typename MutableMatrix<DataType_>::ElementIterator j(sm2.begin_elements()),
                    j_end(sm2.end_elements()) ; j != j_end ; ++i, ++j, ++k)
                {
                    if (j.index() % 7 == 0) 
                        *j = DataType_(2);

                    if (*i != DataType_(0))
                        *k = DataType_(6);                     
                } 

                BandedMatrix<DataType_> & prod(MatrixElementwiseProduct<DataType_>::value(bm1, sm2));

                TEST_CHECK_EQUAL(prod, bm3);
            }

            BandedMatrix<DataType_> bm01(5);
            SparseMatrix<DataType_> sm02(5, 6, 1), sm03(6, 5, 1);

            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(bm01, sm02), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(bm01, sm03), MatrixColumnsDoNotMatch);            
        }
};
BandedMatrixSparseMatrixElementwiseProductTest<float> banded_matrix_sparse_matrix_elementwise_product_test_float("float");
BandedMatrixSparseMatrixElementwiseProductTest<double> banded_matrix_sparse_matrix_elementwise_product_test_double("double");

template <typename DataType_>
class BandedMatrixSparseMatrixElementwiseProductQuickTest :
    public QuickTest
{
    public:
        BandedMatrixSparseMatrixElementwiseProductQuickTest(const std::string & type) :
            QuickTest("banded_matrix_sparse_matrix_elementwise_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(10);
            DenseVector<DataType_> * dv1 (new DenseVector<DataType_>(size, DataType_(3)));                    
            BandedMatrix<DataType_> bm1(size, dv1),  bm3(size); 
            SparseMatrix<DataType_> sm2(size, size, size / 7 + 1);
            typename BandedMatrix<DataType_>::ConstElementIterator i(bm1.begin_elements());            
            typename BandedMatrix<DataType_>::ElementIterator k(bm3.begin_elements());
            for (typename MutableMatrix<DataType_>::ElementIterator j(sm2.begin_elements()),
                j_end(sm2.end_elements()) ; j != j_end ; ++i, ++j, ++k)
            {
                if (j.index() % 7 == 0) 
                    *j = DataType_(2);

                if (*i != DataType_(0))
                    *k = DataType_(6);                     
            } 

            BandedMatrix<DataType_> & prod(MatrixElementwiseProduct<DataType_>::value(bm1, sm2));

            TEST_CHECK_EQUAL(prod, bm3);

            BandedMatrix<DataType_> bm01(5);
            SparseMatrix<DataType_> sm02(5, 6, 1), sm03(6, 5, 1);

            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(bm01, sm02), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(bm01, sm03), MatrixColumnsDoNotMatch);            
        }
};
BandedMatrixSparseMatrixElementwiseProductQuickTest<float> banded_matrix_sparse_matrix_elementwise_product_quick_test_float("float");
BandedMatrixSparseMatrixElementwiseProductQuickTest<double> banded_matrix_sparse_matrix_elementwise_product_quick_test_double("double");

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

            DenseMatrix<DataType_> dm01(2, 3, static_cast<DataType_>(1)), dm02(3, 4, DataType_(1)),
                dm03(2, 4, DataType_(1));

            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(dm03, dm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(dm03, dm02), MatrixColumnsDoNotMatch);
        }
};
DenseMatrixElementwiseProductQuickTest<float> dense_matrix_elementwise_product_quick_test_float("float");
DenseMatrixElementwiseProductQuickTest<double> dense_matrix_elementwise_product_quick_test_double("double");

template <typename DataType_>
class DenseMatrixSparseMatrixElementwiseProductTest :
    public BaseTest
{
    public:
        DenseMatrixSparseMatrixElementwiseProductTest(const std::string & type) :
            BaseTest("dense_matrix_sparse_matrix_elementwise_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size, size + 1, DataType_(2)), dm3(size, size + 1, DataType_(0));
                SparseMatrix<DataType_> sm2(size, size + 1, size / 7 + 1);
                for (typename MutableMatrix<DataType_>::ElementIterator j_end(sm2.end_elements()), j(sm2.begin_elements()), 
                    k(dm3.begin_elements()) ; j != j_end ; ++j, ++k)
                {
                    if (j.index() % 7 == 0) 
                    {
                        *j = DataType_(3);
                        *k = DataType_(6);                        
                    }
                }  
                            
                DenseMatrix<DataType_> & prod(MatrixElementwiseProduct<DataType_>::value(dm1, sm2));

                TEST_CHECK_EQUAL(prod, dm3);
            }

            DenseMatrix<DataType_> dm03(5, 6);
            SparseMatrix<DataType_> sm01(5, 5, 1), sm02(6, 6, 1), sm03(5, 6, 1);

            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(dm03, sm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(dm03, sm02), MatrixColumnsDoNotMatch); 
        }
};
DenseMatrixSparseMatrixElementwiseProductTest<float> dense_matrix_sparse_matrix_elementwise_product_test_float("float");
DenseMatrixSparseMatrixElementwiseProductTest<double> dense_matrix_sparse_matrix_elementwise_product_test_double("double");

template <typename DataType_>
class DenseMatrixSparseMatrixElementwiseProductQuickTest :
    public QuickTest
{
    public:
        DenseMatrixSparseMatrixElementwiseProductQuickTest(const std::string & type) :
            QuickTest("dense_matrix_sparse_matrix_elementwise_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size (11);
            DenseMatrix<DataType_> dm1(size, size + 1, DataType_(2)), dm3(size, size + 1, DataType_(0));
            SparseMatrix<DataType_> sm2(size, size + 1, size / 7 + 1);
            for (typename MutableMatrix<DataType_>::ElementIterator j_end(sm2.end_elements()), j(sm2.begin_elements()), 
                k(dm3.begin_elements()) ; j != j_end ; ++j, ++k)
            {
                if (j.index() % 7 == 0) 
                {
                    *j = DataType_(3);
                    *k = DataType_(6);                        
                }
            }  
                        
            DenseMatrix<DataType_> & prod(MatrixElementwiseProduct<DataType_>::value(dm1, sm2));

            TEST_CHECK_EQUAL(prod, dm3);

            DenseMatrix<DataType_> dm03(5, 6);
            SparseMatrix<DataType_> sm01(5, 5, 1), sm02(6, 6, 1), sm03(5, 6, 1);

            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(dm03, sm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(dm03, sm02), MatrixColumnsDoNotMatch); 
        }
};
DenseMatrixSparseMatrixElementwiseProductQuickTest<float> dense_matrix_sparse_matrix_elementwise_product_quick_test_float("float");
DenseMatrixSparseMatrixElementwiseProductQuickTest<double> dense_matrix_sparse_matrix_elementwise_product_quick_test_double("double");

template <typename DataType_>
class SparseMatrixBandedMatrixElementwiseProductTest :
    public BaseTest
{
    public:
        SparseMatrixBandedMatrixElementwiseProductTest(const std::string & type) :
            BaseTest("sparse_matrix_banded_matrix_elementwise_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseVector<DataType_> * dv2 (new DenseVector<DataType_>(size, DataType_(3)));                    
                BandedMatrix<DataType_> bm2(size, dv2); 
                SparseMatrix<DataType_> sm1(size, size, size / 7 + 1), sm3(size, size, size / 7 + 1);
                typename BandedMatrix<DataType_>::ConstElementIterator j(bm2.begin_elements());                   
                for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()),
                    i_end(sm1.end_elements()), k(sm3.begin_elements()) ; i != i_end ; ++i, ++j, ++k)
                {
                    if (i.index() % 7 == 0) 
                    {
                        *i = DataType_(2);
                        if (*j != DataType_(0))
                            *k = DataType_(6);                        
                    }                       
                } 

                SparseMatrix<DataType_> & prod(MatrixElementwiseProduct<DataType_>::value(sm1, bm2));

                TEST_CHECK_EQUAL(prod, sm3);
            }

            BandedMatrix<DataType_> bm01(5);
            SparseMatrix<DataType_> sm02(5, 6, 1), sm03(6, 5, 1);

            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(sm02, bm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(sm03, bm01), MatrixColumnsDoNotMatch);            
        }
};
SparseMatrixBandedMatrixElementwiseProductTest<float> sparse_matrix_banded_matrix_elementwise_product_test_float("float");
SparseMatrixBandedMatrixElementwiseProductTest<double> sparse_matrix_banded_matrix_elementwise_product_test_double("double");

template <typename DataType_>
class SparseMatrixBandedMatrixElementwiseProductQuickTest :
    public QuickTest
{
    public:
        SparseMatrixBandedMatrixElementwiseProductQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_banded_matrix_elementwise_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(11);
            DenseVector<DataType_> * dv2 (new DenseVector<DataType_>(size, DataType_(3)));                    
            BandedMatrix<DataType_> bm2(size, dv2); 
            SparseMatrix<DataType_> sm1(size, size, size / 7 + 1), sm3(size, size, size / 7 + 1);
            typename BandedMatrix<DataType_>::ConstElementIterator j(bm2.begin_elements());                   
            for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()),
                i_end(sm1.end_elements()), k(sm3.begin_elements()) ; i != i_end ; ++i, ++j, ++k)
            {
                if (i.index() % 7 == 0) 
                {
                    *i = DataType_(2);
                    if (*j != DataType_(0))
                        *k = DataType_(6);                        
                } 
            }

            SparseMatrix<DataType_> & prod(MatrixElementwiseProduct<DataType_>::value(sm1, bm2));

            TEST_CHECK_EQUAL(prod, sm3);

            BandedMatrix<DataType_> bm01(5);
            SparseMatrix<DataType_> sm02(5, 6, 1), sm03(6, 5, 1);

            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(sm02, bm01), MatrixRowsDoNotMatch);
            TEST_CHECK_THROWS(MatrixElementwiseProduct<DataType_>::value(sm03, bm01), MatrixColumnsDoNotMatch);            
        }
};
SparseMatrixBandedMatrixElementwiseProductQuickTest<float> sparse_matrix_banded_matrix_elementwise_product_quick_test_float("float");
SparseMatrixBandedMatrixElementwiseProductQuickTest<double> sparse_matrix_banded_matrix_elementwise_product_quick_test_double("double");

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
                    }
                    if (i.index() % 7 == 0) 
                    {
                        *j = DataType_(3);
                    }
                    if (i.index() % 7 == 0 && i.index() % 10 == 0) 
                    {
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
            unsigned long size(11);
            SparseMatrix<DataType_> sm1(size, size + 1, size / 8 + 1), 
                sm2(size, size + 1, size / 7 + 1), sm3(size, size + 1, size / 8 + 1 );
            for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()), 
                i_end(sm1.end_elements()), j(sm2.begin_elements()), k(sm3.begin_elements()) ; 
                i != i_end ; ++i, ++j, ++k)
            {
                if (i.index() % 10 == 0) 
                {
                    *i = DataType_(2);
                }
                if (i.index() % 7 == 0) 
                {
                    *j = DataType_(3);
                }
                if (i.index() % 7 == 0 && i.index() % 10 == 0) 
                {
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
