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

///< \todo Fix compile errors
#include <libla/dense_matrix.hh>
#include <libla/dense_vector.hh>
#include <libla/sparse_vector.hh>
#include <libla/matrix_vector_product.hh>
#include <libla/matrix_error.cc>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using namespace tests;

template <typename DataType_>
class BandedMatrixDenseVectorProductTest :
    public BaseTest
{
    public:
        BandedMatrixDenseVectorProductTest(const std::string & type) :
            BaseTest("banded_matrix_dense_vector_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseVector<DataType_> * dv1 (new DenseVector<DataType_>(size, DataType_(2)));
                BandedMatrix<DataType_> bm1(size, dv1);
                DenseVector<DataType_> dv4 (size, DataType_(2));
                bm1.band(1)= dv4;
                bm1.band(-2)= dv4;                
                DenseVector<DataType_> dv2(size, DataType_(3)),  dv3(size, DataType_(6 * 3));
                DenseVector<DataType_> prod (MatrixVectorProduct<DataType_>::value(bm1, dv2));

                TEST_CHECK_EQUAL(prod, dv3);
            }

            BandedMatrix<DataType_> bm01(5);
            DenseVector<DataType_> dv01(4, DataType_(1));
            TEST_CHECK_THROWS(MatrixVectorProduct<DataType_>::value(bm01, dv01), MatrixRowsDoNotMatch);
        }
};
BandedMatrixDenseVectorProductTest<float> banded_matrix_dense_vector_product_test_float("float");
BandedMatrixDenseVectorProductTest<double> banded_matrix_dense_vector_product_test_double("double");

template <typename DataType_>
class BandedMatrixDenseVectorProductQuickTest :
    public QuickTest
{
    public:
        BandedMatrixDenseVectorProductQuickTest(const std::string & type) :
            QuickTest("banded_matrix_dense_vector_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseVector<DataType_> * dv1 (new DenseVector<DataType_>(size, DataType_(2)));
            BandedMatrix<DataType_> bm1(size, dv1);
            DenseVector<DataType_> dv4 (size, DataType_(2));
            bm1.band(1)= dv4;     
            bm1.band(-2)= dv4;                   
            DenseVector<DataType_> dv2(size, DataType_(3)),  dv3(size, DataType_(6 * 3));
            DenseVector<DataType_> prod (MatrixVectorProduct<DataType_>::value(bm1, dv2));

            TEST_CHECK_EQUAL(prod, dv3);

            BandedMatrix<DataType_> bm01(5);
            DenseVector<DataType_> dv01(4, DataType_(1));
            TEST_CHECK_THROWS(MatrixVectorProduct<DataType_>::value(bm01, dv01), MatrixRowsDoNotMatch);
        }
};
BandedMatrixDenseVectorProductQuickTest<float> banded_matrix_dense_vector_product_quick_test_float("float");
BandedMatrixDenseVectorProductQuickTest<double> banded_matrix_dense_vector_product_quick_test_double("double");

template <typename DataType_>
class BandedMatrixSparseVectorProductTest :
    public BaseTest
{
    public:
        BandedMatrixSparseVectorProductTest(const std::string & type) :
            BaseTest("banded_matrix_sparse_vector_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseVector<DataType_> * dv1 (new DenseVector<DataType_>(size, DataType_(2)));
                BandedMatrix<DataType_> bm1(size, dv1);
                DenseVector<DataType_> dv4 (size, DataType_(2));
                bm1.band(1)= dv4;      
                bm1.band(-2)= dv4;                              
                SparseVector<DataType_> sv1(size, size / 8 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0) *i = DataType_(3);
                }
                SparseVector<DataType_> sv2(size, size / 8 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv2.begin_elements()), i_end(sv2.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0) *i = DataType_(6 * 3);
                }
                SparseVector<DataType_> prod(MatrixVectorProduct<DataType_>::value(bm1, sv1));

                TEST_CHECK_EQUAL(prod, sv2);
            }

            BandedMatrix<DataType_> bm01(5);
            SparseVector<DataType_> sv01(4, 3);
            TEST_CHECK_THROWS(MatrixVectorProduct<DataType_>::value(bm01, sv01), MatrixRowsDoNotMatch);
        }
};
BandedMatrixSparseVectorProductTest<float> banded_matrix_sparse_vector_product_test_float("float");
BandedMatrixSparseVectorProductTest<double> banded_matrix_sparse_vector_product_test_double("double");

template <typename DataType_>
class BandedMatrixSparseVectorProductQuickTest :
    public QuickTest
{
    public:
        BandedMatrixSparseVectorProductQuickTest(const std::string & type) :
            QuickTest("banded_matrix_sparse_vector_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(20);
            DenseVector<DataType_> * dv1 (new DenseVector<DataType_>(size, DataType_(2)));
            BandedMatrix<DataType_> bm1(size, dv1);
            DenseVector<DataType_> dv4 (size, DataType_(2));
            bm1.band(1)= dv4;   
            bm1.band(-2)= dv4;                                         
            SparseVector<DataType_> sv1(size, size / 8 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 10 == 0) *i = DataType_(3);
            }
            SparseVector<DataType_> sv2(size, size / 8 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv2.begin_elements()), i_end(sv2.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 10 == 0) *i = DataType_(6 * 3);
            }
            SparseVector<DataType_> prod(MatrixVectorProduct<DataType_>::value(bm1, sv1));

            TEST_CHECK_EQUAL(prod, sv2);

            BandedMatrix<DataType_> bm01(5);
            SparseVector<DataType_> sv01(4, 3);
            TEST_CHECK_THROWS(MatrixVectorProduct<DataType_>::value(bm01, sv01), MatrixRowsDoNotMatch);
        }
};
BandedMatrixSparseVectorProductQuickTest<float> banded_matrix_sparse_vector_product_quick_test_float("float");
BandedMatrixSparseVectorProductQuickTest<double> banded_matrix_sparse_vector_product_quick_test_double("double");

template <typename DataType_>
class DenseMatrixDenseVectorProductTest :
    public BaseTest
{
    public:
        DenseMatrixDenseVectorProductTest(const std::string & type) :
            BaseTest("dense_matrix_dense_vector_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size, size + 1, DataType_(2));
                DenseVector<DataType_> dv1(size, DataType_(3)),  dv2(size + 1, DataType_(6 * size));
                DenseVector<DataType_> prod(MatrixVectorProduct<DataType_>::value(dm1, dv1));

                TEST_CHECK_EQUAL(prod, dv2);
            }

            DenseMatrix<DataType_> dm01(3, 4, DataType_(1));
            DenseVector<DataType_> dv01(4, DataType_(1));

            TEST_CHECK_THROWS(MatrixVectorProduct<DataType_>::value(dm01, dv01), MatrixRowsDoNotMatch);
        }
};

DenseMatrixDenseVectorProductTest<float> dense_matrix_dense_vector_product_test_float("float");
DenseMatrixDenseVectorProductTest<double> dense_matrix_dense_vector_product_test_double("double");

template <typename DataType_>
class DenseMatrixDenseVectorProductQuickTest :
    public QuickTest
{
    public:
        DenseMatrixDenseVectorProductQuickTest(const std::string & type) :
            QuickTest("dense_matrix_dense_vector_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseMatrix<DataType_> dm1(size, size + 1, DataType_(2));
            DenseVector<DataType_> dv1(size, DataType_(3)),  dv2(size + 1, DataType_(6 * size));
            DenseVector<DataType_> prod(MatrixVectorProduct<DataType_>::value(dm1, dv1));

            TEST_CHECK_EQUAL(prod, dv2);

            DenseMatrix<DataType_> dm01(3, 4, DataType_(1));
            DenseVector<DataType_> dv01(4, DataType_(1));

            TEST_CHECK_THROWS(MatrixVectorProduct<DataType_>::value(dm01, dv01), MatrixRowsDoNotMatch);

        }
};
DenseMatrixDenseVectorProductQuickTest<float> dense_matrix_dense_vector_product_quick_test_float("float");
DenseMatrixDenseVectorProductQuickTest<double> dense_matrix_dense_vector_product_quick_test_double("double");

template <typename DataType_>
class DenseMatrixSparseVectorProductTest :
    public BaseTest
{
    public:
        DenseMatrixSparseVectorProductTest(const std::string & type) :
            BaseTest("dense_matrix_sparse_vector_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(11) ; size < (1 << 12) ; size <<= 1)
            {
                DenseMatrix<DataType_> dm1(size, size + 1, DataType_(2));
                SparseVector<DataType_> sv1(size, size / 8 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0) *i = DataType_(3);
                }
                SparseVector<DataType_> sv2(size + 1, size / 8 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv2.begin_elements()), i_end(sv2.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = static_cast<DataType_>(6 * (size / 10 + 1));
                }
                DenseVector<DataType_> prod(MatrixVectorProduct<DataType_>::value(dm1, sv1));

                TEST_CHECK_EQUAL(prod, sv2);
            }

            DenseMatrix<DataType_> dm01(3, 4, DataType_(1));
            SparseVector<DataType_> sv01(4, 3);

            TEST_CHECK_THROWS(MatrixVectorProduct<DataType_>::value(dm01, sv01), MatrixRowsDoNotMatch);
        }
};

DenseMatrixSparseVectorProductTest<float> dense_matrix_sparse_vector_product_test_float("float");
DenseMatrixSparseVectorProductTest<double> dense_matrix_sparse_vector_product_test_double("double");

template <typename DataType_>
class DenseMatrixSparseVectorProductQuickTest :
    public QuickTest
{
    public:
        DenseMatrixSparseVectorProductQuickTest(const std::string & type) :
            QuickTest("dense__matrix_sparse_vector_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(20);
            DenseMatrix<DataType_> dm1(size, size + 1, DataType_(2));
            SparseVector<DataType_> sv1(size, size / 8 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 10 == 0) *i = DataType_(3);
            }
            SparseVector<DataType_> sv2(size + 1, size / 8 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv2.begin_elements()), i_end(sv2.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = DataType_(12);
            }
            DenseVector<DataType_> prod(MatrixVectorProduct<DataType_>::value(dm1, sv1));

            TEST_CHECK_EQUAL(prod, sv2);

            DenseMatrix<DataType_> dm01(3, 4, DataType_(1));
            SparseVector<DataType_> sv01(4, 3);

            TEST_CHECK_THROWS(MatrixVectorProduct<DataType_>::value(dm01, sv01), MatrixRowsDoNotMatch);
        }
};

DenseMatrixSparseVectorProductQuickTest<float> dense_matrix_sparse_vector_product_test_quick_float("float");
DenseMatrixSparseVectorProductQuickTest<double> dense_matrix_sparse_vector_product_test_quick_double("double");

template <typename DataType_>
class SparseMatrixDenseVectorProductQuickTest :
    public QuickTest
{
    public:
        SparseMatrixDenseVectorProductQuickTest(const std::string & type) :
            QuickTest("sparse_matrix_dense_vector_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(20);
            SparseMatrix<DataType_> sm1(size, size + 1, size / 8 + 1);
            for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()), 
                i_end(sm1.end_elements()) ; i != i_end ; ++i)
            {
                *i = 2;
            }
            DenseVector<DataType_> dv1(size, DataType_(3)),  dv2(size + 1, DataType_(6 * size));
            SparseVector<DataType_> prod(MatrixVectorProduct<DataType_>::value(sm1, dv1));

            TEST_CHECK_EQUAL(prod, dv2);

            SparseMatrix<DataType_> sm01(3, 4, 1);
            DenseVector<DataType_> dv01(4, DataType_(1));

            TEST_CHECK_THROWS(MatrixVectorProduct<DataType_>::value(sm01, dv01), MatrixRowsDoNotMatch);

        }
};
SparseMatrixDenseVectorProductQuickTest<float> sparse_matrix_dense_vector_product_quick_test_float("float");
SparseMatrixDenseVectorProductQuickTest<double> sparse_matrix_dense_vector_product_quick_test_double("double");

template <typename DataType_>
class SparseMatrixSparseVectorProductQuickTest :
    public QuickTest
{
    public:
        SparseMatrixSparseVectorProductQuickTest(const std::string & type) :
            QuickTest("sparse__matrix_sparse_vector_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(20);
            SparseMatrix<DataType_> sm1(size, size + 1, size / 8 + 1);
            for (typename MutableMatrix<DataType_>::ElementIterator i(sm1.begin_elements()), 
                i_end(sm1.end_elements()) ; i != i_end ; ++i)
            {
                *i = 2;
            }
            SparseVector<DataType_> sv1(size, size / 8 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 10 == 0) *i = DataType_(3);
            }
            SparseVector<DataType_> sv2(size + 1, size / 8 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv2.begin_elements()), i_end(sv2.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = DataType_(12);
            }
            SparseVector<DataType_> prod(MatrixVectorProduct<DataType_>::value(sm1, sv1));

            TEST_CHECK_EQUAL(prod, sv2);

            SparseMatrix<DataType_> sm01(3, 4, 1);
            SparseVector<DataType_> sv01(4, 3);

            TEST_CHECK_THROWS(MatrixVectorProduct<DataType_>::value(sm01, sv01), MatrixRowsDoNotMatch);
        }
};

SparseMatrixSparseVectorProductQuickTest<float> sparse_matrix_sparse_vector_product_test_quick_float("float");
SparseMatrixSparseVectorProductQuickTest<double> sparse_matrix_sparse_vector_product_test_quick_double("double");
