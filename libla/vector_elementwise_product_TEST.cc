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

#include <libla/dense_vector.hh>
#include <libla/sparse_vector.hh>
#include <libla/vector_elementwise_product.hh>
#include <libla/vector_error.cc>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using namespace tests;

template <typename DataType_>
class DenseVectorElementwiseProductTest :
    public BaseTest
{
    public:
        DenseVectorElementwiseProductTest(const std::string & type) :
            BaseTest("dense_vector_elementwise_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size, DataType_(2)), dv2(size, DataType_(3)),
                    dv3(size, DataType_(6));
                DenseVector<DataType_> prod(VectorElementwiseProduct<DataType_>::value(dv1, dv2));

                TEST_CHECK_EQUAL(prod, dv3);
            }

            DenseVector<DataType_> dv01(3, DataType_(1)), dv02(4, DataType_(1));

            TEST_CHECK_THROWS(VectorElementwiseProduct<DataType_>::value(dv02, dv01), VectorSizeDoesNotMatch);
        }
};

DenseVectorElementwiseProductTest<float> dense_vector_elementwise_product_test_float("float");
DenseVectorElementwiseProductTest<double> dense_vector_elementwise_product_test_double("double");

template <typename DataType_>
class DenseVectorElementwiseProductQuickTest :
    public QuickTest
{
    public:
        DenseVectorElementwiseProductQuickTest(const std::string & type) :
            QuickTest("dense_vector_elementwise_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseVector<DataType_> dv1(size, DataType_(2)), dv2(size, DataType_(3)),
                dv3(size, DataType_(6));
            DenseVector<DataType_> prod(VectorElementwiseProduct<DataType_>::value(dv1, dv2));

            TEST_CHECK_EQUAL(prod, dv3);


            DenseVector<DataType_> dv01(3, DataType_(1)), dv02(4, DataType_(1));

            TEST_CHECK_THROWS(VectorElementwiseProduct<DataType_>::value(dv02, dv01), VectorSizeDoesNotMatch);

        }
};
DenseVectorElementwiseProductQuickTest<float> dense_vector_elementwise_product_quick_test_float("float");
DenseVectorElementwiseProductQuickTest<double> dense_vector_elementwise_product_quick_test_double("double");

template <typename DataType_>
class SparseVectorDenseVectorElementwiseProductTest :
    public BaseTest
{
    public:
        SparseVectorDenseVectorElementwiseProductTest(const std::string & type) :
            BaseTest("sparse_vector_dense_vector_elementwise_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                SparseVector<DataType_> sv1(size, size / 7 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 7 == 0) *i = static_cast<DataType_>(2);
                }            
                DenseVector<DataType_> dv2(size, DataType_(3));
                SparseVector<DataType_> sv3(size, size / 8 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv3.begin_elements()), i_end(sv3.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 7 == 0) 
                        *i = static_cast<DataType_>(6);
                }  
                SparseVector<DataType_> prod(VectorElementwiseProduct<DataType_>::value(sv1, dv2));

                TEST_CHECK_EQUAL(prod, sv3);
            }

            SparseVector<DataType_> sv01(3, 1);
            DenseVector<DataType_> dv02(4);

            TEST_CHECK_THROWS(VectorElementwiseProduct<DataType_>::value(sv01, dv02), VectorSizeDoesNotMatch);
        }
};
SparseVectorDenseVectorElementwiseProductTest<float> sparse_vector_dense_vector_elementwise_product_test_float("float");
SparseVectorDenseVectorElementwiseProductTest<double> sparse_vector_dense_vector_elementwise_product_test_double("double");

template <typename DataType_>
class SparseVectorDenseVectorElementwiseProductQuickTest :
    public QuickTest
{
    public:
        SparseVectorDenseVectorElementwiseProductQuickTest(const std::string & type) :
            QuickTest("sparse_vector_dense_vector_elementwise_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(22);
            SparseVector<DataType_> sv1(size, size / 7 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 7 == 0) *i = static_cast<DataType_>(2);
            }            
            DenseVector<DataType_> dv2(size, DataType_(3));
            SparseVector<DataType_> sv3(size, size / 8 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv3.begin_elements()), i_end(sv3.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 7 == 0) 
                    *i = static_cast<DataType_>(6);
            }  
            SparseVector<DataType_> prod(VectorElementwiseProduct<DataType_>::value(sv1, dv2));

            TEST_CHECK_EQUAL(prod, sv3);

            SparseVector<DataType_> sv01(3, 1);
            DenseVector<DataType_> dv02(4);

            TEST_CHECK_THROWS(VectorElementwiseProduct<DataType_>::value(sv01, dv02), VectorSizeDoesNotMatch);
        }
};
SparseVectorDenseVectorElementwiseProductQuickTest<float> sparse_vector_dense_vector_elementwise_product_quick_test_float("float");
SparseVectorDenseVectorElementwiseProductQuickTest<double> sparse_vector_dense_vector_elementwise_product_quick_test_double("double");

template <typename DataType_>
class SparseVectorElementwiseProductTest :
    public BaseTest
{
    public:
        SparseVectorElementwiseProductTest(const std::string & type) :
            BaseTest("sparse_vector_elementwise_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(10) ; size < (1 << 12) ; size <<= 1)
            {
                SparseVector<DataType_> sv1(size, size / 7 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0) *i = static_cast<DataType_>(2);
                }            
                SparseVector<DataType_> sv2(size, size / 8 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv2.begin_elements()), i_end(sv2.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 7 == 0) *i = static_cast<DataType_>(3);
                }  
                SparseVector<DataType_> sv3(size, size / 8 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv3.begin_elements()), i_end(sv3.end_elements()) ;
                        i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0 && i.index() % 7 == 0) 
                        *i = static_cast<DataType_>(6);
                }  
                SparseVector<DataType_> prod(VectorElementwiseProduct<DataType_>::value(sv1, sv2));

                TEST_CHECK_EQUAL(prod, sv3);
            }

            SparseVector<DataType_> sv01(3, 1), sv02(4, 1);

            TEST_CHECK_THROWS(VectorElementwiseProduct<DataType_>::value(sv02, sv01), VectorSizeDoesNotMatch);
        }
};
SparseVectorElementwiseProductTest<float> sparse_vector_elementwise_product_test_float("float");
SparseVectorElementwiseProductTest<double> sparse_vector_elementwise_product_test_double("double");

template <typename DataType_>
class SparseVectorElementwiseProductQuickTest :
    public QuickTest
{
    public:
        SparseVectorElementwiseProductQuickTest(const std::string & type) :
            QuickTest("sparse_vector_elementwise_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(100);
            SparseVector<DataType_> sv1(size, size / 7 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 10 == 0) *i = static_cast<DataType_>(2);
            }            
            SparseVector<DataType_> sv2(size, size / 8 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv2.begin_elements()), i_end(sv2.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 7 == 0) *i = static_cast<DataType_>(3);
            }  
            SparseVector<DataType_> sv3(size, size / 8 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv3.begin_elements()), i_end(sv3.end_elements()) ;
                    i != i_end ; ++i)
            {
                if (i.index() % 10 == 0 && i.index() % 7 == 0) 
                    *i = static_cast<DataType_>(6);
            }  
            SparseVector<DataType_> prod(VectorElementwiseProduct<DataType_>::value(sv1, sv2));

            TEST_CHECK_EQUAL(prod, sv3);

            SparseVector<DataType_> sv01(3, 1), sv02(4, 1);

            TEST_CHECK_THROWS(VectorElementwiseProduct<DataType_>::value(sv02, sv01), VectorSizeDoesNotMatch);
        }
};
SparseVectorElementwiseProductQuickTest<float> sparse_vector_elementwise_product_quick_test_float("float");
SparseVectorElementwiseProductQuickTest<double> sparse_vector_elementwise_product_quick_test_double("double");
