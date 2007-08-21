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
 
#include <libla/dense_vector.hh>
#include <libla/sparse_vector.hh>
#include <libla/vector_norm.hh>
#include <libla/scalar_product.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using namespace tests;

template <typename DataType_>
class DenseScalarProductTest :
    public BaseTest
{
    public:
        DenseScalarProductTest(const std::string & type) :
            BaseTest("dense_scalar_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                DenseVector<DataType_> dv0 (size, DataType_(0)), dv1(size, DataType_(1));

                DataType_ p0(ScalarProduct<>::value(dv1, dv0));
                DataType_ p1(ScalarProduct<>::value(dv1, dv1));
                TEST_CHECK_EQUAL(p0, 0);
                TEST_CHECK_EQUAL(p1, size);

                DenseVector<DataType_> dv2(size);
                for (typename Vector<DataType_>::ElementIterator i(dv2.begin_elements()), i_end(dv2.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
                }

                DataType_ v2(VectorNorm<DataType_, vnt_l_two, false>::value(dv2));
                DataType_ p2(ScalarProduct<>::value(dv2, dv2));
                TEST_CHECK_EQUAL_WITHIN_EPS(v2, p2, sqrt(std::numeric_limits<DataType_>::epsilon()));
            }

            DenseVector<DataType_> dv00(1, DataType_(1));
            DenseVector<DataType_> dv01(2, DataType_(1));

            TEST_CHECK_THROWS(ScalarProduct<>::value(dv00, dv01), VectorSizeDoesNotMatch);
        }
};

DenseScalarProductTest<float> dense_scalar_product_test_float("float");
DenseScalarProductTest<double> dense_scalar_product_test_double("double");

template <typename DataType_>
class DenseScalarProductQuickTest :
    public QuickTest
{
    public:
        DenseScalarProductQuickTest(const std::string & type) :
            QuickTest("dense_scalar_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(22);
            DenseVector<DataType_> dv0 (size, DataType_(0)), dv1(size, DataType_(1));

            DataType_ p0(ScalarProduct<>::value(dv1, dv0));
            DataType_ p1(ScalarProduct<>::value(dv1, dv1));
            TEST_CHECK_EQUAL(p0, 0);
            TEST_CHECK_EQUAL(p1, size);

            DenseVector<DataType_> dv2(size);
            for (typename Vector<DataType_>::ElementIterator i(dv2.begin_elements()), i_end(dv2.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
            }

            DataType_ v2(VectorNorm<DataType_, vnt_l_two, false>::value(dv2));
            DataType_ p2(ScalarProduct<>::value(dv2, dv2));
            TEST_CHECK_EQUAL_WITHIN_EPS(v2, p2, sqrt(std::numeric_limits<DataType_>::epsilon()));

            DenseVector<DataType_> dv00(1, DataType_(1));
            DenseVector<DataType_> dv01(2, DataType_(1));

            TEST_CHECK_THROWS(ScalarProduct<>::value(dv00, dv01), VectorSizeDoesNotMatch);
        }
};
DenseScalarProductQuickTest<float>  dense_scalar_product_quick_test_float("float");
DenseScalarProductQuickTest<double> dense_scalar_product_quick_test_double("double");

template <typename DataType_>
class SparseScalarProductTest :
    public BaseTest
{
    public:
        SparseScalarProductTest(const std::string & type) :
            BaseTest("sparse_scalar_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                DataType_ p1(0);
                SparseVector<DataType_> sv1(size, size / 7 + 1), sv2(size, size / 8 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()), 
                        j(sv2.begin_elements()) ; i != i_end ; ++i)
                {
                    if (i.index() % 10 == 0) 
                    {
                        *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
                    }
                    if (i.index() % 7 == 0) 
                    {
                        *j = static_cast<DataType_>(((i.index() + 1) * 2) / 1.23456789);
                    }
                    if (i.index() % 7 == 0 && i.index() % 10 == 0) 
                    {
                        p1 += *i * *j;
                    }                    
                }

                DataType_ p0(ScalarProduct<>::value(sv1, sv2));
                TEST_CHECK_EQUAL(p0, p1);
             }

            SparseVector<DataType_> sv00(1, 1);
            SparseVector<DataType_> sv01(2, 1);

            TEST_CHECK_THROWS(ScalarProduct<>::value(sv00, sv01), VectorSizeDoesNotMatch);
        }
};
SparseScalarProductTest<float> sparse_scalar_product_test_float("float");
SparseScalarProductTest<double> sparse_scalar_product_test_double("double");

template <typename DataType_>
class SparseScalarProductQuickTest :
    public QuickTest
{
    public:
        SparseScalarProductQuickTest(const std::string & type) :
            QuickTest("sparse_scalar_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size (22);
            DataType_ p1(0);
            SparseVector<DataType_> sv1(size, size / 7 + 1), sv2(size, size / 8 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()), 
                    j(sv2.begin_elements()) ; i != i_end ; ++i)
            {
                if (i.index() % 10 == 0) 
                {
                    *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
                }
                if (i.index() % 7 == 0) 
                {
                    *j = static_cast<DataType_>(((i.index() + 1) * 2) / 1.23456789);
                }
                if (i.index() % 7 == 0 && i.index() % 10 == 0) 
                {
                    p1 += *i * *j;
                }                    
            }

                DataType_ p0(ScalarProduct<>::value(sv1, sv2));
                TEST_CHECK_EQUAL(p0, p1);

            SparseVector<DataType_> sv00(1, 1);
            SparseVector<DataType_> sv01(2, 1);

            TEST_CHECK_THROWS(ScalarProduct<>::value(sv00, sv01), VectorSizeDoesNotMatch);
        }
};
SparseScalarProductQuickTest<float> sparse_scalar_product_quick_test_float("float");
SparseScalarProductQuickTest<double> sparse_scalar_product_quick_double("double");
