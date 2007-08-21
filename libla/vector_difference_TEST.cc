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
#include <libla/vector_difference.hh>
#include <libla/vector_norm.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
class DenseVectorDifferenceTest :
    public BaseTest
{
    public:
        DenseVectorDifferenceTest(const std::string & type) :
            BaseTest("dense_vector_difference_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size);
                DenseVector<DataType_> dv2(size);

                for (typename Vector<DataType_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
                }
                for (typename Vector<DataType_>::ElementIterator i(dv2.begin_elements()), i_end(dv2.end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
                }
                DenseVector<DataType_> difference1(VectorDifference<>::value(dv1, dv2));
                DataType_ v1(VectorNorm<DataType_, vnt_l_one>::value(difference1));
                TEST_CHECK_EQUAL(v1, 0);
            }

            DenseVector<DataType_> dv00(1, DataType_(1));
            DenseVector<DataType_> dv01(5, DataType_(1));
            TEST_CHECK_THROWS(VectorDifference<DataType_>::value(dv00, dv01), VectorSizeDoesNotMatch);
        }
};

DenseVectorDifferenceTest<float> dense_vector_difference_test_float("float");
DenseVectorDifferenceTest<double> dense_vector_difference_test_double("double");

template <typename DataType_>
class DenseVectorDifferenceQuickTest :
    public QuickTest
{
    public:
        DenseVectorDifferenceQuickTest(const std::string & type) :
            QuickTest("dense_vector_difference_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseVector<DataType_> dv1(size);
            DenseVector<DataType_> dv2(size);

            for (typename Vector<DataType_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
            }
            for (typename Vector<DataType_>::ElementIterator i(dv2.begin_elements()), i_end(dv2.end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
            }
            DenseVector<DataType_> difference1(VectorDifference<>::value(dv1, dv2));
            DataType_ v1(VectorNorm<DataType_, vnt_l_one>::value(difference1));
            TEST_CHECK_EQUAL(v1, 0);

            DenseVector<DataType_> dv00(1, DataType_(1));
            DenseVector<DataType_> dv01(5, DataType_(1));

            TEST_CHECK_THROWS(VectorDifference<DataType_>::value(dv00, dv01), VectorSizeDoesNotMatch);
        }
};
DenseVectorDifferenceQuickTest<float>  dense_vector_difference_quick_test_float("float");
DenseVectorDifferenceQuickTest<double> dense_vector_difference_quick_test_double("double");

template <typename DataType_>
class SparseVectorDifferenceTest :
    public BaseTest
{
    public:
        SparseVectorDifferenceTest(const std::string & type) :
            BaseTest("sparse_vector_difference_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                SparseVector<DataType_> sv1(size, size / 8 + 1), sv2(size, size / 7 + 1), sv3(size, size / 8 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()),
                    j(sv2.begin_elements()), k(sv3.begin_elements()) ; i != i_end ; ++i, ++j, ++k)
                {
                    if (i.index() % 10 == 0) 
                    {
                        *i = static_cast<DataType_>(((i.index() +1) * 2) / 1.23456789);
                        *k = static_cast<DataType_>(((i.index() +1) * 2) / 1.23456789);
                    }
                    if (i.index() % 7 == 0) 
                    {
                        *j = static_cast<DataType_>((i.index() +1) / 1.23456789);
                        *k = static_cast<DataType_>((i.index() +1) / -1.23456789);
                    }
                    if (i.index() % 7 == 0 && i.index() % 10 == 0)
                    {
                        *k = static_cast<DataType_>((i.index() +1) / 1.23456789);
                    }                                
                }            

            SparseVector<DataType_> difference1(VectorDifference<>::value(sv1, sv2));
            TEST_CHECK_EQUAL(difference1, sv3);
            }

            SparseVector<DataType_> sv00(1, 1), sv01(5, 1);
            TEST_CHECK_THROWS(VectorDifference<DataType_>::value(sv00, sv01), VectorSizeDoesNotMatch);
        }
};

SparseVectorDifferenceTest<float> sparse_vector_difference_test_float("float");
SparseVectorDifferenceTest<double> sparse_vector_difference_test_double("double");

template <typename DataType_>
class SparseVectorDifferenceQuickTest :
    public QuickTest
{
    public:
        SparseVectorDifferenceQuickTest(const std::string & type) :
            QuickTest("sparse_vector_difference_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(25);
            SparseVector<DataType_> sv1(size, size / 8 + 1), sv2(size, size / 7 + 1), sv3(size, size / 8 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()),
                j(sv2.begin_elements()), k(sv3.begin_elements()) ; i != i_end ; ++i, ++j, ++k)
            {
                if (i.index() % 10 == 0) 
                {
                    *i = static_cast<DataType_>(((i.index() +1) * 2) / 1.23456789);
                    *k = static_cast<DataType_>(((i.index() +1) * 2) / 1.23456789);
                }
                if (i.index() % 7 == 0) 
                {
                    *i = static_cast<DataType_>((i.index() +1) / 1.23456789);
                    *k = static_cast<DataType_>((i.index() +1) / -1.23456789);
                }
                if (i.index() % 7 == 0 && i.index() % 10 == 0)
                {
                    *k = static_cast<DataType_>((i.index() +1) / 1.23456789);
                }                                
            }            

            SparseVector<DataType_> difference1(VectorDifference<>::value(sv1, sv2));
            TEST_CHECK_EQUAL(difference1, sv3);

            SparseVector<DataType_> sv00(1, 1), sv01(5, 1);
            TEST_CHECK_THROWS(VectorDifference<DataType_>::value(sv00, sv01), VectorSizeDoesNotMatch);
        }
};

SparseVectorDifferenceQuickTest<float> sparse_vector_difference_quick_test_float("float");
SparseVectorDifferenceQuickTest<double> sparse_vector_difference_quick_test_double("double");
