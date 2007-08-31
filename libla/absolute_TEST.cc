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
#include <libla/absolute.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace honei;
using namespace tests;

template <typename DataType_>
class DenseVectorAbsoluteValueTest :
    public BaseTest
{
    public:
        DenseVectorAbsoluteValueTest(const std::string & type) :
            BaseTest("dense_vector_absolute_value_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                DenseVector<DataType_> dv1(size), dv2(size);
                for (typename Vector<DataType_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()),
                        j(dv2.begin_elements()) ; i != i_end ; ++i, ++j)
                {
                    *i = DataType_(i.index() / 0.987654321 * (i.index() % 2 == 1 ? -1 : 1));
                    *j = DataType_(i.index() / 0.987654321);
                }

                TEST_CHECK_EQUAL(Absolute<>::value(dv1), dv2);
            }
        }
};

DenseVectorAbsoluteValueTest<float> dense_vector_absolute_value_test_float("float");
DenseVectorAbsoluteValueTest<double> dense_vector_absolute_value_test_double("double");

template <typename DataType_>
class DenseVectorAbsoluteQuickTest :
    public QuickTest
{
    public:
        DenseVectorAbsoluteQuickTest(const std::string & type) :
            QuickTest("dense_vector_absolute_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            DenseVector<DataType_> dv1(size), dv2(size);
            for (typename Vector<DataType_>::ElementIterator i(dv1.begin_elements()), i_end(dv1.end_elements()),
                    j(dv2.begin_elements()) ; i != i_end ; ++i, ++j)
            {
                *i = DataType_(i.index() / 0.987654321 * (i.index() % 2 == 1 ? -1 : 1));
                *j = DataType_(i.index() / 0.987654321);
            }

            TEST_CHECK_EQUAL(Absolute<>::value(dv1), dv2);
        }
};
DenseVectorAbsoluteQuickTest<float>  dense_vector_absolute_quick_test_float("float");
DenseVectorAbsoluteQuickTest<double> dense_vector_absolute_quick_test_double("double");

template <typename DataType_>
class SparseVectorAbsoluteValueTest :
    public BaseTest
{
    public:
        SparseVectorAbsoluteValueTest(const std::string & type) :
            BaseTest("sparse_vector_absolute_value_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                SparseVector<DataType_> sv1(size, size / 8 + 1), sv2(size, size / 9 + 1);
                for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()),
                        j(sv2.begin_elements()) ; i != i_end ; ++i, ++j)
                {
                    if (i.index() % 3 == 0)
                    {
                        *i = DataType_(i.index() / 0.987654321 * (i.index() % 2 == 1 ? -1 : 1));
                        *j = DataType_(i.index() / 0.987654321);
                    }
                }

                TEST_CHECK_EQUAL(Absolute<>::value(sv1), sv2);
            }
        }
};

SparseVectorAbsoluteValueTest<float> sparse_vector_absolute_value_test_float("float");
SparseVectorAbsoluteValueTest<double> sparse_vector_absolute_value_test_double("double");

template <typename DataType_>
class SparseVectorAbsoluteValueQuickTest :
    public QuickTest
{
    public:
        SparseVectorAbsoluteValueQuickTest(const std::string & type) :
            QuickTest("sparse_vector_absolute_value_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(20);
            SparseVector<DataType_> sv1(size, size / 8 + 1), sv2(size, size / 9 + 1);
            for (typename Vector<DataType_>::ElementIterator i(sv1.begin_elements()), i_end(sv1.end_elements()),
                    j(sv2.begin_elements()) ; i != i_end ; ++i, ++j)
            {
                if (i.index() % 3 == 0)
                {
                    *i = DataType_(i.index() / 0.987654321 * (i.index() % 2 == 1 ? -1 : 1));
                    *j = DataType_(i.index() / 0.987654321);
                }
            }

            TEST_CHECK_EQUAL(Absolute<>::value(sv1), sv2);
        }
};

SparseVectorAbsoluteValueQuickTest<float> sparse_vector_absolute_value_quick_test_float("float");
SparseVectorAbsoluteValueQuickTest<double> sparse_vector_absolute_value_quick_test_double("double");
