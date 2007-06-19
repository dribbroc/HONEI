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
#include <libla/vector_element_sum.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <iostream>
#include <vector>
#include <tr1/memory>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
class DenseVectorElementSumTest :
    public BaseTest
{
    public:
        DenseVectorElementSumTest(const std::string & type) :
            BaseTest("dense_vector_element_sum_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                std::tr1::shared_ptr<DenseVector<DataType_> > dv(new DenseVector<DataType_>(size));
                for (typename Vector<DataType_>::ElementIterator i(dv->begin_elements()), i_end(dv->end_elements()) ;
                        i != i_end ; ++i)
                {
                    *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
                }

                DataType_ v1(VectorElementSum<DataType_>::value(*dv));
                DataType_ s1(size * (size + 1) / 2 / 1.23456789);
                // Behavious similar to size^2 * eps
                DataType_ eps1(s1 * 10 * std::numeric_limits<DataType_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS(v1, s1, eps1);
            }
        }
};

DenseVectorElementSumTest<float> dense_vector_element_sum_test_float("float");
DenseVectorElementSumTest<double> dense_vector_element_sum_test_double("double");

template <typename DataType_>
class DenseVectorElementSumQuickTest :
    public QuickTest
{
    public:
        DenseVectorElementSumQuickTest(const std::string & type) :
            QuickTest("dense_vector_element_sum_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            std::tr1::shared_ptr<DenseVector<DataType_> > dv(new DenseVector<DataType_>(size));
            for (typename Vector<DataType_>::ElementIterator i(dv->begin_elements()), i_end(dv->end_elements()) ;
                    i != i_end ; ++i)
            {
                *i = static_cast<DataType_>((i.index() + 1) / 1.23456789);
            }

            DataType_ v1(VectorElementSum<DataType_>::value(*dv));
            DataType_ s1(size * (size + 1) / 2 / 1.23456789);
            // Behavious similar to size^2 * eps
            DataType_ eps1(s1 * 10 * std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(v1, s1, eps1);
        }
};
DenseVectorElementSumQuickTest<float>  dense_vector_element_sum_quick_test_float("float");
DenseVectorElementSumQuickTest<double> dense_vector_element_sum_quick_test_double("double");
