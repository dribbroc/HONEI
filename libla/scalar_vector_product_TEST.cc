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
#include <scalar_vector_product.hh>
#include <libla/vector_norm.hh>
#include <unittest/unittest.hh>

#include <limits>
#include <tr1/memory>
#include <iostream>

using namespace pg512;
using  namespace tests;

template <typename DataType_>
class DenseScalarVectorProductTest :
    public BaseTest
{
    public:
        DenseScalarVectorProductTest(const std::string & type) :
            BaseTest("dense_scalar_vector_product_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            for (unsigned long size(1) ; size < (1 << 14) ; size <<= 1)
            {
                std::tr1::shared_ptr<DenseVector<DataType_> > dv1(new DenseVector<DataType_>(size,
                    static_cast<DataType_>(3)));

                DenseVector<DataType_> prod1(ScalarVectorProduct<DataType_>::value(static_cast<DataType_>(2), *dv1));
                DataType_ v1(VectorNorm<DataType_, vnt_l_one>::value(prod1));
                TEST_CHECK_EQUAL(v1, 6 * size);
            }
        }
};

DenseScalarVectorProductTest<float> dense_scalar_vector_product_test_float("float");
DenseScalarVectorProductTest<double> dense_scalar_vector_product_test_double("double");

template <typename DataType_>
class DenseScalarVectorProductQuickTest :
    public QuickTest
{
    public:
        DenseScalarVectorProductQuickTest(const std::string & type) :
            QuickTest("dense_scalar_vector_product_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            unsigned long size(5);
            std::tr1::shared_ptr<DenseVector<DataType_> > dv1(new DenseVector<DataType_>(size,
                static_cast<DataType_>(3)));

            DenseVector<DataType_> prod1(ScalarVectorProduct<DataType_>::value(static_cast<DataType_>(2), *dv1));
            DataType_ v1(VectorNorm<DataType_, vnt_l_one>::value(prod1));
            TEST_CHECK_EQUAL(v1, 6 * size);
        }
};
DenseScalarVectorProductQuickTest<float>  dense_scalar_vector_product_quick_test_float("float");
DenseScalarVectorProductQuickTest<double> dense_scalar_vector_product_quick_test_double("double");
