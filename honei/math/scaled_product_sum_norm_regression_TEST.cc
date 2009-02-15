/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2009 Markus Geveler <apryde@gmx.de>
 *
 * This file is part of the Math C++ library. LibMath is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * LibMath is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <honei/math/scaled_product_sum_norm.hh>
#include <unittest/unittest.hh>
#include <honei/util/stringify.hh>
#include <iostream>

using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DT_>
class SPSNormRegressionTest:
    public BaseTest
{
    private:
        unsigned long _size;
    public:
        SPSNormRegressionTest(const std::string & tag, unsigned long size) :
            BaseTest("ScaledProductSumNorm regression test<" + tag + ">")
    {
        register_tag(Tag_::name);
        this->_size = size;
    }

        virtual void run() const
        {
            //TODO: implement regression test, as soon as backend-impls are
            //available

            DenseVector<DT_> x(_size, DT_(1));
            DenseVector<DT_> y(_size, DT_(1));
            BandedMatrix<DT_> A(_size);

            A.insert_band(0, x.copy());
            A.insert_band(1, x.copy());
            A.insert_band(-1, x.copy());

            DT_ result = ScaledProductSumNorm<Tag_>::value(DT_(0.1), y, DT_(0.1), A, x);
            TEST_CHECK(true);
        }
};

SPSNormRegressionTest<tags::CPU, float> spsnorm_regr_test_float("float", 10000);
SPSNormRegressionTest<tags::CPU, double> spsnorm_regr_test_double("double", 10000);
