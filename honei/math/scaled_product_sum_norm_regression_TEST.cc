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
            DenseVector<DT_> x(_size, DT_(1.12345));
            DenseVector<DT_> y(_size, DT_(0.23456));

            DenseVector<DT_> non_diag(x.copy());
            DenseVector<DT_> diag(y.copy());

            BandedMatrixQ1<DT_> A(_size, non_diag, non_diag, non_diag, non_diag, diag, non_diag, non_diag, non_diag, non_diag);


            DT_ result_cpu = ScaledProductSumNorm_TUTORIAL<tags::CPU>::value(DT_(1.234), y, DT_(0.1234), A, x);
            DT_ result_sse = ScaledProductSumNorm_TUTORIAL<Tag_>::value(DT_(1.234), y, DT_(0.1234), A, x);

            std::cout << "result_cpu: " << result_cpu << std::endl;
            std::cout << "result_" << Tag_::name <<": " << result_sse << std::endl;
            TEST_CHECK_EQUAL_WITHIN_EPS(result_cpu, result_sse, std::numeric_limits<float>::epsilon() *_size * 3000);
        }
};

#ifdef HONEI_SSE
SPSNormRegressionTest<tags::CPU::SSE, float> spsnorm_regr_test_sse_float("float", 10000);
SPSNormRegressionTest<tags::CPU::SSE, double> spsnorm_regr_test_sse_double("double", 10000);
#endif
