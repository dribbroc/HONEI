/* vim: set number sw=4 sts=4 et nofoldenable : */

/*
 * Copyright (c) 2009 Dirk Ribbrock <dirk.ribbrock@uni-dortmund.de>
 *
 * This file is part of the HONEI C++ library. HONEI is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * HONEI is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */
#include <honei/lbm/tags.hh>
#include <honei/util/unittest.hh>
#include <honei/lbm/lbm_limiter.hh>

using namespace honei;
using namespace tests;
using namespace std;
using namespace lbm;

template <typename DataType_>
class MinModLimiterQuickTest :
    public QuickTest
{
    public:
        MinModLimiterQuickTest(const std::string & type) :
            QuickTest("mm_limiter_quick_test<" + type + ">")
        {
        }

        virtual void run() const
        {
            DataType_ s1(-1.23456);
            DataType_ s2(1.23456);
            DataType_ s3(0.23456);
            TEST_CHECK_EQUAL_WITHIN_EPS(MinModLimiter<tags::CPU>::value(s1), DataType_(0), std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(MinModLimiter<tags::CPU>::value(s2), DataType_(1), std::numeric_limits<DataType_>::epsilon());
            TEST_CHECK_EQUAL_WITHIN_EPS(MinModLimiter<tags::CPU>::value(s3), DataType_(0.23456), std::numeric_limits<DataType_>::epsilon());
        }
};
MinModLimiterQuickTest<float> mm_limiter_quick_test_float("float");
MinModLimiterQuickTest<double> mm_limiter_quick_test_double("double");
