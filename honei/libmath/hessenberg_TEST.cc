/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Danny van Dyk <danny.dyk@tu-dortmund.de>
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

#include <honei/libmath/hessenberg.hh>
#include <unittest/unittest.hh>
#include <honei/util/stringify.hh>

#include <cmath>
#include <limits>

using namespace honei;
using namespace tests;

template <typename Tag_, typename DT_>
class HessenbergTest :
    public QuickTest
{
    private:

    public:
        HessenbergTest(const std::string & type) :
            QuickTest("hessenberg_test<" + type + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            for (unsigned size(3) ; size < 100 ; ++size)
            {
                DenseMatrix<DT_> dm(size, size, DT_(1));

                Hessenberg<Tag_>::value(dm);

                TEST_CHECK_EQUAL_WITHIN_EPS(dm(0, 0), DT_(1.0), DT_(size) * std::numeric_limits<DT_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS(dm(1, 0), -sqrt(DT_(size) - DT_(1.0)), DT_(size * size) * std::numeric_limits<DT_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS(dm(0, 1), -sqrt(DT_(size) - DT_(1.0)), DT_(size * size) * std::numeric_limits<DT_>::epsilon());
                TEST_CHECK_EQUAL_WITHIN_EPS(dm(1, 1), DT_(size) - DT_(1.0), DT_(size * size) * std::numeric_limits<DT_>::epsilon());

                for (unsigned i(2) ; i < size ; ++i)
                {
                    TEST_CHECK_EQUAL_WITHIN_EPS(dm(i, i), DT_(0.0), DT_(size * size) * std::numeric_limits<DT_>::epsilon());
                }
            }
        }
};

HessenbergTest<tags::CPU, float> hessenberg_test_float("float");

