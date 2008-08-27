/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2008 Markus Geveler <apryde@gmx.de>
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

#include <honei/math/prolongation.hh>
#include <unittest/unittest.hh>
#include <honei/util/stringify.hh>

using namespace honei;
using namespace tests;
using namespace std;

template <typename Tag_, typename DT1_>
class ProlongationTest:
    public BaseTest
{
    public:
        ProlongationTest(const std::string & tag) :
            BaseTest("Prolongate test <" + tag + ">")
        {
            register_tag(Tag_::name);
        }

        virtual void run() const
        {
            DenseVector<DT1_> fine(9, DT1_(0));
            DenseVector<DT1_> coarse(4, DT1_(0));
            DenseVector<unsigned long> mask(8);

            Prolongation<Tag_>::value(fine, coarse, mask);
            TEST_CHECK(true);
        }
};
ProlongationTest<tags::CPU, float> prolongate_test("float");
